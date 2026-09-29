"""Command-line entry point. Patient-level input/output stays on the local machine."""

import argparse
import hashlib
import json
import os
from pathlib import Path
import sys
import tempfile

import pandas as pd

from . import QCConfig, run_qc


def _sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def _write_atomic(frame: pd.DataFrame, path: Path) -> None:
    with tempfile.NamedTemporaryFile(
        mode="w", encoding="utf-8", newline="", dir=path.parent,
        prefix=f".{path.name}.", delete=False,
    ) as handle:
        temporary = Path(handle.name)
        try:
            frame.to_csv(handle, index=False)
        except BaseException:
            temporary.unlink(missing_ok=True)
            raise
    os.replace(temporary, path)


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description="Cohort-aware QC for PLINK BED/BIM/FAM data")
    parser.add_argument("--plink-prefix", type=Path, required=True)
    parser.add_argument("--groups", type=Path, required=True, help="CSV with fid,iid,group")
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--case-label", default="case")
    parser.add_argument("--control-label", default="control")
    parser.add_argument("--batch-size", type=int, default=1024)
    parser.add_argument("--max-variant-missing", type=float, default=0.05)
    parser.add_argument("--min-sample-call-rate", type=float, default=0.95)
    parser.add_argument("--min-minor-carriers", type=int, default=2)
    parser.add_argument("--hwe-alpha", type=float, default=0.05)
    parser.add_argument("--exclude-indels", action="store_true")
    parser.add_argument(
        "--exclude-control-maf-band", nargs=2, type=float, metavar=("LOW", "HIGH"),
        help="Optional study-specific exclusion; not a general QC default",
    )
    args = parser.parse_args(argv)
    if args.case_label == args.control_label:
        parser.error("case and control labels must differ")
    config = QCConfig(
        max_variant_missing_rate=args.max_variant_missing,
        min_sample_call_rate=args.min_sample_call_rate,
        min_minor_carriers=args.min_minor_carriers,
        hwe_familywise_alpha=args.hwe_alpha,
        exclude_indels=args.exclude_indels,
        control_maf_exclusion_band=(
            tuple(args.exclude_control_maf_band)
            if args.exclude_control_maf_band is not None else None
        ),
    )
    prefix = args.plink_prefix
    paths = {suffix: Path(f"{prefix}.{suffix}") for suffix in ("bed", "bim", "fam")}
    for path in [*paths.values(), args.groups]:
        if not path.is_file():
            parser.error(f"Missing input: {path}")
    try:
        from pandas_plink import read_plink
    except ImportError as exc:
        parser.error("Install the optional PLINK extra: pip install '.[plink]'")
        raise AssertionError("unreachable") from exc

    bim, fam, bed = read_plink(str(paths["bed"]), verbose=False)
    groups = pd.read_csv(args.groups, dtype={"fid": str, "iid": str})
    result = run_qc(
        bed, bim, fam, groups, config,
        case_label=args.case_label, control_label=args.control_label,
        batch_size=args.batch_size,
    )
    args.output.mkdir(parents=True, exist_ok=True)
    _write_atomic(result.variants, args.output / "variants.csv")
    _write_atomic(result.samples, args.output / "samples.csv")
    metadata = {
        "config": config.as_dict(),
        "case_label": args.case_label,
        "control_label": args.control_label,
        "input_sha256": {
            **{suffix: _sha256(path) for suffix, path in paths.items()},
            "groups": _sha256(args.groups),
        },
        "python_version": sys.version.split()[0],
        "n_variants": len(result.variants),
        "n_samples": len(result.samples),
        "n_preliminary_variants": result.preliminary_variant_count,
        "n_final_variants": int(result.variants["pass_qc"].sum()),
        "n_retained_samples": int(result.samples["pass_qc"].sum()),
        "control_hwe_bonferroni_threshold": result.hwe_threshold,
    }
    output = args.output / "run_metadata.json"
    with tempfile.NamedTemporaryFile(
        mode="w", encoding="utf-8", dir=args.output, prefix=".metadata.", delete=False
    ) as handle:
        temp = Path(handle.name)
        json.dump(metadata, handle, indent=2, sort_keys=True)
        handle.write("\n")
    os.replace(temp, output)
    print(f"QC complete: {metadata['n_final_variants']} variants, "
          f"{metadata['n_retained_samples']} samples retained; output: {args.output}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

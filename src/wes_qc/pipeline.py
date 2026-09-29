"""Chunked, two-pass QC without loading the entire PLINK BED into memory."""

from dataclasses import dataclass
from typing import Any

import numpy as np
import pandas as pd

from .core import QCConfig, mark_variant_pass, variant_statistics


@dataclass
class QCResult:
    variants: pd.DataFrame
    samples: pd.DataFrame
    preliminary_variant_count: int
    hwe_threshold: float


def _batch(array: Any, start: int, end: int) -> np.ndarray:
    block = array[start:end, :]
    if hasattr(block, "compute"):
        block = block.compute()
    return np.asarray(block, dtype=float)


def _metadata(bim: pd.DataFrame, n_variants: int) -> pd.DataFrame:
    columns = ["chrom", "pos", "snp", "a0", "a1"]
    if not set(columns).issubset(bim) or len(bim) != n_variants:
        raise ValueError("BIM fields/rows do not match the genotype matrix")
    result = bim[columns].reset_index(drop=True).copy()
    if result[columns].isna().any().any():
        raise ValueError("BIM contains missing variant identifiers or alleles")
    return result


def _groups(fam: pd.DataFrame, groups: pd.DataFrame, n_samples: int) -> pd.DataFrame:
    keys = ["fid", "iid"]
    if not set(keys).issubset(fam) or len(fam) != n_samples:
        raise ValueError("FAM fields/rows do not match the genotype matrix")
    if not set(keys + ["group"]).issubset(groups):
        raise ValueError("Group CSV needs fid, iid, group columns")
    if fam[keys].isna().any().any() or groups[keys + ["group"]].isna().any().any():
        raise ValueError("Sample IDs and group labels cannot be missing")
    left = fam[keys].astype(str).reset_index(drop=True)
    right = groups[keys + ["group"]].copy()
    right[keys] = right[keys].astype(str)
    if left.duplicated(keys).any() or right.duplicated(keys).any():
        raise ValueError("Sample IDs must be unique in FAM and group CSV")
    aligned = left.merge(right, on=keys, how="left", validate="one_to_one", sort=False)
    if aligned["group"].isna().any() or len(aligned) != n_samples:
        raise ValueError("Every FAM sample needs exactly one group assignment")
    return aligned


def _scan(array: Any, n_variants: int, case: np.ndarray, control: np.ndarray,
          batch_size: int) -> pd.DataFrame:
    parts = [
        variant_statistics(_batch(array, start, min(start + batch_size, n_variants)), case, control)
        for start in range(0, n_variants, batch_size)
    ]
    return pd.concat(parts, ignore_index=True)


def run_qc(
    genotypes: Any,
    bim: pd.DataFrame,
    fam: pd.DataFrame,
    groups: pd.DataFrame,
    config: QCConfig = QCConfig(),
    *,
    case_label: str = "case",
    control_label: str = "control",
    batch_size: int = 1024,
) -> QCResult:
    """Run variant QC, sample QC, then recompute variants on retained samples.

    The input is variants x samples (0=a0/a0, 1=a0/a1, 2=a1/a1, NaN=missing).
    Arrays with a `compute()` method, including Dask arrays, are read by batch.
    """
    if isinstance(batch_size, bool) or not isinstance(batch_size, int) or batch_size < 1:
        raise ValueError("batch_size must be a positive integer")
    if len(genotypes.shape) != 2 or min(genotypes.shape) == 0:
        raise ValueError("genotypes must have variants and samples")
    n_variants, n_samples = map(int, genotypes.shape)
    metadata = _metadata(bim, n_variants)
    sample_table = _groups(fam, groups, n_samples)
    labels = sample_table["group"].astype(str)
    if set(labels) != {case_label, control_label}:
        raise ValueError("Group assignments must use exactly the case and control labels")
    case = (labels == case_label).to_numpy()
    control = (labels == control_label).to_numpy()

    initial = _scan(genotypes, n_variants, case, control, batch_size)
    preliminary, _ = mark_variant_pass(initial, metadata, config)
    passing = preliminary["pass_qc"].to_numpy(dtype=bool)
    n_passing = int(passing.sum())
    if n_passing == 0:
        raise ValueError("No variants passed preliminary QC; cannot assess sample call rates")

    called = np.zeros(n_samples, dtype=np.int64)
    for start in range(0, n_variants, batch_size):
        end = min(start + batch_size, n_variants)
        selected = passing[start:end]
        if selected.any():
            block = _batch(genotypes, start, end)[selected, :]
            called += np.isfinite(block).sum(axis=0)
    sample_table["preliminary_variants"] = n_passing
    sample_table["called_variants"] = called
    sample_table["call_rate"] = called / n_passing
    sample_table["pass_qc"] = sample_table["call_rate"] >= config.min_sample_call_rate
    retained = sample_table["pass_qc"].to_numpy(dtype=bool)
    if not np.any(case & retained) or not np.any(control & retained):
        raise ValueError("Sample QC removed all cases or all controls")

    final_stats = _scan(genotypes, n_variants, case & retained, control & retained, batch_size)
    final, threshold = mark_variant_pass(final_stats, metadata, config)
    return QCResult(
        variants=final,
        samples=sample_table,
        preliminary_variant_count=n_passing,
        hwe_threshold=threshold,
    )

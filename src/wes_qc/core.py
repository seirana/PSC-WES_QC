"""Statistical primitives. Genotypes are variants x samples, coded 0/1/2/NaN."""

from dataclasses import asdict, dataclass
from math import isfinite

import numpy as np
import pandas as pd


@dataclass(frozen=True)
class QCConfig:
    max_variant_missing_rate: float = 0.05
    min_sample_call_rate: float = 0.95
    min_minor_carriers: int = 2
    hwe_familywise_alpha: float = 0.05
    exclude_duplicate_sites: bool = True
    exclude_indels: bool = False
    # Historical analysis excluded 0.001 < control MAF < 0.1. This unusual
    # criterion is opt-in and should only be used for a justified protocol.
    control_maf_exclusion_band: tuple[float, float] | None = None

    def __post_init__(self) -> None:
        if not 0 <= self.max_variant_missing_rate < 1:
            raise ValueError("max_variant_missing_rate must be in [0, 1)")
        if not 0 < self.min_sample_call_rate <= 1:
            raise ValueError("min_sample_call_rate must be in (0, 1]")
        if (isinstance(self.min_minor_carriers, bool)
                or not isinstance(self.min_minor_carriers, int)
                or self.min_minor_carriers < 0):
            raise ValueError("min_minor_carriers must be a non-negative integer")
        if not 0 < self.hwe_familywise_alpha < 1:
            raise ValueError("hwe_familywise_alpha must be in (0, 1)")
        if self.control_maf_exclusion_band is not None:
            low, high = self.control_maf_exclusion_band
            if not (0 <= low < high <= 0.5):
                raise ValueError("control MAF band must satisfy 0 <= low < high <= 0.5")

    def as_dict(self) -> dict:
        value = asdict(self)
        if value["control_maf_exclusion_band"] is not None:
            value["control_maf_exclusion_band"] = list(value["control_maf_exclusion_band"])
        return value


def exact_hwe_pvalue(heterozygotes: int, homozygote_0: int, homozygote_2: int) -> float:
    """Two-sided conditional exact HWE test (heterozygote probability ordering).

    Use on diploid, autosomal control genotypes. Returns NaN for no calls.
    """
    counts = (heterozygotes, homozygote_0, homozygote_2)
    if any(isinstance(x, bool) or int(x) != x or x < 0 for x in counts):
        raise ValueError("genotype counts must be non-negative integers")
    het, hom0, hom2 = map(int, counts)
    n = het + hom0 + hom2
    if n == 0:
        return float("nan")
    rare = 2 * min(hom0, hom2) + het
    mid = rare * (2 * n - rare) // (2 * n)
    if (mid - rare) % 2:
        mid += 1
    probabilities = np.zeros(rare + 1, dtype=float)
    probabilities[mid] = 1.0
    current_het = mid
    rare_hom = (rare - mid) // 2
    common_hom = n - current_het - rare_hom
    while current_het >= 2:
        probabilities[current_het - 2] = (
            probabilities[current_het] * current_het * (current_het - 1)
            / (4 * (rare_hom + 1) * (common_hom + 1))
        )
        current_het -= 2
        rare_hom += 1
        common_hom += 1
    current_het = mid
    rare_hom = (rare - mid) // 2
    common_hom = n - current_het - rare_hom
    while current_het <= rare - 2:
        probabilities[current_het + 2] = (
            probabilities[current_het] * 4 * rare_hom * common_hom
            / ((current_het + 2) * (current_het + 1))
        )
        current_het += 2
        rare_hom -= 1
        common_hom -= 1
    probabilities /= probabilities.sum()
    observed = probabilities[het]
    return float(min(1.0, probabilities[probabilities <= observed + 1e-12].sum()))


def _counts(genotypes: np.ndarray) -> dict[str, np.ndarray]:
    n = genotypes.shape[1]
    called = np.isfinite(genotypes)
    hom0 = np.count_nonzero(genotypes == 0, axis=1)
    het = np.count_nonzero(genotypes == 1, axis=1)
    hom2 = np.count_nonzero(genotypes == 2, axis=1)
    observed = np.count_nonzero(called, axis=1)
    a0_frequency = np.divide(
        2 * hom0 + het, 2 * observed,
        out=np.full(genotypes.shape[0], np.nan), where=observed > 0,
    )
    return {
        "hom0": hom0, "het": het, "hom2": hom2,
        "missing": n - observed, "observed": observed,
        "missing_rate": (n - observed) / n,
        "maf": np.minimum(a0_frequency, 1 - a0_frequency),
        "minor_carriers": np.where(
            a0_frequency <= 0.5, hom0 + het, hom2 + het,
        ),
    }


def variant_statistics(
    genotypes: np.ndarray, case_mask: np.ndarray, control_mask: np.ndarray
) -> pd.DataFrame:
    """Compute counts, call rates, MAF, and control HWE for one variant batch."""
    matrix = np.asarray(genotypes, dtype=float)
    if matrix.ndim != 2 or matrix.shape[1] == 0:
        raise ValueError("genotypes must be a nonempty variants x samples matrix")
    valid = np.isfinite(matrix) & np.isin(matrix, [0, 1, 2])
    if not np.all(valid | np.isnan(matrix)):
        raise ValueError("genotypes must be 0, 1, 2, or NaN")
    cases = np.asarray(case_mask, dtype=bool)
    controls = np.asarray(control_mask, dtype=bool)
    if cases.shape != (matrix.shape[1],) or controls.shape != cases.shape:
        raise ValueError("sample masks do not match genotype columns")
    if np.any(cases & controls) or not np.any(cases) or not np.any(controls):
        raise ValueError("case and control masks must be disjoint and nonempty")
    case = _counts(matrix[:, cases])
    control = _counts(matrix[:, controls])
    all_counts = _counts(matrix[:, cases | controls])
    columns = {}
    for prefix, values in (("case", case), ("control", control), ("all", all_counts)):
        for name, array in values.items():
            columns[f"{prefix}_{name}"] = array
    columns["control_hwe_p"] = [
        exact_hwe_pvalue(h, a, b)
        for h, a, b in zip(control["het"], control["hom0"], control["hom2"])
    ]
    return pd.DataFrame(columns)


def mark_variant_pass(stats: pd.DataFrame, metadata: pd.DataFrame, config: QCConfig) -> tuple[pd.DataFrame, float]:
    """Apply documented filters globally so HWE correction spans all batches."""
    required = {"chrom", "pos", "snp", "a0", "a1"}
    if not required.issubset(metadata) or len(stats) != len(metadata):
        raise ValueError("BIM metadata is missing fields or does not align with genotypes")
    result = pd.concat([metadata.reset_index(drop=True), stats.reset_index(drop=True)], axis=1)
    duplicate = result.duplicated(["chrom", "pos"], keep=False).to_numpy()
    alleles = result[["a0", "a1"]].astype(str)
    indel = ((alleles["a0"].str.len() != 1) | (alleles["a1"].str.len() != 1)).to_numpy()
    chrom = result["chrom"].astype(str).str.removeprefix("chr")
    autosomal = chrom.isin([str(i) for i in range(1, 23)]).to_numpy()
    result["duplicate_site"] = duplicate
    result["indel"] = indel
    result["hwe_applied"] = autosomal
    base = (
        (result["case_missing_rate"].to_numpy() <= config.max_variant_missing_rate)
        & (result["control_missing_rate"].to_numpy() <= config.max_variant_missing_rate)
        & (result["all_minor_carriers"].to_numpy() >= config.min_minor_carriers)
    )
    if config.exclude_duplicate_sites:
        base &= ~duplicate
    if config.exclude_indels:
        base &= ~indel
    if config.control_maf_exclusion_band is not None:
        low, high = config.control_maf_exclusion_band
        maf = result["control_maf"].to_numpy()
        base &= ~((maf > low) & (maf < high))
    tested = base & autosomal
    threshold = config.hwe_familywise_alpha / max(1, int(np.count_nonzero(tested)))
    hwe_pass = ~autosomal | (
        np.isfinite(result["control_hwe_p"].to_numpy())
        & (result["control_hwe_p"].to_numpy() > threshold)
    )
    result["pass_qc"] = base & hwe_pass
    return result, threshold

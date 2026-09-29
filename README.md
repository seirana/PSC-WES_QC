# PSC-WES_QC

Quality control for diploid whole-exome genotypes in PLINK BED/BIM/FAM format.
The maintained Python package is in `src/wes_qc/`. The original 2023 scripts
are retained **unchanged** under `legacy/` for provenance; they are not the
supported entry point.

## What changed

The historical `QC.py` did not run as a normal Python program: it used an
IPython reset, imported non-existent modules, embedded one person's filesystem
paths and cohort counts, and ended with unconditional directory deletion. The
maintained workflow accepts explicit input paths, sample groups and thresholds.
It never deletes files. See [MIGRATION.md](MIGRATION.md) for the scientific and
behavioral differences.

## Install

Python 3.10 or newer:

```bash
python -m venv .venv
source .venv/bin/activate
python -m pip install -e '.[plink]'
```

For the synthetic tests alone, `python -m pip install -e .` is sufficient.
The full command uses `pandas-plink` to read a BED/BIM/FAM trio. This
repository contains no genotype or patient data.

## Prepare sample groups

Supply a CSV with `fid,iid,group` columns and exactly one row for each FAM
sample. The `fid` and `iid` pair must be unique. For example:

```csv
fid,iid,group
F1,S1,case
F2,S2,control
```

These are fictional IDs. Keep the real mapping private. Group sizes and sample
order are derived from the files, not from hard-coded indices or counts.

## Run

```bash
wes-qc \
  --plink-prefix /private/path/cohort \
  --groups /private/path/groups.csv \
  --output /private/path/qc-output
```

The prefix resolves to `cohort.bed`, `cohort.bim`, and `cohort.fam`.
For a study with different group names, use `--case-label` and
`--control-label`. `--batch-size` controls how many variants are read into
memory at once (default 1024).

The program first filters variants using both groups, calculates sample call
rates among preliminary passing variants, excludes samples below the call-rate
threshold, then recomputes variant statistics on the retained samples. It
raises an error when metadata do not align, group assignments are missing,
genotypes are not 0/1/2/NaN, or a whole group disappears. No cleanup or
destructive mutation of inputs occurs.

### Default criteria

| Criterion | Default |
| --- | --- |
| Maximum variant missingness in **each** group | 5% |
| Minimum sample call rate on preliminary passing variants | 95% |
| Minimum carriers of the less common allele across both groups | 2 |
| Duplicate chromosome/position sites | Excluded |
| Indels | Retained and flagged; use `--exclude-indels` to exclude |
| Control HWE exact test on chromosomes 1–22 | Bonferroni `0.05 / number of tested variants` |
| Sex chromosomes and other contigs | HWE skipped and marked `hwe_applied=false` |
| Control MAF band exclusion | Disabled |

The historical script attempted to exclude variants with
`0.001 < control MAF < 0.1`. That unusual study-specific criterion can be
enabled explicitly with `--exclude-control-maf-band 0.001 0.1`; it is not a
general QC default. Thresholds should be selected and documented by the
investigator for the particular cohort, sequencing design, ancestry, and
analysis. No imputation, relatedness, contamination, sex concordance, batch
effects, or ancestry QC is performed here.

Genotypes are assumed to be diploid calls: 0 = a0/a0, 1 = a0/a1, 2 = a1/a1,
NaN = missing. MAF is the lesser of the a0 and a1 frequencies, so it does
not depend on which allele is named first. A duplicate site is any repeated
`chrom,pos` pair, even across input batches.

## Outputs

The private output directory contains:

- `variants.csv`: identifiers, per-group genotype counts, missing rates, MAF,
  control HWE p-values, duplicate/indel/HWE flags, and `pass_qc`;
- `samples.csv`: FID, IID, group, call rate and `pass_qc`;
- `run_metadata.json`: thresholds, counts, input SHA-256 hashes and Python
  version.

Outputs contain sample identifiers. The repository ignores common genotype,
group-map and output paths, but **do not commit or upload patient-level
files**. Keep results in an access-controlled location and review any other
files before pushing. See [PRIVACY.md](PRIVACY.md).

## Test

```bash
python -m pip install -e .
python -m unittest discover -s tests -v
```

CI tests statistical primitives, input validation, cross-batch duplicate
detection and a chunked two-pass run on synthetic genotypes. It does not test
the reader against a real cohort or validate scientific thresholds for a
particular study. A local pilot run with approved data and an independent
comparison to established QC tools remain necessary before interpreting
research results.

## License

No license has been granted in this repository. Public visibility by itself
does not grant reuse rights.

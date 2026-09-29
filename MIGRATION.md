# Migration from the 2023 scripts

The original `QC.py` and `codes/` are archived under `legacy/` without
editing their research record. Do not run `legacy/QC.py` on a working
directory: its final step deletes hard-coded output directories.

| Historical stage | Maintained behavior |
| --- | --- |
| `stp0` and `stp1_2` | Read BED by bounded variant batches; no intermediate chromosome files |
| `stp1_1` | Use BIM chromosome, position and alleles directly; flag duplicate sites and indels |
| `stp2` | Count 0/1/2/NaN genotypes from explicit case/control assignments |
| `stp3` | Calculate MAF from allele counts, without rewriting IDs or swapping alleles |
| `stp4` and `stp5` | Calculate rates and MAF with observed group denominators |
| `stp6` and `stp7` | Apply variant filters and exact control HWE with a global Bonferroni denominator |
| `stp8` | Calculate sample call rates from preliminary passing variants |
| `stp9` | Recompute variant statistics after sample filtering and produce final flags |
| `stp10` | Write aligned variant and sample tables after successful computation |
| `stp11` | Removed; the maintained code never deletes data |

The old driver also referred to module names that did not exist and omitted
parentheses when invoking the final SNP stage. Its HWE function and valid-SNP
function were named opposite to their filenames. The old code assumed 883
cases and 4509 controls, used local paths and mixed chromosome-local with
global variant indices. Its MAF code stored allele frequency rather than the
lesser allele frequency and could divide by zero. The new package makes these
assumptions and computations explicit.

These changes can alter which variants pass. They are correctness and
reproducibility changes, not a claim that historical study results have been
replicated. Preserve the old outputs and compare a small, approved pilot
cohort before replacing an existing analysis protocol.

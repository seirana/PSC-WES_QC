# Reproducibility

Keep the `run_metadata.json` file beside the local outputs. It records input
hashes, filter settings, group labels, and counts. For a reported analysis,
also record the Git commit, package environment, sequencing and genotype
calling pipeline, PLINK conversion procedure, reference build, cohort
definition, relatedness/ancestry handling, and approval for using the data.

Variant filters depend on sample filters. This implementation uses:

1. Initial variant statistics and global HWE correction on all assigned
   samples.
2. Sample call rates among preliminary passing variants.
3. Final variant statistics and global HWE correction on retained samples.

The second pass does not repeatedly refilter samples to convergence. Both
the preliminary count and final flags are written, so this sequence can be
audited. HWE is tested only on autosomal control genotypes; sex chromosomes
and other contigs carry a `hwe_applied=false` flag.

Input data and output tables are not committed. A SHA-256 hash identifies
the input bytes without embedding genotype data in the repository. Hashes
alone do not document biological provenance; keep the source and processing
history separately.

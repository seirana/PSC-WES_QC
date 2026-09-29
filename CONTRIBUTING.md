# Contributing

Keep maintained code in `src/wes_qc/` and preserve `legacy/` as an
unchanged historical snapshot. New filters need an explicit configuration
option, documented denominator and behavior for missing values, and a
synthetic regression test. Avoid hard-coded cohort sizes, personal paths,
implicit deletions, and patient-level test fixtures.

Run `python -m unittest discover -s tests -v` before proposing changes.
Document any change to the scientific selection criteria in the pull request
and update the migration/reproducibility notes where relevant.

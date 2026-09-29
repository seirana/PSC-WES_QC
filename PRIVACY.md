# Handling cohort data

This public repository contains source code and synthetic tests only.
PLINK BED/BIM/FAM files, group mappings, variant outputs and sample outputs
may contain sensitive patient-level information. Store them outside the
repository in an approved access-controlled location.

`.gitignore` helps prevent accidental additions, but it does not protect a
file that was already committed or force-added. Before each push, inspect
`git status` and `git diff --cached --name-only`, and confirm that no
patient data, credentials or private paths are staged. Seek the appropriate
data governance review before sharing any cohort-level results.

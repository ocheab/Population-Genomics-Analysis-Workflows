# Reproducibility guidance

A reproducible use of the code requires:

1. An independently verifiable definition of sample and variant inclusion.
2. Genomic reference version and coordinate convention.
3. Source data permissions, file paths and checksums.
4. Captured software versions, R `sessionInfo()`, and applicable Python dependencies.
5. Consistent seeds, fold assignments and site selection rules.
6. A clear record of intermediate outputs and upstream requirements.
7. An audit of missing calls, ambiguous alleles and unsupported features.

This public code-only release intentionally excludes restricted or unpublished research data and analyses. The original scripts are kept byte-for-byte unchanged, which means they retain application-specific input assumptions. Therefore, independent end-to-end reproduction is not guaranteed from the repository alone. This should be clearly distinguished from *source-code availability*.

For a local syntax check where R is available:

```bash
for file in R/*.R; do Rscript -e 'invisible(parse(file=commandArgs(TRUE)[1]))' "$file" || exit 1; done
```

Parsing is not an integration test; scientific reproducibility also requires the correct data and comparisons.

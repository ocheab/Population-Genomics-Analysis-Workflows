# Inputs and environment

## Inputs are not included

The scripts expect locally authorised variant/genotype inputs and metadata. Depending on the analysis step, inputs may include compressed VCFs and indexes, sample-level quality files, Fws estimates, group labels, reference sequences, gene annotation files and previously generated RDS objects. Their precise paths are specified in `R/00_config.R` and accessed by `R/01_helpers.R`.

Do **not** commit actual sample identifiers, private metadata, genotype files, study results or assay outputs without separate disclosure and data-use permission.

## R setup

Run all R commands from the repository root to preserve relative `source('R/...')` statements. Install R and `ggplot2` as appropriate for the local operating system. Plotting modules may skip exports if `ggplot2` is absent.

```bash
Rscript --version
Rscript -e 'install.packages("ggplot2", repos="https://cloud.r-project.org")'
```

Specify the local project root with `PI4K_ROOT`, if different from the default configured by the original scripts. The original configuration also contains hard-coded file names and expected data properties; using unrelated data may require a separate compatible configuration. Changing the original scripts is outside the scope of this exact-source release.

## Python setup

Python scripts can be imported byte-for-byte from a local project using `tools/import_python_from_local_project.sh`. Dependency installation must follow the requirements of the original imported scripts; no Python workflow or dependency lockfile was supplied in the present source archive.

## Security and provenance

Record checksums, software versions and data-source versions in private research records. Source files do not carry a data-use licence, and no authorization to redistribute primary datasets is implied by public code availability.

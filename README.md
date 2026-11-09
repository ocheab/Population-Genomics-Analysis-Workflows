# Population Genomics Analysis Workflows (R and Python)

This repository publishes **analysis source code** and general technical documentation for genotype quality assessment, population-genetic comparisons, and exploratory marker-based prediction. It is a methods/software resource, **not** a dataset, thesis-results repository, or validated classification system.

## What is included

- `R/`: exact, unmodified R source files from the archived analysis workflow.
- `python/`: location for exact Python processing scripts, to be imported from the original project using the supplied local helper; Python script files were not present in the source archive used to assemble this package.
- `docs/`: general input guidance, method descriptions, execution notes, limitations, and reproducibility instructions.
- `tools/`: local Python-import and release-audit helpers (these are packaging utilities, not the original analysis scripts).

**No observed genotype datasets, individual-level records, tabulated findings, model predictions, figure outputs, or thesis results are distributed.**

## Scope and workflow

The code provides methods for:

1. Variant-file ingestion and sample/site quality evaluation.
2. Haploid-consensus allele handling and missing/ambiguous-call tracking.
3. Allele frequencies, polymorphism summaries, rarefaction, and within-group diversity.
4. Exploratory within-region LD, group differentiation, permutation testing, and bootstrap summaries.
5. Train/test separation, feature selection, nested cross-validation, and benchmark model comparisons.
6. Sensitivity analyses, assay-preparation aids, figure export, and provenance auditing.

Not every stage runs independently. Scripts have file dependencies and must be executed from their original expected working directory, with suitable input data available locally. The R scripts are **preserved exactly**, including original file paths, study-specific checks, variable names, and assertions; they have **not** been converted into a general-purpose package or validated on arbitrary external datasets.

## Requirements

- R with base/recommended packages; `ggplot2` is used for plotting.
- Python 3 for preparation utilities; original Python analysis dependencies depend on the scripts imported later.
- Locally accessible variant data, sample metadata, reference resources and upstream intermediate files compatible with the configuration.
- A working R environment, correct reference coordinates and an appropriate data-use licence.

See [input and environment guide](docs/INPUTS_AND_ENVIRONMENT.md) for setup assumptions, [code map](docs/CODE_MAP.md) for the original R script responsibilities, and [statistical methods](docs/STATISTICAL_METHODS.md) for general interpretation.

## Getting started (with your own authorised input files)

```bash
# From the repository root
Rscript --version
Rscript -e 'cat("ggplot2:", requireNamespace("ggplot2", quietly=TRUE), "\n")'

# Provide the root directory expected by the unmodified R configuration.
# The scripts also use hard-coded data paths and checks defined in R/00_config.R.
export PI4K_ROOT=/absolute/path/to/local/analysis/project

# Inspect the configuration and requirements first.
sed -n '1,120p' R/00_config.R

# Run data preflight only after all required inputs exist.
Rscript R/02_preflight_qc.R
```

The displayed commands are operational examples, **not** claims that the repository ships the required data or is directly runnable without project setup. See [RUN_ORDER](docs/RUN_ORDER.md) before running downstream scripts.

## Python scripts

The supplied archive contained R sources and research outputs, not Python source files. To copy the **original, byte-identical** `.py` scripts from an existing local project into this clean repository, run:

```bash
bash tools/import_python_from_local_project.sh /absolute/path/to/local/project/scripts
python3 tools/check_release_contents.py
```

Review each imported file for credentials, restricted information, and embedded findings before making the repository public. The import script never edits the `.py` contents.

## Reproducibility and interpretation

Code availability alone does not establish reproducibility: input data access, reference versions, sample definitions, software dependencies and complete processing provenance must also be documented for an executable study replication. The presence of a statistical test or a prediction routine in source code does not validate its assumptions or generalisability. Consult [REPRODUCIBILITY](docs/REPRODUCIBILITY.md) and [LIMITATIONS](docs/LIMITATIONS.md).

## Code availability and citation

Exact source copies are available in this repository. No empirical results or third-party data are released here. Versioned releases and software citations can be added after authorship, institutional ownership, and licensing have been confirmed. **No open-source licence is granted by this repository unless a licence file is deliberately added.**

## Responsible release

Some scripts retain the original application-specific names, constants and paths. **A generic README does not anonymise the actual code.** Review all source files and Git history before publication. See [PUBLIC_RELEASE_CHECKLIST](docs/PUBLIC_RELEASE_CHECKLIST.md).

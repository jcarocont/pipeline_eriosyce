# pipeline-landscape-eriosyce

Landscape genomics pipeline for *Eriosyce* subgen. *Eriosyce* — from raw DArT SNPs to gene flow, population structure, and genotype-environment association.

![Pipeline overview](pipeline_diagram.png)

---

## Overview

The pipeline integrates eight analytical steps orchestrated by a single Fish shell script (`pipeline-thesis.fish`). Each step is an independent R script callable via `Rapp`. Steps 6 and 8 call external CLI tools (IQ-TREE2, BayesAss3).

The turnover/gradient forest analysis (step 4) uses a custom R package developed from `turnover_functions.R`.

---

## Dependencies

### R packages
- `dartRverse`, `dartR.base`, `dartR.popgenomics`
- `LEA`, `tess3r`
- `ranger`, `vegan`, `ecodist`
- `tidyverse`, `furrr`, `ggplot2`, `patchwork`, `sf`, `rnaturalearth`
- `eulerr`, `fields`

### External tools
- [IQ-TREE2](http://www.iqtree.org/) — must be in `$PATH` as `iqtree2`
- [BayesAss3](https://github.com/stevemussmann/BayesAss3-SNPs) — must be in `$PATH` as `BA3`
- [PLINK](https://www.cog-genomics.org/plink/) — binary expected at `./plink`

### Local
- `turnover_functions.R` — source of the GF turnover package (included in repo)

---

## Inputs

All inputs go in the working directory. Edit the `CONFIG` and `INPUTS` block at the top of `pipeline-thesis.fish`:

| Variable | Description |
|---|---|
| `raw_dart.csv` | Raw DArT SNP file |
| `covariables.csv` | Individual metadata / population assignments |
| `envar_final-corto.csv` | Environmental variables per locality (rows = localities, cols = variables) |
| `coords.csv` | Geographic coordinates (`id`, `lat`, `lon`) |

---

## Running

```fish
# Full pipeline
fish pipeline-thesis.fish

# Individual step (example)
Rapp lectura_et_filtrado.r \
    --filename    raw_dart.csv \
    --covfilename covariables.csv \
    --outname     results/run1 \
    --savepath    results/
```

Each script exposes its arguments as `--param` flags (parsed via `Rapp`'s `#| description:` header convention).

---

## Pipeline steps

| Step | Script(s) | What it does |
|---|---|---|
| 01 | `lectura_et_filtrado.r` | Read depth, MAF, call rate, LD, imputation filters on DArT data |
| 02 | `run_tess3r.r` + `plot_tess3r.r` | Spatial population structure via TESS3; cross-entropy K selection |
| 03 | `run_lfmm2.r` | Genotype–environment association with latent factor mixed models |
| 04 | `random_forest_model.R` + `composite_turnover_full.R` | Per-SNP random forests; genomic turnover curves along environmental gradients |
| 05 | `selective_loci_lfmm-gf.r` + `format_str_by_category.r` | Classify loci (LFMM ∩/∪ GF, neutral); export per-category STR and FASTA |
| 06 | `iqtree2` (CLI) | Maximum-likelihood phylogenetics per loci category |
| 07 | `multivariate_comparation_tests.r` | Mantel, MRM, db-RDA, variance partitioning (LFMM / GF / geography) |
| 08 | `BA3` (CLI) + `plot_bayesass_panel.r` | Contemporary gene flow estimation; comparative migration maps |

---

## Outputs

The pipeline writes to step-specific subdirectories:

```
step1_gene/      # filtered VCF, STR, FASTA
step2_figs/      # TESS3 structure plots
step3_textfiles/ # LFMM results, SNP lists
step4_figs/      # turnover curves, diagnostics
step4_textfiles/ # SNP classifications, variable importance
step5_gene/      # per-category STR and FASTA
step6_iqtree/    # trees and distance matrices
step7_figs/      # Euler variance partitioning diagram
step8_figs/      # BayesAss migration panel
```

> **Note:** RDS objects, model files, and result CSVs are excluded from version control (`.gitignore`). The repo tracks code and inputs only.

---

## Repository structure

```
.
├── pipeline-thesis.fish          # master orchestrator
├── lectura_et_filtrado.r
├── run_tess3r.r
├── plot_tess3r.r
├── run_lfmm2.r
├── random_forest_model.R
├── composite_turnover_full.R
├── turnover_functions.R          # GF turnover package source
├── selective_loci_lfmm-gf.r
├── format_str_by_category.r
├── multivariate_comparation_tests.r
├── plot_bayesass_panel.r
└── pipeline_diagram.png
```

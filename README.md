# Cellular Senescence Classification from scRNA-seq/snRNA-seq Data

This repository contains the core analysis workflow used to identify and classify senescent cells.

This workflow is adapted from https://github.com/BeckmannLab/brainCellularSenescenceAndStructure 

Original study:
Lund et al. *Establishing the relationship between brain cellular senescence and brain structure.* Cell (2026).  


## Inputs

The workflow expects:

- a Seurat object with raw RNA counts
- cell type annotations
- sample IDs
- a curated canonical senescence gene list
- a threshold grid defining candidate percentiles to test

## Main scripts

- `scripts/01_score_senescence.R`  
  Compute AUCell senescence scores for each cell

- `scripts/02_threshold_cells.R`  
  Apply candidate thresholds within each cell type and generate binary senescence labels

- `scripts/03_make_pseudobulk.R`  
  Aggregate counts into pseudobulk profiles by senescence label and sample

- `scripts/04_run_de.R`  
  Run differential expression using edgeR/dreamlet

- `scripts/05_select_best_threshold.R`  
  Test enrichment of canonical senescence genes and select the best threshold per cell type

- `scripts/06_pull_final_outputs.R`  
  Assemble the final output objects from the selected thresholds

## Expected outputs

Successful implementation of this workflow generates:

1. continuous senescence activity scores per cell
2. optimized threshold-based senescence classifications within each annotated cell population
3. cell type-specific senescence differential expression signatures
4. estimated proportions of senescent cells for each cell type


## Core dependencies

Main R packages used in this workflow:

- Seurat
- AUCell
- dreamlet
- edgeR
- limma
- BiocParallel
- GSEABase
- dplyr
- tidyr
- stringr
- foreach
- doParallel

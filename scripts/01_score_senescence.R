##Preparation
#1.	Load necessary packages and functions.

library(Seurat)
library(AUCell)
library(dplyr)
library(BiocParallel)
library(GSEABase)
library(dreamlet)
library(foreach)
library(doParallel)
library(stringr)
library(tidyr)
library(limma)
source("https://raw.githubusercontent.com/BeckmannLab/brainCellularSenescenceAndStructure/main/functions/dreamletCompareClusters_edgeR.R") #load in customized dreamlet function for DE analysis that is robust to sparse data
source("https://raw.githubusercontent.com/BeckmannLab/brainCellularSenescenceAndStructure/main/functions/processOneAssay_edgeR.R")# load in second customized dreamlet function for DE analysis that is robust to sparse data

#2. Load your Seurat object
Seurat_object = readRDS("PATH/your_seurat_object.rds")

#3.	Extract raw counts
data = GetAssayData(Seurat_object, slot = "counts", assay = "RNA")

#4.	Define user-specified canonical senescence gene set (literature derived).
sen_genes <-c("4-HNE", "AXL", "BCL2", "CCL2", "CCL3", "CCL4", "CCL5", "CDKN1A", "CDKN2A", "CDKN2B", "CDKN2D", "CSF1", "CSF2RA", "CXCL1", "CXCL8", "GLB1", "H2AX", "HMGB1", "IGF1", "IL1A", "IL1B", "IL27", "IL6", "LGALS3", "LGALS3BP", "LMNB1", "MACROH2A1", "MIF", "MMP12", "MMP3", "MTOR", "PCNA", "PLAUR", "SA-β-Gal", "SATB1", "SERPINE1", "SPP1", "STING1", "TGFB1", "TIMP2", "TNF", "TP53","ACVR1B", "ANG", "ANGPT1","ANGPTL4", "AREG", "AXL", "BEX3", "BMP2", "BMP6", "C3", "CCL1", "CCL13", "CCL16", "CCL2", "CCL20", "CCL24", "CCL26", "CCL3", "CCL3L1", "CCL4", "CCL5", "CCL7", "CCL8", 
"CD55", "CD9", "CSF1", "CSF2", "CSF2RB", "CST4", "CTNNB1","CTSB", "CXCL1", "CXCL10", "CXCL12", "CXCL16","CXCL2", "CXCL3", "CXCL8", "CXCR2", "DKK1", "EDN1","EGF", "EGFR", "EREG", "ESM1", "ETS2", "FAS", "FGF1","FGF2", "FGF7", "GDF15", "GEM", "GMFG", "HGF", "HMGB1","ICAM1", "ICAM3", "IGF1", "IGFBP1", "IGFBP2", "IGFBP3", "IGFBP4", "IGFBP5", "IGFBP6", "IGFBP7", "IL10","IL13", "IL15", "IL18", "IL1A", "IL1B", "IL2", "IL32", "IL6", "IL6ST", "IL7", "INHA", "IQGAP2", "ITGA2", "ITPKA","JUN", "KITLG", "LCP1", "MIF", "MMP1", "MMP10","MMP12", "MMP13", "MMP14", "MMP2", "MMP3", "MMP9", "NAP1L4", "NRG1", "PAPPA", "PECAM1", "PGF", "PIGF", "PLAT","PLAU", "PLAUR", "PTBP1", "PTGER2", "PTGES", "RPS6KA5","SCAMP4", "SELPLG", "SEMA3F", "SERPINB4", "SERPINE1", "SERPINE2", "SPP1", "SPX", "TIMP2", "TNF","TNFRSF10C", "TNFRSF11B", "TNFRSF1A", "TNFRSF1B", "TUBGCP2", "VEGFA", "VEGFC", "VGF", "WNT16", "WNT2")
all_genes = unique(c(sen_genes))

geneSets = GeneSet(all_genes, setName = "sen_list")

#5.	Define the threshold_ grid containing the potential senescence proportion thresholds to be tested for each cell type for downstream workflow. These threshold values and cell type combinations can be customized based on the biological question and dataset characteristics.

# Define the AUCell threshold grid used to classify senescence-positive cells.
# Each row corresponds to a unique combination of:
# (1) percent: descriptive label for the proportion of highest-scoring cells retained,
# (2) celltype: cell type in which the threshold was applied, and
# (3) threshold: AUCell quantile cutoff used to define positive cells.
#
# Thresholds were evaluated across five cutoffs:
# 0.99 (top 1%), 0.95 (top 5%), 0.90 (top 10%),
# 0.80 (top 20%), and 0.70 (top 30%).
#
# The resulting 'threshold_grid' object was used in downstream analyses to iteratively define AUCell-positive cells within each cell type.

threshold_grid <- structure(list(percent = c("1percent","1percent","1percent","1percent","1percent","1percent","1percent","5percent","5percent","5percent","5percent","5percent","5percent","5percent","10percent","10percent","10percent","10percent","10percent","10percent","10percent","20percent","20percent","20percent","20percent","20percent","20percent","20percent","30percent","30percent","30percent","30percent","30percent","30percent","30percent"), celltype = c("MG","Exc","Oli","Int","NonNeu","Ast","OPC","MG","Exc","Oli","Int","NonNeu","Ast","OPC","MG","Exc","Oli","Int","NonNeu","Ast","OPC","MG","Exc","Oli","Int","NonNeu","Ast","OPC","MG","Exc","Oli","Int","NonNeu","Ast","OPC"), threshold = c(0.99,0.99,0.99,0.99,0.99,0.99,0.99,0.95,0.95,0.95,0.95,0.95,0.95,0.95,0.90,0.90,0.90,0.90,0.90,0.90,0.90,0.80,0.80,0.80,0.80,0.80,0.80,0.80,0.70,0.70,0.70,0.70,0.70,0.70,0.70)), row.names = c(NA, -35L), class = "data.frame")

##Step-by-step details [section]

#1. Run AUCell to compute continuous senescence activity scores for each cell using a raw gene expression count matrix and a curated senescence gene set. The resulting AUCX2 object contains per-cell AUCell enrichment scores, where higher scores indicate greater senescence-associated gene activity.

AUCX2 <- AUCell_run(data,geneSets = geneSets,BPPARAM  = BiocParallel::MulticoreParam(30)) #this uses 30 cores on the cluster, for non-parallel processing, use SerialParam()

#2 Save the AUCell object for downstream threshold testing, pseudobulk aggregation, and differential expression analyses. 
saveRDS(AUCX2, "AUCX2_sen_list.rds")

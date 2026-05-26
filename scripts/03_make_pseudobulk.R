#Generate pseudobulk expression profiles
percentages = c("1percent","5percent","10percent", "20percent", "30percent")

# Step 10. Define threshold-specific senescence labels and extract Sample_ID from cell barcodes.
>all_sene = final_df %>%
        mutate(Sample_ID = str_extract(cell, "^[^_]+"))

folder = "path"
# Step 11. Convert the Seurat object to a SingleCellExperiment #object.
>sce = as.SingleCellExperiment(Seurat_object)

foreach(i = percentages) %do% {
# Step 12. Subset cells for the current threshold and construct cell type-specific senescence labels.
    sene = all_sene[which(all_sene$percent == i),]
    sene$cell_sen = paste0(sene$celltype,"_",sene$sen)
    sen3 = sene[,c("cell_sen","cell","Sample_ID"), drop = FALSE]
    pb <- list()
# Step 13. Match AUCell-scored cells to the SingleCellExperiment object and assign metadata.
    matches <- match(colnames(sce), sen3$cell)
    valid_matches <- !is.na(matches)
    sce2 <- sce[, valid_matches]
    colData(sce2)$sen_auc <- as.character(sen3[matches[valid_matches], "cell_sen", drop = TRUE])
    colData(sce2)$Sample_ID <- as.character(sen3[matches[valid_matches], "Sample_ID", drop = TRUE])
    print(i)
    file = paste0(folder, paste0("expression/pseudobulk_", i, ".RDS"))
    if (!file.exists(file)){
# Step 14. Aggregate raw counts into pseudobulk expression profiles by senescence label and sample.
        pb = aggregateToPseudoBulk(sce2,
            assay = "counts", 
            cluster_id = "sen_auc",
            sample_id = "Sample_ID",
            BPPARAM = MulticoreParam(40))
        saveRDS(pb, file)
# Step 15. Save the threshold-specific pseudobulk object.
     }else{
# Step 15. Load the pseudobulk object if it already exists.
         pb=readRDS(file)
     }
}

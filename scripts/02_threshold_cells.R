#Define threshold AUCell scores within each cell type
final_df = c()
resfinal=foreach(i = 1:nrow(threshold_grid)) %do% {
# Step 1. Subset the Seurat object according to the cell type specified in the threshold_grid grid.
    subset = subset(Seurat_object, celltype == threshold_grid[i,2])
# Step 2. Extract AUCell senescence scores for cells within the selected cell type.
    subset_aucx2 = AUCX2[,colnames(AUCX2) %in% colnames(subset) ]
    AUC_mat <- t(getAUC(subset_aucx2))
    AUC_mat = as.data.frame(AUC_mat)
# Step 3. Compute the AUCell score quantile threshold
    # specified for the current condition.
    threshold <- quantile(AUC_mat$sen_list, threshold_grid[i,3])
# Step 4. Classify cells as senescent (TRUE) or
    # non-senescent (FALSE).
    AUC_mat$sen <- AUC_mat$sen_list >= threshold
# Step 5. Annotate cells with the corresponding cell type, threshold label, and cell barcode.
    AUC_mat$celltype = threshold_grid[i,2]
    AUC_mat$percent = threshold_grid[i,1]
    AUC_mat$cell = rownames(AUC_mat)
    print(i)
# Step 6. Save the threshold-specific AUCell classification results.
    saveRDS(AUC_mat, paste0("path", threshold_grid[i,2],"_146_AUC_mat_top_", threshold_grid[i,1],".RDS"))  
# Step 7. Append results into a cumulative table for downstream pseudobulk aggregation and DE analysis.
    final_df = rbind(AUC_mat, final_df)
}

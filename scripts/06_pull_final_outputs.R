# Step 30. Use the selected threshold for each cell type to generate the final workflow outputs.

# Extract AUCell classifications corresponding to the selected threshold for each cell type.
selected_auc <- do.call(rbind,lapply(seq_len(nrow(best)), function(i) {subset(final_df,celltype == best$celltype[i] & percent == best$percent[i])}))

# Load the DE results corresponding to the selected threshold for each cell type.
selected_de <- do.call(rbind,lapply(seq_len(nrow(best)), function(i){readRDS(paste0("../DE/correct_edger_",best$celltype[i],"_",best$percent[i],"_fit.RDS"))}))

# Compute the estimated proportion of senescent cells for each cell type.
estimated_prop <- aggregate(
    sen ~ celltype,
    data = selected_auc,
    FUN = mean
)

colnames(estimated_prop)[2] <-
    "estimated_proportion_senescent"

# Store final workflow outputs.
final_outputs <- list(
    continuous_scores = AUCX2,
    selected_classification = selected_auc,
    selected_de = selected_de,
    estimated_proportions = estimated_prop,
    best_thresholds = best
)

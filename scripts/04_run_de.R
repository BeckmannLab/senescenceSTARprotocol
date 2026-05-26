#16.	Define grid for downstream DE analyses
folder = "path”

cts  = c("Ast","Exc","Int","MG","NonNeu","Oli","OPC")
ctsl = c("ast","exc","int","mg","noneu","oli","opc")
pcts = c("1percent","5percent","20percent","30percent","10percent")

de_grid <- do.call(rbind, lapply(pcts, function(p) {
  data.frame(
    V1 = rep(paste0("pseudobulk_", p, ".RDS"), length(cts)),
    V2 = rep(p, length(cts)),
    V3 = rep(paste0("pb_", p), length(cts)),
    V4 = if (p == "5percent") c("Ast_FALSE", paste0(cts[-1], "_TRUE")) else paste0(cts, "_TRUE"),
    V5 = if (p == "5percent") c("Ast_TRUE",  paste0(cts[-1], "_FALSE")) else paste0(cts, "_FALSE"),
    V6 = ctsl,
    V7 = rep(0.4, length(cts)),
    stringsAsFactors = FALSE
  )
}))

foreach(i = 1:nrow(de_grid)) %do% {
# Step 17. Check whether the DE result for this threshold/cell type combination already exists; skip if it does.
      if(length(which(paste0("correct_edger_", de_grid[i,6],"_", de_grid[i,2], "_fit.RDS")==list.files("path",recursive=TRUE)))==0){
# Step 18. Load the threshold-specific pseudobulk object.
    pb = readRDS(paste0(folder, de_grid[i,1]))
# Step 19. Extract the number of cells contributing to each pseudobulk profile.
    ncells = as.data.frame(int_colData(pb)$n_cells)
    ncells$celltype = rownames(ncells) 
    ncells2 = ncells %>%
      pivot_longer(!celltype,names_to = "id", values_to = "count")
    ncells2 = as.data.frame(ncells2)
    ncells2$id = gsub("\\.", "_", ncells2$id)
    colnames(ncells2) = c("sen_auc", "Sample_ID", "ncell")
# Step 20. Add subject identifiers and cell-count covariates to pseudobulk metadata.
    metadata(pb)$aggr_means$Indv_ID = metadata(pb)$aggr_means$Sample_ID 
    metadata(pb)$aggr_means$Indv_ID=gsub("R|L", "", metadata(pb)$aggr_means$Indv_ID)
    metadata(pb)$aggr_means$ncells = as.numeric(ncells2$ncell)
# Step 21. Define the senescent versus non-senescent comparison for the current cell type.
    ct.pairs <- c(de_grid[i,4], de_grid[i,5])
# Step 22. Run differential expression using edgeR via dreamlet.
    try(fit <- dreamletCompareClusters_edgeR(pb, ct.pairs, method = "fixed",min.prop = 0.4,useProcessOneAssay_edgeR = TRUE))
# Step 23. Extract the full DE table and save the result.
    fit_top <- as.data.frame(topTags(fit, n = Inf))
    print(i)
    #save
saveRDS(fit_top,paste0("path/correct_edger_", de_grid[i,6],"_", de_grid[i,2], "_fit.RDS"))

}else{
    print(paste(i, "done"))
}}

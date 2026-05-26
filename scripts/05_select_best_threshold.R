# Step 24. For each cell type and threshold condition,load the corresponding DE results and combine them into a cumulative DE table.
all2 = c()
for (i in 1:nrow(de_grid)){
    df = readRDS(paste0("../DE/correct_edger_",de_grid[i,6],"_", de_grid[i,2], "_fit.RDS"))
    df$celltype = de_grid[i,6]
    df$test = de_grid[i,2]
    df$symbol = rownames(df)
    all2 = rbind(df,all2)
}
# Step 25. Define significant DE genes using an FDR threshold of 0.05.
all2$sig = all2$FDR <= 0.05

# Step 26. Format the canonical senescence gene set as a one-column annotation table and label these genes as #canonical.
all_genes = as.data.frame(unique(all_genes)) #all_genes from above canonical senescent genes
colnames(all_genes) = "symbol"
all_genes$sen_gene = TRUE

# Step 27. Merge the combined DE table with the canonical senescence gene annotation table to classify each gene as canonical (TRUE) or non-canonical (FALSE).
sc_with_sen=merge(all2,all_genes, by = "symbol", all.x = T)
sc_with_sen$sen_gene[is.na(sc_with_sen$sen_gene)]= FALSE

# Step 28. Test enrichment of canonical senescence genes among significant DEGs for each cell type and threshold condition.

threshold_enrichment_results <-do.call(rbind,lapply(split(sc_with_sen,interaction(sc_with_sen$celltype,sc_with_sen$test,drop = TRUE)),function(d) {

            # Create contingency table of canonical senescence gene membership versus DEG significance.
            tab <- table(d$sen_gene, d$sig)

            # Skip comparisons lacking sufficient dimensions.
            if (nrow(tab) < 2 || ncol(tab) < 2) return(NULL)

            # Perform Fisher's exact test.
            ft <- try(fisher.test(tab), silent = TRUE)

            if (inherits(ft, "try-error")) return(NULL)

            # Store enrichment statistics for the current threshold and cell type.
            data.frame(
                celltype = unique(d$celltype),
                percent = unique(d$test),
                p_value = ft$p.value,
                odds_ratio = unname(ft$estimate),
                total_expressed = nrow(d),
                degs = sum(d$sig),
                kegg_true = sum(d$sen_gene & d$sig),
                total_expressed_kegg = sum(d$sen_gene),
                stringsAsFactors = FALSE
            )}))

# Step 29. Select the optimal threshold per cell type using the strongest enrichment of canonical senescence genes.
threshold_enrichment_results$odds_ratio <-
    as.numeric(threshold_enrichment_results$odds_ratio)

threshold_enrichment_results$p_value <-
    as.numeric(threshold_enrichment_results$p_value)

best <- do.call(rbind, lapply(
  split(threshold_enrichment_results, threshold_enrichment_results$celltype),
  function(d) {

    # If no thresholds are significant for this cell type,
    # return a placeholder message.
    if (all(d$p_value >= 0.05 | is.na(d$p_value))) {
      return(data.frame(
        celltype = unique(d$celltype),
        percent = NA,
        odds_ratio = NA,
        p_value = NA,
        result = "No senescence identified"
      ))
    }

    # Keep only significant thresholds.
    d <- d[d$p_value < 0.05 & !is.na(d$p_value), ]

    # Rank thresholds by highest odds ratio and,
    # when tied, by lowest enrichment p-value.
    d <- d[order(-d$odds_ratio, d$p_value), ]

    # Select the top-ranked threshold.
    out <- d[1, c("celltype", "percent", "odds_ratio", "p_value"), drop = FALSE]
    out$result <- "Significant enrichment"
    out
  }
))

rownames(best) <- NULL

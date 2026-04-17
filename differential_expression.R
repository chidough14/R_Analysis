# --------------------------------------------------------------------------
# Script: 02_differential_expression.r
# Objective: Run DESeq2 and create a ranked gene list for GSEA/GScore
# --------------------------------------------------------------------------

library(DESeq2)
library(dplyr)
library(ggplot2)
library(EnhancedVolcano) # Great for high-quality DE visualization

# 1. Load the filtered counts from script 01
dds <- DESeqDataSet(readRDS("outputs/raw_counts_filtered.rds"), design = ~ definition)

# 2. Run Differential Expression Analysis
message("--- Running DESeq2 (Wald Test) ---")
# This performs estimation of size factors, dispersions, and the GLM fit
dds <- DESeq(dds)

# 1. Extract results specifically comparing Tumor to Normal
# Contrast format: c("factor_name", "numerator_level", "denominator_level")
res <- results(dds, 
               contrast = c("definition", "Primary solid Tumor", "Solid Tissue Normal"))

# 2. Run Shrinkage
# Since "Solid Tissue Normal" isn't the reference, we use the 'ashr' or 'normal' 
# shrinkage type which supports arbitrary contrasts (apeglm only likes coefficients)
message("--- Running LFC Shrinkage ---")
res_shrunken <- lfcShrink(dds, 
                          contrast = c("definition", "Primary solid Tumor", "Solid Tissue Normal"), 
                          type = "ashr")

# 3. Create the Ranked Gene List
message("--- Creating Ranked Gene List ---")
res_df <- as.data.frame(res) %>%
  mutate(entrez_id = rowData(dds)$entrez_id) %>%
  filter(!is.na(entrez_id)) %>% 
  filter(!is.na(padj))

# Calculate ranking metric: -log10(pvalue) * sign(log2FoldChange)
res_df <- res_df %>%
  mutate(rank_metric = -log10(pvalue) * sign(log2FoldChange)) %>%
  arrange(desc(rank_metric))

# Create the named vector for GSEA
ranked_list <- res_df$rank_metric
names(ranked_list) <- res_df$entrez_id

# 4. Save and Export
saveRDS(res_df, "outputs/de_results_full.rds")
saveRDS(ranked_list, "outputs/ranked_gene_list.rds")
write.csv(res_df, "outputs/tables/differential_expression_results.csv")

# 6. Visualization: Volcano Plot
message("--- Generating Volcano Plot ---")
png("outputs/figures/volcano_plot.png", width = 1000, height = 800)
EnhancedVolcano(res_df,
                lab = rownames(res_df),
                x = 'log2FoldChange',
                y = 'padj',
                title = 'Tumor vs Normal: Colorectal Cancer',
                pCutoff = 1e-05,
                FCcutoff = 1.0,
                pointSize = 2.0,
                labSize = 4.0)
dev.off()

message

#Summary table for report
summary_table <- data.frame(
  Category = c("Total Genes Analyzed", "Significantly Up-regulated", "Significantly Down-regulated", "Non-Significant"),
  Criteria = c("After low-count filter", "LFC > 1 & padj < 0.05", "LFC < -1 & padj < 0.05", "padj >= 0.05 or |LFC| <= 1"),
  Count = c(
    nrow(res_df),
    nrow(res_df %>% filter(padj < 0.05 & log2FoldChange > 1)),
    nrow(res_df %>% filter(padj < 0.05 & log2FoldChange < -1)),
    nrow(res_df %>% filter(padj >= 0.05 | abs(log2FoldChange) <= 1))
  )
)

print(summary_table)


#Handling infinite ranks

# 1. Load your ranked list
ranked_list <- readRDS("outputs/ranked_gene_list.rds")

# 2. Check for Infinite values
inf_count <- sum(!is.finite(ranked_list))
message("Number of infinite values found: ", inf_count)

# 3. Replace Infinite values
if(inf_count > 0){
  # Find the maximum absolute finite value
  max_finite <- max(abs(ranked_list[is.finite(ranked_list)]))
  
  # Replace Inf with max_finite * 1.1 (10% higher)
  ranked_list[ranked_list == Inf] <- max_finite * 1.1
  ranked_list[ranked_list == -Inf] <- -max_finite * 1.1
  
  message("Infinite values replaced with: +/- ", max_finite * 1.1)
}

# 4. Remove any NA names (GSEA will fail if a gene ID is missing)
ranked_list <- ranked_list[!is.na(names(ranked_list))]


# We want to keep the duplicate with the highest absolute rank (most significant)
message("Original list size: ", length(ranked_list))

# Sort by absolute value descending
ranked_list <- ranked_list[order(abs(ranked_list), decreasing = TRUE)]

# Keep only the first occurrence of each Entrez ID
ranked_list <- ranked_list[!duplicated(names(ranked_list))]

# 3. Final Cleaning
# Sort back into standard descending order for GSEA
ranked_list <- sort(ranked_list, decreasing = TRUE)

message("Deduplicated list size: ", length(ranked_list))

# 4. Save the clean list
saveRDS(ranked_list, "outputs/ranked_gene_list.rds")
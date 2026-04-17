# --------------------------------------------------------------------------
# Script: 04_method_comparison.r
# Objective: Calculate Similarity Metrics (Jaccard, Overlap, Spearman)
# --------------------------------------------------------------------------

library(dplyr)
library(tidyr)
library(tibble)

orards <- readRDS("outputs/ora_results.rds")
gseards <- readRDS("outputs/gsea_results.rds")
gscorerds <- readRDS("outputs/gscore_results.rds")
# 1. Load Data
ora_k <- as.data.frame(orards$kegg)
gsea_k <- as.data.frame(gseards$kegg)
gscore_k <- gscorerds$kegg  # Already a dataframe

# 2. Extract Significant Pathways (p.adjust < 0.05)
# Note: For GScore we use 'padj', for others 'p.adjust'
sig_ora <- ora_k %>% filter(p.adjust < 0.05) %>% pull(ID)
sig_gsea <- gsea_k %>% filter(p.adjust < 0.05) %>% pull(ID)
sig_gscore <- gscore_k %>% filter(padj < 0.05) %>% pull(ID)

# 3. Metric Calculation Functions
# These formulas will be the ones you cite in your Methods section
calc_jaccard <- function(a, b) {
  intersection <- length(intersect(a, b))
  union <- length(union(a, b))
  return(if (union == 0) 0 else intersection / union)
}

calc_overlap <- function(a, b) {
  intersection <- length(intersect(a, b))
  min_size <- min(length(a), length(b))
  return(if (min_size == 0) 0 else intersection / min_size)
}

# 4. Generate Similarity Table
similarity_results <- data.frame(
  Comparison = c("ORA vs GSEA", "GSEA vs GScore", "ORA vs GScore"),
  Jaccard = c(
    calc_jaccard(sig_ora, sig_gsea),
    calc_jaccard(sig_gsea, sig_gscore),
    calc_jaccard(sig_ora, sig_gscore)
  ),
  Overlap_Coefficient = c(
    calc_overlap(sig_ora, sig_gsea),
    calc_overlap(sig_gsea, sig_gscore),
    calc_overlap(sig_ora, sig_gscore)
  )
)

# 5. Spearman Rank Correlation (Ranking Consistency)
all_ids <- unique(c(ora_k$ID, gsea_k$ID, gscore_k$ID))

rank_df <- data.frame(ID = all_ids) %>%
  left_join(ora_k %>% select(ID, p_ora = pvalue), by = "ID") %>%
  left_join(gsea_k %>% select(ID, p_gsea = pvalue), by = "ID") %>%
  left_join(gscore_k %>% select(ID, p_gscore = pvalue), by = "ID") %>%
  replace_na(list(p_ora = 1, p_gsea = 1, p_gscore = 1))

# Calculate Spearman Rho
spearman_matrix <- cor(rank_df %>% select(starts_with("p_")), method = "spearman")

# 6. Save Results
write.csv(similarity_results, "outputs/tables/similarity_metrics_kegg.csv", row.names = FALSE)
saveRDS(spearman_matrix, "outputs/tables/spearman_matrix_kegg.rds")

message("--- Comparison Metrics Calculated ---")
print(similarity_results)
print(spearman_matrix)



# Calculate Similarity Metrics for GO Biological Process

# 1. Prepare Data
ora_go_df <- as.data.frame(ora_go)
gsea_go_df <- as.data.frame(gsea_go)
gscore_go_df <- gscore_go


# 2. Extract Significant IDs (p.adjust < 0.05)
sig_ora_go <- ora_go_df %>% filter(p.adjust < 0.05) %>% pull(ID)
sig_gsea_go <- gsea_go_df %>% filter(p.adjust < 0.05) %>% pull(ID)
sig_gscore_go <- gscore_go_df %>% filter(padj < 0.05) %>% pull(ID)


# 3. Calculate GO Metrics
go_similarity <- data.frame(
  Comparison = c("ORA vs GSEA", "GSEA vs GScore", "ORA vs GScore"),
  Jaccard = c(
    calc_jaccard(sig_ora_go, sig_gsea_go),
    calc_jaccard(sig_gsea_go, sig_gscore_go),
    calc_jaccard(sig_ora_go, sig_gscore_go)
  ),
  Overlap_Coefficient = c(
    calc_overlap(sig_ora_go, sig_gsea_go),
    calc_overlap(sig_gsea_go, sig_gscore_go),
    calc_overlap(sig_ora_go, sig_gscore_go)
  )
)

write.csv(go_similarity, "outputs/tables/similarity_metrics_go.csv", row.names = FALSE)
message("--- GO Comparison Metrics Calculated ---")
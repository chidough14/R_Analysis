# --------------------------------------------------------------------------
# Script: 05_comparison_visualizations.r
# --------------------------------------------------------------------------
library(ggvenn)
library(pheatmap)

# --- FIGURE 1: VENN DIAGRAM (KEGG) ---
venn_data_kegg <- list(
  ORA = sig_ora,
  GSEA = sig_gsea,
  GScore = sig_gscore
)

png("outputs/figures/comparison_venn_kegg.png", width = 800, height = 800)
print(ggvenn(venn_data_kegg, fill_color = c("#00AFBB", "#E7B800", "#FC4E07"),
             stroke_size = 0.5, set_name_size = 6))
dev.off()

# --- FIGURE 2: SPEARMAN MATRIX HEATMAP (KEGG) ---
# Use the spearman_matrix we calculated in script 04
png("outputs/figures/comparison_heatmap_kegg.png", width = 700, height = 600)
pheatmap(spearman_matrix, 
         display_numbers = TRUE, 
         color = colorRampPalette(c("white", "#fdd49e", "#ef6548", "#7f0000"))(100),
         main = "Spearman Rank Correlation (KEGG P-values)")
dev.off()

# --- FIGURE 3: OVERLAP BARPLOT ---
# This visualizes the Similarity Table from Script 04
library(ggplot2)
sim_plot_data <- similarity_results %>%
  pivot_longer(cols = c("Jaccard", "Overlap_Coefficient"), names_to = "Metric", values_to = "Value")

png("outputs/figures/comparison_metrics_barplot.png", width = 800, height = 500)
ggplot(sim_plot_data, aes(x = Comparison, y = Value, fill = Metric)) +
  geom_bar(stat = "identity", position = "dodge") +
  theme_minimal() +
  scale_fill_manual(values = c("#2E9FDF", "#FC4E07")) +
  labs(title = "Method Similarity Metrics (KEGG)", y = "Score (0 to 1)")
dev.off()

message("--- Comparison Figures Generated ---")



venn_data_go <- list(
  ORA = sig_ora_go,
  GSEA = sig_gsea_go,
  GScore = sig_gscore_go
)

png("outputs/figures/comparison_venn_go.png", width = 800, height = 800)
print(ggvenn(venn_data_go, fill_color = c("#440154FF", "#21908CFF", "#FDE725FF"),
             stroke_size = 0.5, set_name_size = 6))
dev.off()

# --- 2. SPEARMAN MATRIX HEATMAP (GO) ---
# First, align the data for GO
all_ids_go <- unique(c(ora_go_df$ID, gsea_go_df$ID, gscore_go_df$ID))
rank_df_go <- data.frame(ID = all_ids_go) %>%
  left_join(ora_go_df %>% select(ID, p_ora = pvalue), by = "ID") %>%
  left_join(gsea_go_df %>% select(ID, p_gsea = pvalue), by = "ID") %>%
  left_join(gscore_go_df %>% select(ID, p_gscore = pvalue), by = "ID") %>%
  replace_na(list(p_ora = 1, p_gsea = 1, p_gscore = 1))

spearman_matrix_go <- cor(rank_df_go %>% select(starts_with("p_")), method = "spearman")



png("outputs/figures/comparison_heatmap_go.png", width = 700, height = 600)
pheatmap(spearman_matrix_go, 
         display_numbers = TRUE, 
         color = colorRampPalette(c("white", "#d1e5f0", "#67a9cf", "#2166ac"))(100),
         main = "Spearman Rank Correlation (GO P-values)")
dev.off()

# TOP PATHWAY CONSISTENCY PLOT ---
top_30_ids <- gscore_go_df %>% slice_min(pvalue, n = 30) %>% pull(ID)

comp_plot_data <- rank_df_go %>%
  filter(ID %in% top_30_ids) %>%
  left_join(gscore_go_df %>% select(ID, Description), by = "ID") %>%
  pivot_longer(cols = starts_with("p_"), names_to = "Method", values_to = "Pval") %>%
  mutate(Method = gsub("p_", "", Method)) %>%
  mutate(logP = -log10(Pval))

png("outputs/figures/go_top_comparison_refined.png", width = 1000, height = 800)

ggplot(comp_plot_data, aes(x = reorder(Description, logP), y = logP, fill = Method)) +
  geom_bar(stat = "identity", position = position_dodge(width = 0.8), width = 0.7) +
  coord_flip() +
  scale_fill_manual(values = c("gscore" = "#E64B35FF", "gsea" = "#00A087FF", "ora" = "#4DBBD5FF")) +
  theme_bw() +
  labs(title = "Comparative Significance of Top GScore GO Terms",
       subtitle = "Red bars indicate pathways identified primarily through gene co-expression coordination",
       x = NULL, 
       y = "-log10(Adjusted P-value)") +
  theme(axis.text.y = element_text(size = 9),
        legend.position = "top")

dev.off()



sim_plot_data_go <- go_similarity %>%
  pivot_longer(cols = c("Jaccard", "Overlap_Coefficient"), 
               names_to = "Metric", 
               values_to = "Value")

# Generate the Plot
png("outputs/figures/comparison_metrics_barplot_go.png", width = 800, height = 500)

ggplot(sim_plot_data_go, aes(x = Comparison, y = Value, fill = Metric)) +
  geom_bar(stat = "identity", position = "dodge", color = "black", size = 0.2) +
  theme_minimal() +
  scale_fill_manual(values = c("#440154FF", "#FDE725FF")) + # Using Viridis-style colors for GO
  ylim(0, 1) + # Ensures the scale is consistent (0 to 1)
  labs(title = "Method Similarity Metrics: Gene Ontology (BP)",
       subtitle = "Comparing Intersection (Jaccard) vs. Subset Overlap (Overlap)",
       x = "Method Comparison",
       y = "Coefficient Score",
       fill = "Metric Type") +
  theme(axis.text.x = element_text(face = "bold"),
        legend.position = "bottom")

dev.off()

message("--- GO Comparison Figures Successfully Generated ---")

# Create a summary of the Jaccard and Overlap for both
master_comparison <- rbind(
  similarity_results %>% mutate(Database = "KEGG"),
  go_similarity %>% mutate(Database = "GO")
)

write.csv(master_comparison, "outputs/tables/final_comparison_summary.csv", row.names = FALSE)
print(master_comparison)

# Find pathways significant in GScore (padj < 0.05) but NOT in ORA or GSEA
gscore_exclusive <- rank_df_go %>%
  filter(p_gscore < 0.05 & p_ora > 0.05 & p_gsea > 0.05) %>%
  left_join(gscore_go_df %>% select(ID, Description), by = "ID") %>%
  arrange(p_gscore) %>%
  head(5)

write.csv(gscore_exclusive, "outputs/tables/gscore_exclusive_discoveries.csv")
print(gscore_exclusive)

library(clusterProfiler)
library(org.Hs.eg.db)
library(AnnotationDbi)
library(dplyr)

# 1. Load Data from Script 02
res_df <- readRDS("outputs/de_results_full.rds")
ranked_list <- readRDS("outputs/ranked_gene_list.rds")

# Define "Significant Genes" for ORA (Adjust p < 0.05 and |LFC| > 1)
sig_genes <- res_df %>% 
  filter(padj < 0.05 & abs(log2FoldChange) > 1) %>% 
  pull(entrez_id)

# 2. OVER-REPRESENTATION ANALYSIS (ORA)
message("--- Running ORA ---")
ora_go <- enrichGO(gene = sig_genes, OrgDb = org.Hs.eg.db, ont = "BP", pvalueCutoff = 0.05)
ora_kegg <- enrichKEGG(gene = sig_genes, organism = 'hsa', pvalueCutoff = 0.05)

# Total Number of Significant Results (p.adj < 0.05)=================
ora_results <- readRDS("outputs/ora_results.rds")
ora_go <- ora_results$go
ora_kegg <- ora_results$kegg
total_go_sig <- nrow(as.data.frame(ora_go))
total_kegg_sig <- nrow(as.data.frame(ora_kegg))

# Extract Top 10 Terms for the Text and Table
top_10_go <- as.data.frame(ora_go) %>% head(10)
top_10_kegg <- as.data.frame(ora_kegg) %>% head(10)

message("Total Significant GO BP: ", total_go_sig)
message("Total Significant KEGG: ", total_kegg_sig)
print(top_10_go$Description)
print(top_10_kegg$Description)
#============================================================

# 3. GENE SET ENRICHMENT ANALYSIS (GSEA)
message("--- Running GSEA ---")
set.seed(42)
gsea_go <- gseGO(geneList = ranked_list, OrgDb = org.Hs.eg.db, ont = "BP", pvalueCutoff = 0.05, eps = 0, verbose = FALSE)
gsea_kegg <- gseKEGG(geneList = ranked_list, organism = 'hsa', pvalueCutoff = 0.05)

#'=============Delete==============================
gseares <- readRDS("outputs/gsea_results.rds")
gsea_go <- gseares$go
gsea_kegg <- gseares$kegg
total_gsea_go <- nrow(as.data.frame(gsea_go))
total_gsea_kegg <- nrow(as.data.frame(gsea_kegg))
cat("Total GSEA GO BP:", total_gsea_go, "\n")
cat("Total GSEA KEGG:", total_gsea_kegg, "\n")

# Extract Top 10 Terms for the Text and Table
top_10_gseago <- as.data.frame(gsea_go) %>% head(10)
top_10_gseakegg <- as.data.frame(gsea_kegg) %>% head(10)

gsea_go_clean <- as.data.frame(gsea_go) %>%
  select(ID, Description, NES, pvalue, p.adjust) %>%
  head(10)

gsea_kegg_clean <- as.data.frame(gsea_kegg) %>%
  select(ID, Description, NES, pvalue, p.adjust) %>%
  head(10)
#'=============Delete======================================

library(org.Hs.eg.db)
library(dplyr)
library(reshape2)

# 1. Prepare Data
tumor_samples <- colData(vsd)$barcode[vsd$definition == "Primary solid Tumor"]
tumor_counts <- assay(vsd)[, tumor_samples]

# IMPORTANT: Ensure rownames of counts are Entrez IDs
# We fetch the mapping from the vsd object we created in script 01
entrez_ids <- rowData(vsd)$entrez_id
rownames(tumor_counts) <- entrez_ids

# Filter out rows that didn't map to an Entrez ID
tumor_counts <- tumor_counts[!is.na(rownames(tumor_counts)), ]

# Filter for significant genes only (p < 0.05)
sig_entrez <- res_df %>% filter(padj < 0.05) %>% pull(entrez_id)
tumor_counts_sig <- tumor_counts[rownames(tumor_counts) %in% sig_entrez, ]

message("Genes in Correlation Matrix: ", nrow(tumor_counts_sig))

# 2. Calculate Correlation Matrix
cor_matrix <- cor(t(tumor_counts_sig), method = "pearson")

message("Sample IDs in cor_matrix: ", paste(head(rownames(cor_matrix)), collapse=", "))


kegg_data <- clusterProfiler::download_KEGG("hsa")
kegg_table <- kegg_data$KEGGPATHID2EXTID

kegg_sets <- split(as.character(kegg_table[,2]), as.character(kegg_table[,1]))

# extract GO BP Sets (Fixed logic)
go_annots <- AnnotationDbi::select(org.Hs.eg.db, 
                                   keys = keys(org.Hs.eg.db), 
                                   columns = c("ENTREZID", "GO"), 
                                   keytype = "ENTREZID")

# Filter for BP and create a clean list
go_bp_only <- go_annots %>% filter(!is.na(GO)) 
go_sets <- split(as.character(go_bp_only$ENTREZID), 
                 as.character(go_bp_only$GO))


# 4. Improved GScore Engine
run_gscore_analysis <- function(gene_set_list, cor_mat, threshold = 0.5) {
  all_genes <- rownames(cor_mat)
  if(length(all_genes) < 2) return(NULL)
  
  total_pairs <- (length(all_genes) * (length(all_genes) - 1)) / 2
  coexp_universe <- sum(abs(cor_mat[lower.tri(cor_mat)]) > threshold, na.rm = TRUE)
  
  message("Universe total pairs: ", total_pairs)
  message("Universe co-expressed pairs: ", coexp_universe)
  
  results <- lapply(names(gene_set_list), function(gs_name) {
    gs_genes <- intersect(as.character(gene_set_list[[gs_name]]), all_genes)
    
    # Lowered threshold to 5 genes to see if we get any results initially
    if(length(gs_genes) < 5) return(NULL) 
    
    gs_cor <- cor_mat[gs_genes, gs_genes]
    gs_pairs <- (length(gs_genes) * (length(gs_genes) - 1)) / 2
    gs_coexp <- sum(abs(gs_cor[lower.tri(gs_cor)]) > threshold, na.rm = TRUE)
    
    p_val <- phyper(gs_coexp - 1, gs_pairs, total_pairs - gs_pairs, coexp_universe, lower.tail = FALSE)
    
    data.frame(ID = gs_name, 
               Pathway_Size = length(gs_genes),
               CoExp_Pairs = gs_coexp, 
               pvalue = p_val)
  })
  
  res <- do.call(rbind, results)
  if(!is.null(res) && nrow(res) > 0) {
    res$padj <- p.adjust(res$pvalue, method = "BH")
    return(res %>% arrange(pvalue))
  }
  return(NULL)
}

# 4. Run Analysis with a lower gene requirement
# Sometimes pathways in your data are smaller than the standard 10
gscore_kegg <- run_gscore_analysis(kegg_sets, cor_matrix, threshold = 0.5)
gscore_go <- run_gscore_analysis(go_sets, cor_matrix, threshold = 0.5)

if(is.null(gscore_kegg)) {
  message("KEGG is still NULL. Checking overlap...")
  overlap_check <- intersect(unlist(kegg_sets), rownames(cor_matrix))
  message("Number of KEGG genes found in your data: ", length(overlap_check))
}

kegg_names <- as.data.frame(kegg_data$KEGGPATHID2NAME)
colnames(kegg_names) <- c("ID", "Description")

# Merge names into your results
gscore_kegg <- gscore_kegg %>%
  left_join(kegg_names, by = "ID")


library(GO.db)
go_names <- AnnotationDbi::select(GO.db, 
                                  keys = as.character(gscore_go$ID), 
                                  columns = "TERM", 
                                  keytype = "GOID")
gscore_go <- gscore_go %>%
  left_join(go_names, by = c("ID" = "GOID")) %>%
  mutate(Description = TERM)

#  Save Results
saveRDS(list(go = gscore_go, kegg = gscore_kegg), "outputs/gscore_results.rds")
write.csv(gscore_kegg, "outputs/tables/gscore_kegg_results.csv")
write.csv(gscore_go, "outputs/tables/gscore_go_results.csv")

message("--- GScore Analysis Complete ---")

#================Delete========================
gscoreres <- readRDS("outputs/gscore_results.rds")
gscore_go <- gscoreres$go
gscore_kegg <- gscoreres$kegg

total_gscore_go <- nrow(as.data.frame(gscore_go))
total_gscore_kegg <- nrow(as.data.frame(gscore_kegg))
cat("Total GSCORE GO BP:", total_gscore_go, "\n")
cat("Total GSCORE KEGG:", total_gscore_kegg, "\n")

# Extract Top 10 Terms for the Text and Table
top_10_gscorego <- as.data.frame(gscore_go) %>% head(10)
top_10_gscorekegg <- as.data.frame(gscore_kegg) %>% head(10)

gscore_go_clean <- as.data.frame(gscore_go) %>%
  select(ID, Description, pvalue, padj) %>%
  head(10)

gsscore_kegg_clean <- as.data.frame(gscore_kegg) %>%
  select(ID, Description, pvalue, padj) %>%
  head(10)
#==============================================


#  Save Objects for Comparison
saveRDS(list(go=ora_go, kegg=ora_kegg), "outputs/ora_results.rds")
write.csv(ora_kegg, "outputs/tables/ora_kegg_results.csv")
write.csv(ora_go, "outputs/tables/ora_go_results.csv")

saveRDS(list(go=gsea_go, kegg=gsea_kegg), "outputs/gsea_results.rds")
write.csv(gsea_kegg, "outputs/tables/gsea_kegg_results.csv")
write.csv(gsea_go, "outputs/tables/gsea_go_results.csv")


#  Quick Visualization for ORA
png("outputs/figures/ora_kegg_dotplot.png", width = 800, height = 600)
dotplot(ora_kegg, showCategory=20) + ggtitle("ORA KEGG: Dotplot")
dev.off()

png("outputs/figures/ora_kegg_barplot.png", width = 800, height = 600)
print(barplot(ora_kegg, showCategory=20) + ggtitle("ORA KEGG: Barplot"))
dev.off()

png("outputs/figures/ora_go_dotplot.png", width = 800, height = 600)
print(dotplot(ora_go, showCategory=20) + ggtitle("ORA GO: Dotplot"))
dev.off()
 
png("outputs/figures/ora_go_barplot.png", width = 800, height = 600)
print(barplot(ora_go, showCategory=20) + ggtitle("ORA GO: Barplot"))
dev.off()

# Quick Visualization for GSEA
# KEGG
png("outputs/figures/gsea_kegg_dotplot.png", width = 800, height = 600)
print(dotplot(gsea_kegg, showCategory=20) + ggtitle("GSEA KEGG: Dotplot"))
dev.off()

png("outputs/figures/gsea_kegg_barplot.png", width = 800, height = 600)
# For GSEA, barplots usually show the Normalized Enrichment Score (NES)
print(df <- as.data.frame(gsea_kegg) %>% slice_max(abs(NES), n=20) %>%
        ggplot(aes(x=reorder(Description, NES), y=NES, fill=p.adjust)) +
        geom_bar(stat="identity") + coord_flip() + theme_minimal() +
        ggtitle("GSEA KEGG: Barplot (Top NES)"))
dev.off()

# GO
png("outputs/figures/gsea_go_dotplot.png", width = 800, height = 600)
print(dotplot(gsea_go, showCategory=20) + ggtitle("GSEA GO: Dotplot"))
dev.off()

png("outputs/figures/gsea_go_barplot.png", width = 800, height = 600)
print(df <- as.data.frame(gsea_go) %>% slice_max(abs(NES), n=20) %>%
        ggplot(aes(x=reorder(Description, NES), y=NES, fill=p.adjust)) +
        geom_bar(stat="identity") + coord_flip() + theme_minimal() +
        ggtitle("GSEA GO: Barplot (Top NES)"))
dev.off()

#. GSCORE PLOTS (Custom ggplot logic) ---

library(stringr)
library(ggplot2)

plot_gscore <- function(data, type="dot", title="GScore") {
  # 1. Filter for the top 20 most significant pathways
  top_data <- data %>% 
    slice_min(pvalue, n=20) %>%
    # Ensure we use 'Description' for the label, fallback to ID if Description is missing
    mutate(Label = ifelse(is.na(Description), ID, Description)) %>%
    # Shorten very long names so they fit on the plot
    mutate(Label = str_trunc(Label, 50))
  
  if(type == "dot") {
    ggplot(top_data, aes(x = CoExp_Pairs / Pathway_Size, 
                         y = reorder(Label, pvalue), 
                         size = Pathway_Size, 
                         color = padj)) +
      geom_point() + 
      scale_color_gradient(low="red", high="blue") +
      theme_minimal() + 
      labs(title = title, x = "Co-expression Density", y = "Pathway Description")
  } else {
    ggplot(top_data, aes(x = CoExp_Pairs, 
                         y = reorder(Label, CoExp_Pairs), 
                         fill = padj)) +
      geom_bar(stat="identity") + 
      scale_fill_gradient(low="red", high="blue") +
      theme_minimal() + 
      labs(title = title, x = "Number of Co-expressed Pairs", y = "Pathway Description")
  }
}

# KEGG
png("outputs/figures/gscore_kegg_dotplot.png", width = 800, height = 600)
print(plot_gscore(gscore_kegg, "dot", "GScore KEGG: Dotplot"))
dev.off()

png("outputs/figures/gscore_kegg_barplot.png", width = 800, height = 600)
print(plot_gscore(gscore_kegg, "bar", "GScore KEGG: Barplot"))
dev.off()

# GO
png("outputs/figures/gscore_go_dotplot.png", width = 800, height = 600)
print(plot_gscore(gscore_go, "dot", "GScore GO: Dotplot"))
dev.off()

png("outputs/figures/gscore_go_barplot.png", width = 800, height = 600)
print(plot_gscore(gscore_go, "bar", "GScore GO: Barplot"))
dev.off()
message("--- Enrichment Analysis Complete ---")
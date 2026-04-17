
# 1. Load Libraries
library(TCGAbiolinks)
library(SummarizedExperiment)
library(DESeq2)
library(biomaRt)
library(org.Hs.eg.db)
library(dplyr)
library(tibble)

# Create directory for outputs
dir.create("outputs/tables", recursive = TRUE, showWarnings = FALSE)
dir.create("outputs/figures", recursive = TRUE, showWarnings = FALSE)

# 2. Data Acquisition (TCGA-COAD and TCGA-READ)
message("--- Downloading TCGA CRC Data ---")
query_crc <- GDCquery(
  project = c("TCGA-COAD", "TCGA-READ"),
  data.category = "Transcriptome Profiling",
  data.type = "Gene Expression Quantification",
  workflow.type = "STAR - Counts"
)

GDCdownload(
  query_crc, 
  method = "api", 
  files.per.chunk = 20, 
  directory = "GDCdata"
)
message("--- Preparing Data (This may take 5-10 mins) ---")
data_crc <- GDCprepare(query_crc, directory = "GDCdata")

# 3. ID Mapping (Ensembl to Gene Symbol & Entrez)
message("--- Mapping Gene Identifiers ---")
gene_metadata <- as.data.frame(rowData(data_crc))

# Remove version numbers from Ensembl IDs for better mapping
gene_metadata$ensembl_clean <- gsub("\\..*$", "", gene_metadata$gene_id)

# Map to Entrez IDs (Essential for KEGG/GO enrichment later)
mapping <- mapIds(org.Hs.eg.db,
                  keys = gene_metadata$ensembl_clean,
                  column = "ENTREZID",
                  keytype = "ENSEMBL",
                  multiVals = "first")

rowData(data_crc)$entrez_id <- mapping

# 4. Pre-filtering
keep <- rowSums(assay(data_crc) >= 10) >= 3
data_filtered <- data_crc[keep, ]

# 5. Normalization: Variance Stabilizing Transformation (VST)
message("--- Performing VST Normalization ---")
dds <- DESeqDataSet(data_filtered, design = ~ definition)
vsd <- vst(dds, blind = FALSE)

# 6. Save Outputs
saveRDS(data_filtered, "outputs/raw_counts_filtered.rds")
saveRDS(vsd, "outputs/vst_normalized_data.rds")

# Export a preview table
write.csv(assay(vsd)[1:100, 1:10], "outputs/tables/vst_normalized_preview.csv")

library(ggplot2)

#Total samples in your analysis
total_samples <- ncol(vsd)

#Count of samples by group (Definition)
sample_counts <- table(vsd$definition)

# Display results
message("Total Samples: ", total_samples)
print(sample_counts)

# 7. Quality Control Figure: PCA Plot
pca_data <- plotPCA(vsd, intgroup = "definition", returnData = FALSE) 

# Generate and save the plot
png("outputs/figures/pca_plot_preprocessing.png", width = 800, height = 600)
print(
  plotPCA(vsd, intgroup = "definition") + 
    theme_minimal() + 
    labs(title = "PCA of CRC Samples",
         subtitle = "Comparison of Primary Tumor vs. Solid Tissue Normal",
         caption = "Top 500 features by variance") +
    theme(plot.title = element_text(face = "bold"))
)
dev.off()

message("--- Preprocessing Complete. Files saved in /outputs ---")
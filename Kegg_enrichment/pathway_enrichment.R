
# Load necessary libraries
library(tidyverse)
library(readr)
library(clusterProfiler)
library(org.Hs.eg.db)
library(biomaRt)

# Load DESeq2 results from CSV
paired_end_deseq2 <- read_csv("significant_DEGs_Deseq2.csv", col_names = TRUE, show_col_types = FALSE)

# Remove the first column if it contains row names or unwanted data
paired_end_deseq2 <- paired_end_deseq2[,-1]

# Connect to Ensembl
ensembl <- useMart("ensembl", dataset = "hsapiens_gene_ensembl")

# Extract the Ensembl IDs and map them to Entrez IDs
gene_info <- getBM(attributes = c("ensembl_gene_id", "entrezgene_id"),
                   filters = "ensembl_gene_id",
                   values = unique(paired_end_deseq2$ensembl_gene_id),
                   mart = ensembl)

# Perform the left join using ensembl_gene_id as the key
paired_end_deseq2_entrez <- left_join(paired_end_deseq2, gene_info, by = "ensembl_gene_id")

# Remove rows without Entrez IDs
paired_end_deseq2_entrez <- paired_end_deseq2_entrez %>% filter(!is.na(entrezgene_id))

# Create a named vector with Entrez IDs as names and log2 fold changes as values
entrez_ids_deseq2 <- setNames(paired_end_deseq2_entrez$log2FoldChange, paired_end_deseq2_entrez$entrezgene_id)

# Perform KEGG enrichment analysis for DESeq2
kegg_enrichment_deseq2 <- enrichKEGG(
  gene = names(entrez_ids_deseq2),
  organism = "hsa",
  pAdjustMethod = "BH",
  qvalueCutoff = 0.05
)

# Convert the KEGG enrichment results to a data frame
kegg_enrichment_df_deseq2 <- as.data.frame(kegg_enrichment_deseq2)

# Save the DESeq2 KEGG enrichment results to a CSV file
write.csv(kegg_enrichment_df_deseq2, file = "kegg_enrichment_deseq2.csv", row.names = TRUE)

# Load EdgeR results from CSV
paired_end_edger <- read_csv("significant_DEGs_edger.csv", col_names = TRUE, show_col_types = FALSE)

# Remove the first column if it contains row names or unwanted data
paired_end_edger <- paired_end_edger[,-1]

# Extract the Ensembl IDs and map them to Entrez IDs for EdgeR
gene_info_edger <- getBM(attributes = c("ensembl_gene_id", "entrezgene_id"),
                         filters = "ensembl_gene_id",
                         values = unique(paired_end_edger$ensembl_gene_id),
                         mart = ensembl)

# Perform the left join using ensembl_gene_id as the key for EdgeR
paired_end_edger_entrez <- left_join(paired_end_edger, gene_info_edger, by = "ensembl_gene_id")

# Remove rows without Entrez IDs for EdgeR
paired_end_edger_entrez <- paired_end_edger_entrez %>% filter(!is.na(entrezgene_id))

# Create a named vector with Entrez IDs as names and log2 fold changes as values for EdgeR
entrez_ids_edger <- setNames(paired_end_edger_entrez$logFC, paired_end_edger_entrez$entrezgene_id)

# Perform KEGG enrichment analysis for EdgeR
kegg_enrichment_edger <- enrichKEGG(
  gene = names(entrez_ids_edger),
  organism = "hsa",
  pAdjustMethod = "BH",
  qvalueCutoff = 0.05
)

# Convert the KEGG enrichment results to a data frame for EdgeR
kegg_enrichment_df_edger <- as.data.frame(kegg_enrichment_edger)

# Save the EdgeR KEGG enrichment results to a CSV file
write.csv(kegg_enrichment_df_edger, file = "kegg_enrichment_edger.csv", row.names = TRUE)

# Plot top 20 KEGG pathways for DESeq2
ggplot(kegg_enrichment_df_deseq2 %>% top_n(20, wt = -pvalue), aes(x = reorder(Description, -pvalue), y = -log10(pvalue))) +
  geom_bar(stat = "identity", fill = "steelblue") +
  coord_flip() +
  labs(title = "Top 20 KEGG Enriched Pathways (DESeq2)",
       x = "Pathway",
       y = "-Log10 P-value") +
  theme_minimal()

# Plot top 20 KEGG pathways for EdgeR
ggplot(kegg_enrichment_df_edger %>% top_n(20, wt = -pvalue), aes(x = reorder(Description, -pvalue), y = -log10(pvalue))) +
  geom_bar(stat = "identity", fill = "red") +
  coord_flip() +
  labs(title = "Top 20 KEGG Enriched Pathways (EdgeR)",
       x = "Pathway",
       y = "-Log10 P-value") +
  theme_minimal()

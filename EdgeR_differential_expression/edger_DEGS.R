# EdgeR Analysis
# Author: Mercy Ndavi
# Date: 2024-08-20

# Install and load required packages
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
BiocManager::install("edgeR")

library(edgeR)
library(dplyr)
library(readxl)
library(ggplot2)
library(biomaRt)

# Read the data file and set column names
data <- read.csv("feature_counts", header = FALSE, sep = "\t")
colnames(data) <- data[2, ]
cleaded_data <- data[, c("Geneid", "chr", "Start", "End", "Strand", "Length", "SRR...")]

# Process count data
count_data_subset <- as.data.frame(lapply(cleaded_data[, 6:ncol(cleaded_data)], as.numeric))
rownames(count_data_subset) <- cleaded_data$Geneid

# Read and prepare metadata
sample_metadata <- read_excel("sample_metadata.xlsx", sheet = 1)
rownames(sample_metadata) <- sample_metadata$sample
group <- factor(sample_metadata$condition)

# Create DGEList and filter
dge <- DGEList(counts = count_data_subset, group = group)
keep <- rowSums(dge$counts >= 10) >= 1
dge <- dge[keep, ]
dge <- calcNormFactors(dge)

# Model design and fitting
design <- model.matrix(~ group)
dge <- estimateCommonDisp(dge)
dge <- estimateTagwiseDisp(dge)
fit <- glmFit(dge, design)
lrt <- glmLRT(fit)

# Adjust for multiple testing
degs <- topTags(lrt, n = Inf)$table

# Volcano plot visualization
ggplot(degs, aes(x = logFC, y = -log10(FDR))) +
  geom_point(aes(color = FDR < 0.05), alpha = 0.8) +
  labs(title = "Volcano Plot", x = "Log Fold Change", y = "-Log10(FDR)") +
  theme_minimal()

# Annotate genes with biomaRt
ensembl <- useMart("ensembl", dataset = "hsapiens_gene_ensembl")
gene_info <- getBM(attributes = c('ensembl_gene_id', 'hgnc_symbol'),
                   filters = 'ensembl_gene_id',
                   values = rownames(degs), mart = ensembl)
degs$ensembl_gene_id <- rownames(degs)
adjusted_with_symbols <- merge(degs, gene_info, by.x = "ensembl_gene_id")

# Save top 50 genes
top_50_genes <- head(adjusted_with_symbols[order(adjusted_with_symbols$FDR), ], 50)
write.csv(top_50_genes, file = "top_50_genes_edger.csv", row.names = TRUE)

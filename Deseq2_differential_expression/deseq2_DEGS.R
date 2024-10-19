#!/usr/bin/env Rscript

# Load necessary libraries
library(dplyr)
library(readxl)
library(DESeq2)
library(ggplot2)
library(biomaRt)

# Load feature counts data
# Replace 'your_feature_counts_file.csv' with the path to your counts data file
data <- read.csv("your_feature_counts_file.csv", header = FALSE, sep = "\t")  # Adjust separator if needed

# Set the second row as column names and clean the data
colnames(data) <- data[2, ]
data <- data[-c(1, 2), ]
rownames(data) <- NULL  # Reset row names

# Specify the desired order of columns (update as necessary)
new_order <- c("Geneid", "chr", "Start", "End", "Strand", "Length", "Sample1", "Sample2")  # Replace sample names

# Create a cleaned DataFrame with existing columns
existing_columns <- intersect(new_order, colnames(data))
cleaned_data <- data[, existing_columns, drop = FALSE]

# Load sample metadata
# Replace 'your_sample_metadata_file.xlsx' with the path to your sample metadata file
sample_metadata <- read_excel("your_sample_metadata_file.xlsx", sheet = 1)

# Set row names to the sample IDs in the metadata
sample_metadata <- as.data.frame(sample_metadata)
rownames(sample_metadata) <- sample_metadata$sample  # Adjust the sample column name as necessary
sample_metadata$sample <- NULL  # Remove redundant sample column

# Verify that the column names match with sample metadata
if (!all(colnames(cleaned_data) %in% rownames(sample_metadata))) {
  stop("Column names in count data do not match sample metadata.")
}

# Prepare count data for DESeq2
count_data_subset <- cleaned_data[, !(colnames(cleaned_data) %in% c("Geneid", "Start", "End", "Strand", "Length"))]
rownames(count_data_subset) <- cleaned_data$Geneid
count_data_subset <- as.data.frame(lapply(count_data_subset, as.numeric))

# Create DESeq2 dataset
dds <- DESeqDataSetFromMatrix(countData = count_data_subset, colData = sample_metadata, design = ~ condition)  # Adjust the design formula as needed

# Remove low-quality reads
keep <- rowSums(counts(dds)) >= 10
dds <- dds[keep, ]

# Set the factor level for conditions
dds$condition <- relevel(dds$condition, ref = "normal")  # Adjust reference level if necessary

# Run DESeq2 analysis
dds <- DESeq(dds)

# View results
res <- results(dds)
summary(res)

# Filter results for significant genes
res_filtered <- res[!is.na(res$padj), ]
res_filtered <- res_filtered[res_filtered$padj < 0.05 & abs(res_filtered$log2FoldChange) > 1, ]

# Create a volcano plot
res_filtered$significance <- ifelse(res_filtered$padj < 0.05 & abs(res_filtered$log2FoldChange) > 1, "Significant",
                                     ifelse(res_filtered$padj < 0.05, "p-value significant",
                                            ifelse(abs(res_filtered$log2FoldChange) > 1, "Fold change significant", "Not significant")))

ggplot(res_filtered, aes(x = log2FoldChange, y = -log10(padj), color = significance)) +
  geom_point(alpha = 0.4, size = 1.75) +
  theme_minimal() +
  labs(title = "Volcano Plot of Differential Expression",
       x = "Log2 Fold Change",
       y = "-Log10 Adjusted P-value") +
  geom_vline(xintercept = c(-1, 1), col = "red", linetype = "dashed")  # Add threshold lines

# Save significant results to a CSV file
write.csv(as.data.frame(res_filtered), "significant_results.csv", row.names = TRUE)

# Use biomaRt for gene annotation (optional)
mart <- useMart("ensembl", dataset = "hsapiens_gene_ensembl")  # Change to your species as needed

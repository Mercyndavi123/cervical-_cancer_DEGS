#!/bin/bash

# Add FastQC to the PATH
export PATH=$PATH:/path/to/FastQC

# Define the directory containing FASTQ files
FASTQ_DIR="path/to/your/fastq_files"
# Define the output directory for FastQC results
OUTPUT_DIR="path/to/your/output_directory"

# Create the output directory if it does not exist
mkdir -p "$OUTPUT_DIR"

# Run FastQC on each FASTQ file in the directory
for FASTQ_FILE in "$FASTQ_DIR"/*.fastq
do
    echo "Running FastQC on $FASTQ_FILE"
    fastqc "$FASTQ_FILE" -o "$OUTPUT_DIR"
done

echo "FastQC analysis complete. Results are in $OUTPUT_DIR"

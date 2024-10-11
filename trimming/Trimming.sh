#!/bin/bash

# Activate the environment for Trim Galore
source activate trim_galore_env

# Check if trim_galore is available
if ! command -v trim_galore &> /dev/null; then
    echo "Error: trim_galore could not be found."
    exit 1
fi

# Define the input directory containing the FASTQ files
INPUT_DIR="."

# Create output directory for trimmed files
OUTPUT_DIR="./trimmed_reads"

# Create output directory if it doesn't exist
mkdir -p "$OUTPUT_DIR"

# Process paired-end reads
echo "Processing paired-end reads..."
for file1 in ${INPUT_DIR}/*_1.fastq; do
    file2="${file1/_1.fastq/_2.fastq}"
    if [[ -f "$file2" ]]; then
        echo "Trimming $file1 and $file2"
        trim_galore --paired --quality 25 --length 36 --output_dir "$OUTPUT_DIR" "$file1" "$file2"
    fi
done

# Process single-end reads
echo "Processing single-end reads..."
for file in ${INPUT_DIR}/*.fastq; do
    if [[ ! $file =~ _1.fastq$ && ! $file =~ _2.fastq$ ]]; then
        echo "Trimming $file"
        trim_galore --quality 25 --length 36 --fastqc --output_dir "$OUTPUT_DIR" "$file"
    fi
done

echo "Trimming complete. Trimmed files are in $OUTPUT_DIR."

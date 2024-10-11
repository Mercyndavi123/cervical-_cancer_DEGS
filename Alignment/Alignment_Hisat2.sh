#!/bin/bash

# Set directories (replace with your own paths)
REFERENCE_DIR="path/to/reference/genome"
FASTQ_DIR="path/to/trimmed/fastq/files"
OUTPUT_DIR="path/to/alignment/output"

# Paths to required tools (replace with your own paths)
hisat2_path="path/to/hisat2"
samtools_path="path/to/samtools"

# HISAT2 index base name
INDEX_BASE="${REFERENCE_DIR}/genome"

# Number of threads
num_threads=8  # Adjust this based on your system's capacity

# Create output directory if it doesn't exist
mkdir -p "${OUTPUT_DIR}"

# Get all trimmed files (both paired and single-end)
for file in ${FASTQ_DIR}/*_trimmed.fq; do
  # Extract base prefix without the extension
  PREFIX=$(basename "$file" | sed 's/_R[12]_trimmed\.fq\|_trimmed\.fq//')

  # Define paired-end files
  R1="${FASTQ_DIR}/${PREFIX}_R1_trimmed.fq"
  R2="${FASTQ_DIR}/${PREFIX}_R2_trimmed.fq"
  SAM_FILE="${OUTPUT_DIR}/${PREFIX}.sam"
  SINGLE="${FASTQ_DIR}/${PREFIX}_trimmed.fq"

  # Check if both paired-end reads exist
  if [[ -f "$R1" && -f "$R2" ]]; then
    echo "Aligning paired-end reads for ${PREFIX}..."

    # Run HISAT2 for paired-end reads with multiple threads
    ${hisat2_path} -x ${INDEX_BASE} -1 ${R1} -2 ${R2} -S ${SAM_FILE} -p ${num_threads}

  # Check if single-end read exists
  elif [[ -f "$SINGLE" ]]; then
    echo "Aligning single-end read for ${PREFIX}..."

    # Run HISAT2 for single-end reads with multiple threads
    ${hisat2_path} -x ${INDEX_BASE} -U ${SINGLE} -S ${SAM_FILE} -


#!/bin/bash

# Define the path to the annotation file (GTF or GFF)
ANNOTATION_FILE="path/to/annotation_file.gtf"

# Define the output directory for single-end and paired-end results
SINGLE_END_OUTPUT_DIR="path/to/single_end_output"
PAIRED_END_OUTPUT_DIR="path/to/paired_end_output"

# Define the number of threads to use
THREADS=4

# Create the output directories if they don't exist
mkdir -p "$SINGLE_END_OUTPUT_DIR"
mkdir -p "$PAIRED_END_OUTPUT_DIR"

# Define the input directories containing BAM files
SINGLE_END_INPUT_DIR="path/to/single_end_bam_files"
PAIRED_END_INPUT_DIR="path/to/paired_end_bam_files"

# Find all BAM files in the input directories
SINGLE_END_FILES=("$SINGLE_END_INPUT_DIR"/*.bam)
PAIRED_END_FILES=($(find "$PAIRED_END_INPUT_DIR" -name "*.bam"))

# Define output file names for the results
SINGLE_END_OUTPUT_FILE="$SINGLE_END_OUTPUT_DIR/single_end_feature_counts.txt"
PAIRED_END_OUTPUT_FILE="$PAIRED_END_OUTPUT_DIR/paired_end_feature_counts.txt"

# Run featureCounts for single-end BAM files
featureCounts -T "$THREADS" -a "$ANNOTATION_FILE" -o "$SINGLE_END_OUTPUT_FILE" "${SINGLE_END_FILES[@]}"

# Print a message when done with single-end files
echo "featureCounts has finished processing single-end files. The results are in $SINGLE_END_OUTPUT_FILE."

# Run featureCounts for paired-end BAM files
featureCounts -T "$THREADS" -p -a "$ANNOTATION_FILE" -o "$PAIRED_END_OUTPUT_FILE" "${PAIRED_END_FILES[@]}"

# Print a message when done with paired-end files
echo "featureCounts has finished processing paired-end files. The results are in $PAIRED_END_OUTPUT_FILE."

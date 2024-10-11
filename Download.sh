#!/bin/bash

# Define the directory where you want to save the downloaded files
DOWNLOAD_DIR="/home/frank.onyambu/mercy/data/cervical_cancer_work/fastqfiles"

# Add sratoolkit to PATH
export PATH=$PATH:~/mercy/data/sratoolkit.3.1.1-ubuntu64/bin

# Create the download directory if it doesn't exist
mkdir -p "$DOWNLOAD_DIR"

# List of SRA accession IDs to download
SRA_ACCESSIONS=(
    )

# Change to the download directory
cd "$DOWNLOAD_DIR"

# Download each SRA accession
for ACCESSION in "${SRA_ACCESSIONS[@]}"; do
    echo "Downloading SRA data for accession: $ACCESSION"
    prefetch "$ACCESSION"
done

# Convert downloaded SRA files to FASTQ format
# Ensure that fasterq-dump is available and configured
for ACCESSION in "${SRA_ACCESSIONS[@]}"; do
    echo "Converting SRA data to FASTQ for accession: $ACCESSION"
    fasterq-dump "$ACCESSION"
done
echo "Data download and conversion completed!"


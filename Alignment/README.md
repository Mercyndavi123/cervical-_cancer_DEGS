The alignment was done using HISAT2, a tool designed for efficiently aligning RNA sequencing (RNA-seq) data to a reference genome.

The script sets up directories for the reference genome, trimmed FASTQ files, and output results, ensuring the output directory exists. It then processes the trimmed FASTQ files, executing HISAT2 for both paired-end and single-end reads.

For paired-end reads, the script utilizes multiple threads to enhance efficiency, specifying the reference genome index and output files. For single-end reads, it adapts the command accordingly. After each alignment, the script confirms completion for that sample, streamlining RNA-seq data analysis.




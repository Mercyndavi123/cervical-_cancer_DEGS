This script is used to clean up sequencing data by trimming FASTQ files. It works with both paired-end and single-end reads to remove low-quality parts from the sequences.

The script looks for paired-end FASTQ files and trims them based on quality and length. It also processes single-end reads with the same quality checks. The cleaned files are saved in a specific output folder for later use.

Make sure Trim Galore is installed and that you have activated the right environment before running the script. You can change the input and output folder paths in the script to match where your files are located.

After running this script, you will have high-quality FASTQ files ready for your analysis.


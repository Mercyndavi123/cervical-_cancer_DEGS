This script is designed to evaluate the quality of FASTQ files using the FastQC tool. It starts by adding FastQC to the system's PATH, which allows the script to run the FastQC program from any location on the system.

The script defines two key directories: the first is the input directory where the FASTQ files are stored, and the second is the output directory where the results of the FastQC analysis will be saved. If the output directory does not already exist, the script automatically creates it to ensure a smooth workflow.

Once the directories are set up, the script proceeds to process each FASTQ file found in the input directory. It runs the FastQC program on each file, generating quality reports that summarize important metrics about the sequencing data, such as the overall quality scores, the presence of adapter sequences, and other potential issues.

When the script has completed its execution, users can find all the FastQC reports in the designated output directory. These reports provide crucial insights into the quality of the sequencing data, helping researchers identify any problems and make informed decisions about further data processing steps.

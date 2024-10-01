#!/bin/bash

# Define the output directory where FastQC results are saved
output_dir="/stor/work/Bio321G_RY_Spring2024/Summer2024/glassfrogs/results"

# Check if there are FastQC output files in the output directory
if ls $output_dir/*.zip 1> /dev/null 2>&1; then
    echo "Found FastQC output files in $output_dir. Running MultiQC..."

    # Run MultiQC in the output directory
    multiqc $output_dir > multiqc.log 2>&1

    echo "MultiQC analysis completed. Report generated in $output_dir."
else
    echo "No FastQC output files found in $output_dir. Please ensure FastQC has completed successfully before running MultiQC."
fi


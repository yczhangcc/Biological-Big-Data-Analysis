#!/bin/bash

# Define input and output directories
input_dir="/stor/work/Bio321G_RY_Spring2024/Summer2024/glassfrogs/trim_results"
output_dir="/stor/work/Bio321G_RY_Spring2024/Summer2024/glassfrogs/again_fastqc"

# Check if there are .fastq files in the input directory
if ls "$input_dir"/*.fastq 1> /dev/null 2>&1; then
    echo "Found .fastq files in $input_dir. Running FastQC..."

    # Run FastQC on all .fastq files in the input directory
    nohup fastqc -o "$output_dir" "$input_dir"/*.fastq > fastqc.log 2>&1 &

    # Wait for all FastQC processes to finish
    wait

    # Check if FastQC completed successfully
    if [ $? -eq 0 ]; then
        echo "FastQC analysis completed. Output saved in $output_dir."
    else
        echo "Error: FastQC analysis failed. Please check fastqc.log for details."
    fi
else
    echo "No .fastq files found in $input_dir. Please check the path and try again."
fi

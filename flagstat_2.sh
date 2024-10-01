#!/bin/bash

# Directory containing the BAM files
bam_dir="/stor/work/FRI_321G_RY_Spring2024/Summer2024/glassfrogs/Cecilia_Summer/bwa/newbwa_result/"
reference="/stor/work/FRI_321G_RY_Spring2024/Summer2024/glassfrogs/Cecilia_Summer/bwa/GlassFrog_TrinityTranscriptome_27March24.fasta"

echo "Using reference: $reference"

for file in "$bam_dir"/*sorted.bam; do
    output_file="$(basename "$file").flagstat"
    
    # Generate flagstat report
    samtools flagstat "$file" > "$output_file"
    
    # Extract mapping percentage
    mapping_percentage=$(grep "mapped (" "$output_file" | awk -F '[()%]' '{print $2}')
    
    # Handle case where mapping_percentage might be empty
    if [ -z "$mapping_percentage" ]; then
        mapping_percentage="N/A"
    fi
    
    echo "File: $file"
    echo "Mapping Percentage: $mapping_percentage"
    echo "Reference used: $reference"
done

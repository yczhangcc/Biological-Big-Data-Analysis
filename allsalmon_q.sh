#!/bin/bash

# Define paths
transcriptome_reference="/stor/work/Bio321G_RY_Spring2024/Summer2024/glassfrogs/data/transcriptome/GlassFrog_TrinityTranscriptome_27March24.fasta"
bam_dir="/stor/work/Bio321G_RY_Spring2024/Summer2024/glassfrogs/Cecilia/salmon"
output_dir="/stor/work/Bio321G_RY_Spring2024/Summer2024/glassfrogs/Cecilia/salmon_results"

# Create output directory if it does not exist
mkdir -p "$output_dir"

# List of samples
samples=(
    "Ep_B-4352"
    # Add more samples if needed
)

# Quantify transcript abundance for each sample
for sample in "${samples[@]}"; do
    bam_file="${bam_dir}/${sample}_combined_sorted.bam"
    output_sample_dir="${output_dir}/${sample}"

    # Check if BAM file exists
    if [[ -f "$bam_file" ]]; then
        echo "Running Salmon quantification on $bam_file..."
        salmon quant -t "$transcriptome_reference" -l A -a "$bam_file" -o "$output_sample_dir"
        echo "Salmon quantification completed for $sample."
    else
        echo "Error: BAM file $bam_file does not exist. Skipping $sample."
    fi
done

echo "All samples processed."

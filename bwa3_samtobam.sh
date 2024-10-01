#!/bin/bash
#Step 3: Convert SAM to BAM, Sort, and Index BAM Files
# Define the directory containing SAM files
sam_dir="/stor/work/Bio321G_RY_Spring2024/Summer2024/glassfrogs/Cecilia/bwa/newbwa_result"

# Ensure the 'sam' directory exists and has SAM files
if [ ! -d "$sam_dir" ]; then
    echo "Error: SAM directory $sam_dir not found."
    exit 1
fi

# Move into the SAM directory
cd "$sam_dir" || exit

# Iterate through each SAM file in the directory
for file in *.sam; do
    # Extract sample name from the SAM file name
    sampleName=$(basename "$file" .sam)

    # Convert SAM to BAM
    samtools view -bS "$file" > "$sampleName".bam

    # Sort BAM file
    samtools sort -o "$sampleName"_sorted.bam "$sampleName".bam

    # Index sorted BAM file
    samtools index "$sampleName"_sorted.bam
done
samtools view -bS Ep_N-4043a_combined.sam > Ep_N-4043a_combined.bam
samtools sort -o Ep_N-4043a_combined_sorted.bam Ep_N-4043a_combined.bam
samtools index Ep_N-3711_combined_sorted.bam

samtools view -bS N-3568_combined.sam > N-3568_combined.bam
samtools sort -o N-3568_combined_sorted.bam  N-3568_combined.bam
samtools index N-3568_combined_sorted.bam

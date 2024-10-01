#!/bin/bash

# Define directories and paths
input_dir="/stor/work/Bio321G_RY_Spring2024/Summer2024/glassfrogs/data/original_sequences_macrogen"
output_dir="/stor/work/Bio321G_RY_Spring2024/Summer2024/glassfrogs/trimfix_results"
adapter_path="/stor/work/Bio321G_RY_Spring2024/Summer2024/glassfrogs/Cecilia/TRIM/Trimmomatic-0.39/adapters/TruSeq3-PE.fa"
trimmomatic_jar="/stor/work/Bio321G_RY_Spring2024/Summer2024/glassfrogs/Cecilia/TRIM/Trimmomatic-0.39/trimmomatic-0.39.jar"

# Create output directory if it does not exist
mkdir -p $output_dir

# List of samples (pair-end samples)
paired_end_samples=(
    "B-4033"
    "B-4161"
    "N-3568"
    "N-3713"
)

# Function to process paired-end samples
process_paired_end_sample() {
    sample=$1
    echo "Processing paired-end sample: $sample"
    
    input_file_r1="${input_dir}/${sample}_1.fastq.gz"
    input_file_r2="${input_dir}/${sample}_2.fastq.gz"
    paired_output_file_r1="${output_dir}/${sample}_1_paired.fastq"
    unpaired_output_file_r1="${output_dir}/${sample}_1_unpaired.fastq"
    paired_output_file_r2="${output_dir}/${sample}_2_paired.fastq"
    unpaired_output_file_r2="${output_dir}/${sample}_2_unpaired.fastq"

    java -jar $trimmomatic_jar PE -phred33 \
        $input_file_r1 $input_file_r2 \
        $paired_output_file_r1 $unpaired_output_file_r1 \
        $paired_output_file_r2 $unpaired_output_file_r2 \
        ILLUMINACLIP:$adapter_path:2:30:10:8 \
        LEADING:3 TRAILING:3 SLIDINGWINDOW:4:20 MINLEN:36 HEADCROP:10
    
    echo "Finished processing paired-end sample: $sample"
}

# Process all paired-end samples
for sample in "${paired_end_samples[@]}"
do
    process_paired_end_sample $sample
done

echo "All samples processed."

# nohup ./trimv3.sh > output.log 2>&1 &

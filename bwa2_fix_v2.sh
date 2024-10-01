#!/bin/bash
# Step 2: Perform Alignment for Each Paired FASTQ File

# Define paths
transcriptome_reference="/stor/work/Bio321G_RY_Spring2024/Summer2024/glassfrogs/Cecilia/bwa/Trinity.fasta"
fastq_dir1="/stor/work/Bio321G_RY_Spring2024/Summer2024/glassfrogs/R1_TRIM"
fastq_dir2="/stor/work/Bio321G_RY_Spring2024/Summer2024/glassfrogs/R2_TRIM"
working_dir="/stor/work/Bio321G_RY_Spring2024/Summer2024/glassfrogs/Cecilia/bwa/newbwa_result"
log_file="${working_dir}/command2_2_history.log"

# Create the working directory if it doesn't exist
mkdir -p "$working_dir"
cd "$working_dir"

# Initialize the log file
echo "Command History Log" > "$log_file"
echo "===================" >> "$log_file"

# List of sample base names
paired_end_samples=(
  
   "B-4033"
   "B-4161"
   "N-3713"
   "N-3568"
   
)

# Perform alignment for each paired sample
echo "Performing alignments..."
for sample in "${paired_end_samples[@]}"; do
    read1="${fastq_dir1}/${sample}_1_paired.fastq"
    read2="${fastq_dir2}/${sample}_2_paired.fastq"

    if [[ -f "$read1" && -f "$read2" ]]; then
        echo "Running BWA MEM on $read1 and $read2..."
        
        # Run BWA MEM
        bwa mem "$transcriptome_reference" "$read1" "$read2" > "${working_dir}/${sample}_combined.sam" 2>> "$log_file"

        # Check if SAM file is empty
        if [[ ! -s "${working_dir}/${sample}_combined.sam" ]]; then
            echo "Warning: SAM file ${sample}_combined.sam is empty."
            echo "Warning: SAM file ${sample}_combined.sam is empty." >> "$log_file"
        fi
    else
        echo "Paired file for $sample is missing, skipping."
        echo "Error: Paired file for $sample is missing, skipping." >> "$log_file"
    fi
done

echo "Alignment process completed."

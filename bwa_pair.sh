#!/bin/bash

# Define the paths for the input data and working directories
transcriptome_reference="/stor/work/Bio321G_RY_Spring2024/Summer2024/glassfrogs/data/transcriptome/GlassFrog_TrinityTranscriptome_27March24.fasta"
fastq_dir="/stor/work/Bio321G_RY_Spring2024/Summer2024/glassfrogs/trim_results"
working_dir="/stor/work/Bio321G_RY_Spring2024/Summer2024/glassfrogs/working_dir"
log_file="${working_dir}/command_history.log"

# Create the working directory if it doesn't exist
mkdir -p "$working_dir"
cd "$working_dir"

# Initialize the log file
echo "Command History Log" > "$log_file"
echo "===================" >> "$log_file"

# Function to log commands and their output
log_command() {
    echo "Command: $@" >> "$log_file"
    "$@" >> "$log_file" 2>&1
    if [ $? -ne 0 ]; then
        echo "Error: Command failed - $@" >> "$log_file"
    else
        echo "Success: Command completed - $@" >> "$log_file"
    fi
    echo "-------------------" >> "$log_file"
}

# Index the transcriptome if not already indexed
if [[ ! -f "${transcriptome_reference}.bwt" ]]; then
    echo "Indexing the transcriptome..."
    log_command bwa index -a bwtsw "$transcriptome_reference"
fi

# Ensure the indexed reference files are in the working directory
echo "Ensuring indexed reference files are in the working directory..."
log_command cp "${transcriptome_reference}."* "$working_dir"

# Update the reference path to point to the working directory
working_reference="${working_dir}/$(basename "$transcriptome_reference")"

# Verify the presence of indexed files
for ext in amb ann bwt pac sa; do
    if [[ ! -f "${working_reference}.${ext}" ]]; then
        echo "Error: Indexed file ${working_reference}.${ext} is missing."
        exit 1
    fi
done

# Perform alignment for each paired FASTQ file
echo "Performing alignments..."
for read1 in "$fastq_dir"/*_R_1_paired.fastq; do
    base=$(basename "$read1" _R_1_paired.fastq)
    read2="${fastq_dir}/${base}_R_2_paired.fastq"

    if [[ -f "$read2" ]]; then
        echo "Running BWA MEM on $read1 and $read2..."
        output_base="${base}_combined"
        
        # Synchronize read pairs if needed
        log_command repair.sh in1="$read1" in2="$read2" out1="${fastq_dir}/${base}_R_1_paired_repaired.fastq" out2="${fastq_dir}/${base}_R_2_paired_repaired.fastq"

        # Run BWA MEM
        log_command bwa mem "$working_reference" "${fastq_dir}/${base}_R_1_paired_repaired.fastq" "${fastq_dir}/${base}_R_2_paired_repaired.fastq" > "${working_dir}/${output_base}.sam"

        # Check if the SAM file is empty
        if [[ ! -s "${working_dir}/${output_base}.sam" ]]; then
            echo "Warning: SAM file ${output_base}.sam is empty."
            echo "Warning: SAM file ${output_base}.sam is empty." >> "$log_file"
            echo "-------------------" >> "$log_file"
        fi
    else
        echo "Paired file for $read1 or $read2 is missing, skipping."
        echo "Error: Paired file for $read1 or $read2 is missing, skipping." >> "$log_file"
        echo "-------------------" >> "$log_file"
    fi
done

# Create a directory for SAM files and move the generated SAM files there
echo "Creating directory for SAM files and moving SAM files..."
mkdir -p "${working_dir}/samFiles"
log_command mv "${working_dir}"/*.sam "${working_dir}/samFiles/"

# Convert SAM to BAM, sort, and index BAM files
echo "Processing SAM files to BAM, sorting and indexing..."
for sam_file in "${working_dir}/samFiles"/*.sam; do
    base=$(basename "$sam_file" .sam)
    
    # Check if SAM file is not empty
    if [[ -s "$sam_file" ]]; then
        # Convert SAM to BAM
        log_command samtools view -bS "$sam_file" > "${working_dir}/samFiles/${base}.bam"

        # Sort BAM file
        log_command samtools sort -o "${working_dir}/samFiles/${base}_sorted.bam" "${working_dir}/samFiles/${base}.bam"

        # Index BAM file
        log_command samtools index "${working_dir}/samFiles/${base}_sorted.bam"
    else
        echo "Warning: SAM file $sam_file is empty, skipping BAM conversion."
        echo "Warning: SAM file $sam_file is empty, skipping BAM conversion." >> "$log_file"
        echo "-------------------" >> "$log_file"
    fi
done

echo "Pipeline completed successfully."

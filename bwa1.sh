#!/bin/bash
#Step 1: Indexing the Transcriptome Reference
# Define paths
transcriptome_reference="/stor/work/Bio321G_RY_Spring2024/Summer2024/glassfrogs/Cecilia/bwa/Trinity.fasta"
working_dir="/stor/work/Bio321G_RY_Spring2024/Summer2024/glassfrogs/Cecilia/bwa/newbwa_result"
log_file="${working_dir}/command_history.log"

# Create the working directory if it doesn't exist
mkdir -p "$working_dir"
cd "$working_dir"

# Initialize the log file
echo "Command History Log" > "$log_file"
echo "===================" >> "$log_file"

# Index the transcriptome if not already indexed
if [[ ! -f "${transcriptome_reference}.bwt" ]]; then
    echo "Indexing the transcriptome..."
    bwa index -a bwtsw "$transcriptome_reference" >> "$log_file" 2>&1

    # Check if indexing was successful
    if [[ $? -eq 0 ]]; then
        echo "Transcriptome indexing completed successfully."
    else
        echo "Error: Transcriptome indexing failed. See log for details."
        exit 1
    fi
else
    echo "Transcriptome reference already indexed."
fi

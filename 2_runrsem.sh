#!/bin/bash

# Define paths and files
data_dir="/stor/work/Bio321G_RY_Spring2024/Summer2024/glassfrogs/Cecilia/bowtie2"
bam_file="$data_dir/nodiscordant_no_partial_sorted.bam"
rsem_ref="/stor/work/Bio321G_RY_Spring2024/Summer2024/glassfrogs/Cecilia/rsem_ali/rsem_ref_index"
output_dir="/stor/work/Bio321G_RY_Spring2024/Summer2024/glassfrogs/Cecilia/rsem_ali/output"

# create
mkdir -p "$output_dir"

# run RSEM
nohup rsem-calculate-expression -p 8 --paired-end \
                               --bam \
                               --estimate-rspd \
                               --append-names \
                               "$bam_file" \
                               "$rsem_ref" "$output_dir/alignment_result" > "$output_dir/rsem.log" 2>&1 &

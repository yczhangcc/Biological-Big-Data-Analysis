#!/bin/bash

# Define paths
SALMON_INDEX_DIR="/stor/work/Bio321G_RY_Spring2024/Summer2024/glassfrogs/Cecilia/salmon/salmon_index"
BAM_FILE="/stor/work/Bio321G_RY_Spring2024/Summer2024/glassfrogs/Cecilia/bwa/newbwa_result/B-4033_combined_sorted.bam"
OUTPUT_DIR="/stor/work/Bio321G_RY_Spring2024/Summer2024/glassfrogs/Cecilia/salmon/newresult"

# Run Salmon quant with BAM file
salmon quant -t $SALMON_INDEX_DIR -l A -a $BAM_FILE -p 2 -o $OUTPUT_DIR/E4033_quant

# Run Salmon quant with reduced threads
#nohup salmon quant -t $SALMON_INDEX_DIR -l A -a $BAM_FILE -p 4 -o $OUTPUT_DIR/Ep_quant_fix > salmon_quant.log 2>&1 &

echo "Salmon quantification completed successfully."
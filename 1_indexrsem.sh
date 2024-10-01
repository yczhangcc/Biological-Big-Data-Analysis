#!/bin/bash

# Define paths
GENOME_FASTA="/stor/work/Bio321G_RY_Spring2024/Summer2024/glassfrogs/Cecilia/rsem_ali/GlassFrog_TrinityTranscriptome_27March24.fasta"
#ANNOTATION_GTF="/path/to/annotation_file.gtf"
OUTPUT_DIR="/stor/work/Bio321G_RY_Spring2024/Summer2024/glassfrogs/Cecilia/rsem_ali/rsem_index"
#https://deweylab.github.io/RSEM/rsem-prepare-reference.html

# Create the RSEM reference
# rsem-prepare-reference --gtf $ANNOTATION_GTF $GENOME_FASTA $OUTPUT_DIR/glassfrog_rsem_reference
 rsem-prepare-reference -- $GENOME_FASTA $OUTPUT_DIR/glassfrog_rsem_reference

echo "RSEM reference indexing completed successfully."

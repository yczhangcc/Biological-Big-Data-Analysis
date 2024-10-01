# !/bin/bash
#sudo apt-get update
#sudo apt-get install -y build-essential zlib1g-dev \
#    libncurses5-dev libbz2-dev liblzma-dev \
#    libboost-all-dev cmake git
#Sorry, user yz29548 is not allowed to execute '/usr/bin/apt-get update' as root on educcomp01.ccbb.utexas.edu.
#git clone https://github.com/trinityrnaseq/trinityrnaseq.git
#cd trinityrnaseq
#Trinity --version

# Define variables
TRINITY_DIR="/stor/work/Bio321G_RY_Spring2024/Summer2024/glassfrogs/Cecilia/TRINITY/trinityrnaseq"
TRIMMED_DIR="/stor/work/Bio321G_RY_Spring2024/Summer2024/glassfrogs/trimfix1_results"
OUTPUT_DIR="/stor/work/Bio321G_RY_Spring2024/Summer2024/glassfrogs/trinity_assembly1"

# Check Trinity installation
$TRINITY_DIR/Trinity --version

# Concatenate all R1 and R2 files for input
LEFT_READS=$(ls $TRIMMED_DIR/*_1_paired.fastq | tr '\n' ',')
RIGHT_READS=$(ls $TRIMMED_DIR/*_2_paired.fastq | tr '\n' ',')

# Remove the trailing comma
LEFT_READS=${LEFT_READS%,}
RIGHT_READS=${RIGHT_READS%,}



# Run Trinity with the concatenated FASTQ files
nohup $TRINITY_DIR/Trinity --seqType fq --max_memory 50G --CPU 8 \
    --left "$LEFT_READS" \
    --right "$RIGHT_READS" \
    --output "$OUTPUT_DIR" &

if [ $? -eq 0 ]; then
    echo "Trinity assembly completed successfully."
else
    echo "Trinity assembly failed. Please check the error messages above."
fi



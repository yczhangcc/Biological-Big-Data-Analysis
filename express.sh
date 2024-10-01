# !/bin/bash
#https://pachterlab.github.io/eXpress/manual.html#running
#express -o salmon/ /stor/work/Bio321G_RY_Spring2024/Summer2024/glassfrogs/becca_trinity/trinity_out_dir/GlassFrog_reassembly_11July24.fasta /stor/work/Bio321G_RY_Spring2024/Summer2024/glassfrogs/Cecilia/bwa/newbwa_result/Ep_B-4352_combined_sorted.bam


# Set the paths
output_base_dir="/stor/work/Bio321G_RY_Spring2024/Summer2024/glassfrogs/Cecilia/express/express_result/"
fasta_file="/stor/work/Bio321G_RY_Spring2024/Summer2024/glassfrogs/becca_trinity/trinity_out_dir/GlassFrog_reassembly_11July24.fasta"
bam_dir="/stor/work/Bio321G_RY_Spring2024/Summer2024/glassfrogs/Cecilia/bwa/newbwa_result/"

# Create the base output directory if it doesn't exist
mkdir -p "$output_base_dir"

# Loop through all combined_sorted.bam files in the bam_dir
for bam_file in "$bam_dir"/*combined_sorted.bam; do
  # Extract the base name of the bam file (e.g., "Ep_B-4352_combined_sorted")
  base_name=$(basename "$bam_file" .bam)
  
  # Set the specific output directory for this bam file
  output_dir="${output_base_dir}${base_name}"
  
  # Create the specific output directory if it doesn't exist
  mkdir -p "$output_dir"
  
  # Run eXpress for each BAM file with paired-end parameters
  express -o "$output_dir"  "$fasta_file" "$bam_file"
  
  # Rename the output files by moving them to the specific directory with a new name
  for result_file in "$output_dir"/*; do
    mv "$result_file" "$output_dir/${base_name}_$(basename "$result_file")"
  done
done

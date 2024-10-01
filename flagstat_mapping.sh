#!/bin/bash

for file in /stor/work/Bio321G_RY_Spring2024/Summer2024/glassfrogs/working_dir/*.sorted.bam;do
        echo "samtools flagstat $file >$(basename $file).flagstat"; 

mapping_percentage=$(samtools flagstat your_input_file.bam | grep "mapped (" | awk '{print $5}')
echo "Mapping Percentage: $mapping_percentage";

done
# the higher the percentage the better it is , it will generate result of mapping number as many line text 
#the flagstat and the salmon can be done seperately

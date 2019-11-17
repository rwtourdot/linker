#!/bin/bash

CHR_list=('chr22' 'chr21' 'chr20' 'chr19' 'chr18' 'chr17' 'chr16' 'chr15' 'chr14' 'chr13' 'chr12' 'chr11' 'chr10' 'chr9' 'chr8' 'chr7' 'chr6' 'chr5' 'chr4' 'chr3' 'chr2' 'chr1' 'chrX')

input_het_cov="./output/filtered_coverage_mar12_k562_"
input_hic_bam_file="/czlab/Data/K-562/HiC/Bulk/HiC069.bam"
#input_hic_bam_file="/czlab/Data/K-562/HiC/Bulk/HiC071.bam"

string_id="mar13_K562_69"
nohup_file="nohup_hic_phase_"$string_id"_"

for chr_choice in "${CHR_list[@]}"
do
        loop_nohup_file=$nohup_file$chr_choice".out"
        loop_het_cov_file=$input_het_cov$chr_choice".dat"
        echo $loop_nohup_file
        echo $string_id
	echo $chr_choice
        nohup ./linker hic_phase -c $chr_choice -i $input_hic_bam_file -m $loop_het_cov_file -n $string_id > $loop_nohup_file
done



#!/bin/bash

CHR_list=('chr22' 'chr21' 'chr20' 'chr19' 'chr18' 'chr17' 'chr16' 'chr15' 'chr14' 'chr13' 'chr12' 'chr11' 'chr10' 'chr9' 'chr8' 'chr7' 'chr6' 'chr5' 'chr4' 'chr3' 'chr2' 'chr1' 'chrX')

tech="illumina"
string_id="jan22_RPE1_illumina"  #_filter
nohup_file="nohup_filter_normal_"$string_id"_"$tech"_"
cov_str_normal="jan22_RPE1_illumina"

for chr_choice in "${CHR_list[@]}"
do
        input_normal_cov="./output/het_coverage_"$cov_str_normal"_"$chr_choice".dat"
        loop_nohup_file=$nohup_file$chr_choice".out"
        echo $loop_nohup_file
        echo $string_id"_"$chr_choice
        nohup ./linker filter -m $input_normal_cov -n $string_id"_"$chr_choice > $loop_nohup_file
done


#cov_str_tumor="apr17_tumor_illumina"
#input_tumor_cov="./output/het_coverage_"$cov_str_tumor"_"$chr_choice".dat"

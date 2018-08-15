#!/bin/bash

CHR_list=('chr22' 'chr21' 'chr20' 'chr19' 'chr18' 'chr17' 'chr16' 'chr15' 'chr14' 'chr13' 'chr12' 'chr11' 'chr10' 'chr9' 'chr8' 'chr7' 'chr6' 'chr5' 'chr4' 'chr3' 'chr2' 'chr1' 'chrX')

tech="tenx"
string_id="june18_tumor_EA_tenx"
nohup_file="nohup_filter_tumor_"$string_id"_"$tech"_"
cov_str_normal="june18_normal_EA_tenx"
cov_str_tumor="june18_tumor_EA_tenx"

for chr_choice in "${CHR_list[@]}"
do
        input_tumor_cov="./output/het_coverage_"$cov_str_tumor"_"$chr_choice".dat"
        input_normal_cov="./output/het_coverage_"$cov_str_normal"_"$chr_choice".dat"
        loop_nohup_file=$nohup_file$chr_choice".out"
        echo $loop_nohup_file
        echo $string_id"_"$chr_choice
        nohup ./linker filter -t $input_tumor_cov -m $input_normal_cov -n $string_id"_"$chr_choice > $loop_nohup_file
done




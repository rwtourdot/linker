#!/bin/bash

#CHR_list=('chr1' 'chr2' 'chr3' 'chr4' 'chr5' 'chr6' 'chr7' 'chr8' 'chr9' 'chr10' 'chr11' 'chr12' 'chr13' 'chr14' 'chr15' 'chr16' 'chr17' 'chr18' 'chr19' 'chr20' 'chr21' 'chr22' 'chrX')
CHR_list=('chr22' 'chr21' 'chr20' 'chr19' 'chr18' 'chr17' 'chr16' 'chr15' 'chr14' 'chr13' 'chr12' 'chr11' 'chr10' 'chr9' 'chr8' 'chr7' 'chr6' 'chr5' 'chr4' 'chr3' 'chr2' 'chr1' 'chrX')

binsize=10000
string_id_base="june12_PCRFree_tumor"
nohup_file_base="nohup_cn_phase_"$string_id_base
hap_file_string="hap_solution_may22_tenx"
#cov_file_string="filtered_coverage_june11_PCRFree_normal"
cov_file_string="filtered_coverage_june11_PCRFree_tumor"
#cov_file_string="filtered_coverage_june11_illumina_normal"
#cov_file_string="filtered_coverage_june11_illumina_tumor"

#cov_file_string="filtered_coverage_may23_tumor"
#cov_file_string="filtered_coverage_may23_normal"

for chr_choice in "${CHR_list[@]}"
do
	input_hap_file="./output/"$hap_file_string"_"$chr_choice".dat"
	input_coverage="./output/"$cov_file_string"_"$chr_choice".dat"
        loop_nohup_file=$nohup_file_base"_"$chr_choice".out"
	string_id=$string_id_base"_"$chr_choice
        echo $loop_nohup_file
        echo $string_id
        nohup ./linker cn_phase -s $input_hap_file -d $input_coverage -b $binsize -n $string_id > $loop_nohup_file
done


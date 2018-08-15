#!/bin/bash

CHR_list=('chr22' 'chr21' 'chr20' 'chr19' 'chr18' 'chr17' 'chr16' 'chr15' 'chr14' 'chr13' 'chr12' 'chr11' 'chr10' 'chr9' 'chr8' 'chr7' 'chr6' 'chr5' 'chr4' 'chr3' 'chr2' 'chr1' 'chrX')

input_vcf_file="./sample_data/BL1954_PCRFree.hets.recalibrated.vcf"
#input_bam_file="/czlab/Data/HCC1954_10X_New/GRCh38/BL1954_10xG.bam"
#input_bam_file="/czlab/Data/HCC1954_10X_New/GRCh38/HCC1954_10xG.bam"
input_bam_file="/czlab/Data/HCC1954_10X_New/GRCh38/HCC1954_10xG_EA.bam"
tech="tenx"

string_id="june18_tumor_EA"
nohup_file="nohup_bx_bin_"$string_id"_"$tech"_"
binsize=10000

for chr_choice in "${CHR_list[@]}"
do
        loop_nohup_file=$nohup_file$chr_choice".out"
        echo $loop_nohup_file
        echo $string_id
        nohup ./linker bx_bin -c $chr_choice -i $input_bam_file -e $tech -b $binsize -n $string_id > $loop_nohup_file
done



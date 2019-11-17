#!/bin/bash

CHR_list=('chr22' 'chr21' 'chr20' 'chr19' 'chr18' 'chr17' 'chr16' 'chr15' 'chr14' 'chr13' 'chr12' 'chr11' 'chr10' 'chr9' 'chr8' 'chr7' 'chr6' 'chr5' 'chr4' 'chr3' 'chr2' 'chr1' 'chrX')

input_vcf_file="./sample_data/K562_het_variants.rmcent.vcf"
#input_bam_file="/czlab/Data/K-562-10X/bam1/phased_possorted_bam.bam"
input_bam_file="/czlab/Data/NA12878/10X/NA12878_WGS_v2.bam"
tech="tenx"  #tech="illumina"

string_id="may20_NA12878"
nohup_file="nohup_coverage_"$string_id"_"$tech"_"

for chr_choice in "${CHR_list[@]}"
do
        loop_nohup_file=$nohup_file$chr_choice".out"
        echo $loop_nohup_file
        echo $string_id
        nohup ./linker coverage -c $chr_choice -i $input_bam_file -v $input_vcf_file -e $tech -n $string_id > $loop_nohup_file
done


#input_bam_file="/czlab/Data/HCC1954_10X_New/GRCh38/HCC1954_10xG.bam"
#input_vcf_file="./sample_data/BL1954_PCRFree.hets.recalibrated.vcf"

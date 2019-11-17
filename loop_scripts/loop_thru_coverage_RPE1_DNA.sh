#!/bin/bash

CHR_list=('chr22' 'chr21' 'chr20' 'chr19' 'chr18' 'chr17' 'chr16' 'chr15' 'chr14' 'chr13' 'chr12' 'chr11' 'chr10' 'chr9' 'chr8' 'chr7' 'chr6' 'chr5' 'chr4' 'chr3' 'chr2' 'chr1' 'chrX')

#input_vcf_file="./sample_data/RPE-1.hets.vcf"
#input_bam_file="/singlecellcenter/RPE-1/RPE_bulk.bam"
#input_vcf_file="./sample_data/RPE-1.hets.chr1-X.BA.SNP.vcf"
input_vcf_file="./sample_data/RPE-1.hets.chr1-X.BA.SNP.rmcent.vcf"
input_bam_file="/czlab/Data/RPE-1/bulk_standard/RPE-1_bulk.bam"

tech="illumina"

string_id="jan22_RPE1"
nohup_file="nohup_coverage_"$string_id"_"$tech"_"

for chr_choice in "${CHR_list[@]}"
do
        loop_nohup_file=$nohup_file$chr_choice".out"
        echo $loop_nohup_file
        echo $string_id
        nohup ./linker coverage -c $chr_choice -i $input_bam_file -v $input_vcf_file -e $tech -n $string_id > $loop_nohup_file
done


#input_vcf_file="./sample_data/K562_het_variants.vcf"
#input_bam_file="/czlab/Data/K-562/DNA/K562.standard.bam"

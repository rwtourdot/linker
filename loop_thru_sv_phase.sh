#!/bin/bash

tech="tenx"
string_id="may25_interchromosomal"
#string_id="apr26_intrachromosomal"
nohup_file_base="nohup_sv_phase_"$string_id
bam_file_string="/czlab/Data/HCC1954_10X_New/GRCh38/HCC1954_10xG.bam"
vcf_file_string="./sample_data/BL1954_PCRFree.hets.recalibrated.vcf"
sv_file_string="./sample_data/Interchromosomal.txt"
#sv_file_string="./sample_data/IntraLongRange.txt"

loop_nohup_file=$nohup_file_base".out"
echo $loop_nohup_file
echo $string_id
nohup ./linker sv_phase -i $bam_file_string -v $vcf_file_string -s $sv_file_string -e $tech -n $string_id > $loop_nohup_file







#CHR_list=('chr1' 'chr2' 'chr3' 'chr4' 'chr5' 'chr6' 'chr7' 'chr8' 'chr9' 'chr10' 'chr11' 'chr12' 'chr13' 'chr14' 'chr15' 'chr16' 'chr17' 'chr18' 'chr19' 'chr20' 'chr21' 'chr22' 'chrX')
#CHR_list=('chr22' 'chr21' 'chr20' 'chr19' 'chr18' 'chr17' 'chr16' 'chr15' 'chr14' 'chr13' 'chr12' 'chr11' 'chr10' 'chr9' 'chr8' 'chr7' 'chr6' 'chr5' 'chr4' 'chr3' 'chr2' 'chr1' 'chrX')

#for chr_choice in "${CHR_list[@]}"
#do
#"_"$chr_choice".out"
#"_"$chr_choice

#done


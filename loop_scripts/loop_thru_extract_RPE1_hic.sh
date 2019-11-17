#!/bin/bash

#CHR_list=('chr1' 'chr2' 'chr3' 'chr4' 'chr5' 'chr6' 'chr7' 'chr8' 'chr9' 'chr10' 'chr11' 'chr12' 'chr13' 'chr14' 'chr15' 'chr16' 'chr17' 'chr18' 'chr19' 'chr20' 'chr21' 'chr22' 'chrX')
CHR_list=('chr22' 'chr21' 'chr20' 'chr19' 'chr18' 'chr17' 'chr16' 'chr15' 'chr14' 'chr13' 'chr12' 'chr11' 'chr10' 'chr9' 'chr8' 'chr7' 'chr6' 'chr5' 'chr4' 'chr3' 'chr2' 'chr1' 'chrX' 'chrY')

#input_vcf_file="./sample_data/BL1954_PCRFree.hets.recalibrated.rmcent.vcf"
#input_bam_file="/czlab/Data/HCC1954_10X_New/GRCh38/BL1954_10xG.bam"
input_vcf_file="./vcf_data/RPE-1.hets.chr1-X.BA.SNP.rmcent.vcf"
#input_bam_file="/czlab/Data/RPE-1/bulk_in-situ-Hi-C/RPE1_WT_SRS1045722_in_situ_Hi-C.bam"
input_bam_file="/czlab/Data/RPE-1/bulk_in-situ-Hi-C/RPE1_WT_SRS1045724_in_situ_Hi-C.bam"
#input_bam_file2="/czlab/Data/HCC1954_10X_New/GRCh38/HCC1954_10xG.bam"
#input_vcf_file="/czlab/Data/NA12878/10X/NA12878_WGS_v2.vcf"
#input_bam_file="/czlab/Data/NA12878/10X/NA12878_WGS_v2.bam"
#input_bam_file2="/czlab/Data/NA12878/10X/NA12878_WGS_210.bam"
tech="hic"

string_id="nov7_RPE1_24"
nohup_file="nohup_extract_"$string_id"_"$tech"_"


for chr_choice in "${CHR_list[@]}"
do
        loop_nohup_file=$nohup_file$chr_choice".out"
        echo $loop_nohup_file
        echo $string_id
        #nohup ./linker phase -c $chr_choice -i $input_bam_file -v $input_vcf_file -e $tech -n $string_id > $loop_nohup_file &
        nohup ./linker extract -c $chr_choice -i $input_bam_file -v $input_vcf_file -e $tech -n $string_id > $loop_nohup_file &
        #nohup ./linker extract -c $chr_choice -i $input_bam_file -s $input_bam_file2 -v $input_vcf_file -e $tech -n $string_id > $loop_nohup_file &
done



#input_vcf_file="./sample_data/BL1954_PCRFree.hets.recalibrated.vcf"

#!/bin/bash

#CHR_list=('chr1' 'chr2' 'chr3' 'chr4' 'chr5' 'chr6' 'chr7' 'chr8' 'chr9' 'chr10' 'chr11' 'chr12' 'chr13' 'chr14' 'chr15' 'chr16' 'chr17' 'chr18' 'chr19' 'chr20' 'chr21' 'chr22' 'chrX')
#CHR_list=('chr22' 'chr21' 'chr20' 'chr19' 'chr18' 'chr17' 'chr16' 'chr15' 'chr14' 'chr13' 'chr12' 'chr11' 'chr10' 'chr9' 'chr8' 'chr7' 'chr6' 'chr5' 'chr4' 'chr3' 'chr2' 'chr1' 'chrX')

#input_graph_file="./output/graph_variant_jun10_BL1954_tenx_"
#input_graph_file="./output/graph_variant_oct23_BL1954_all_hic_"
#input_hap_file="./output/hap_solution_oct21_HCC1954_"
#input_hap_file="./output/hap_solution_oct24_BL1954_"
#string_id="nov5_BL1954_700"
#nohup_file="nohup_scaffold_"$string_id"_"

#block_cutoff=-1800
#block_cutoff=-1700
#block_cutoff=-2000
#block_cutoff=-700

cn_file_path="./../cpp_linker_v4.0_somatic_cn_phasing/allelic_counts_tumor_svaba__plotACN_pdf"
cn_files=$( ls $cn_file_path/EAC-*.CN.tsv )

binsize=25000

for cfile in $cn_files
do
	sample_id="$(echo $cfile | cut -d\/ -f5 | cut -d\- -f2,3 | cut -d\. -f1 )"
	chrom="$(echo $cfile | cut -d\/ -f5 | cut -d\. -f3 )"
	loop_label=$sample_id"_"$chrom"_nov14"
	loop_nohup_file="nohup_cn_phase_"$loop_label".out"
	echo $sample_id
	echo $loop_label
	echo $loop_nohup_file
        nohup ./linker cn_phase -m $cfile -n $loop_label -b $binsize > $loop_nohup_file 
done




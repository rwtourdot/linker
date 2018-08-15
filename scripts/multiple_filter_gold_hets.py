import os
import sys

chrom_list = ["chr1","chr2","chr3","chr4","chr5","chr6","chr7","chr8","chr9","chr10","chr11","chr12","chr13","chr14","chr15","chr16","chr17","chr18","chr19","chr20","chr21","chr22","chrX"]

##########################################################
het_file_folder = "./../output"
het_file_tenx = "filtered_coverage_may23_normal_tenx_"
het_file_pcrfree = "filtered_coverage_june11_PCRFree_normal_"
het_file_illumina = "filtered_coverage_june11_illumina_normal_"

het_file_tenx_tumor = "filtered_coverage_may23_tumor_tenx_"
het_file_pcrfree_tumor = "filtered_coverage_june11_PCRFree_tumor_"
het_file_illumina_tumor = "filtered_coverage_june11_illumina_tumor_"

output_het_file_tenx = "gold_" + het_file_tenx
output_het_file_pcrfree = "gold_" + het_file_pcrfree
output_het_file_illumina = "gold_" + het_file_illumina

output_het_file_tenx_tumor = "gold_" + het_file_tenx_tumor
output_het_file_pcrfree_tumor = "gold_" + het_file_pcrfree_tumor
output_het_file_illumina_tumor = "gold_" + het_file_illumina_tumor

##########################################################
def load_het_cov(hfile,chr,chr_het_dict):
    lines = [line.split() for line in open(hfile)]
    chr_het_dict[chr] = [];
    for l in lines:
        position = int(l[1]);
        chr_het_dict[chr].append(position);
    return chr_het_dict

def filter_het_union(hfile,chr,udict,ofile):
    fid = open(ofile,'w') 
    lines = [line.split() for line in open(hfile)]
    for l in lines:
        position = int(l[1]);
        if position in udict:
            fid.write("\t".join(l)+"\n")   #print("\t".join(l))
    fid.close()


##########################################################
chr_het_dict_tenx = {}; chr_het_dict_illumina = {}; chr_het_dict_pcrfree = {}
for chr in chrom_list:
    het_tenx = het_file_folder + "/" + het_file_tenx + chr + ".dat"
    het_pcrfree = het_file_folder + "/" + het_file_pcrfree + chr + ".dat"
    het_illumina = het_file_folder + "/" + het_file_illumina + chr + ".dat"
    chr_het_dict_tenx = load_het_cov(het_tenx,chr,chr_het_dict_tenx)
    chr_het_dict_illumina = load_het_cov(het_illumina,chr,chr_het_dict_illumina)
    chr_het_dict_pcrfree = load_het_cov(het_pcrfree,chr,chr_het_dict_pcrfree)


##########################################################
union_chr_dict = {}
for chr in chrom_list:
    set_tenx_chr = set(chr_het_dict_tenx[chr])
    set_illumina_chr = set(chr_het_dict_illumina[chr])
    set_pcrfree_chr = set(chr_het_dict_pcrfree[chr])
    union_pos_chr = list(set_tenx_chr & set_illumina_chr & set_pcrfree_chr)
    union_chr_dict[chr] = union_pos_chr
    print(chr,len(union_pos_chr))


for chr in chrom_list:
    #het_tenx = het_file_folder + "/" + het_file_tenx + chr + ".dat"
    #het_pcrfree = het_file_folder + "/" + het_file_pcrfree + chr + ".dat"
    #het_illumina = het_file_folder + "/" + het_file_illumina + chr + ".dat"
    #ofile_het_tenx = het_file_folder + "/" + output_het_file_tenx + chr + ".dat"
    #ofile_het_pcrfree = het_file_folder + "/" + output_het_file_pcrfree + chr + ".dat"
    #ofile_het_illumina = het_file_folder + "/" + output_het_file_illumina + chr + ".dat"
    het_tenx = het_file_folder + "/" + het_file_tenx_tumor + chr + ".dat"
    het_pcrfree = het_file_folder + "/" + het_file_pcrfree_tumor + chr + ".dat"
    het_illumina = het_file_folder + "/" + het_file_illumina_tumor + chr + ".dat"
    ofile_het_tenx = het_file_folder + "/" + output_het_file_tenx_tumor + chr + ".dat"
    ofile_het_pcrfree = het_file_folder + "/" + output_het_file_pcrfree_tumor + chr + ".dat"
    ofile_het_illumina = het_file_folder + "/" + output_het_file_illumina_tumor + chr + ".dat"
    filter_het_union(het_tenx,chr,union_chr_dict[chr],ofile_het_tenx)
    filter_het_union(het_pcrfree,chr,union_chr_dict[chr],ofile_het_pcrfree)
    filter_het_union(het_illumina,chr,union_chr_dict[chr],ofile_het_illumina)






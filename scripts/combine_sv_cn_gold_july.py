import os
import sys

chrom_list = ["chr1","chr2","chr3","chr4","chr5","chr6","chr7","chr8","chr9","chr10","chr11","chr12","chr13","chr14","chr15","chr16","chr17","chr18","chr19","chr20","chr21","chr22","chrX"]
cn_phased_folder = "./../output"
sv_file_long = "./../output/sv_phased_july18_long.dat"
sv_file_short = "./../output/sv_phased_july18_short.dat"
##########################################################
cn_file_string = "hap_gold_solution_july_"
sv_output_file = "./../output/phased_sv_sites_gold_july.dat"


##########################################################
class sv_junction:
    def __init__(self,chr1,chr2,pos1,pos2,ty):
        self.chr1 = chr1;
        self.chr2 = chr2;
        self.pos1 = pos1;
        self.pos2 = pos2;
        self.haplotype1 = {};
        self.haplotype2 = {};
        self.type_sv = ty;

    def add_het_site1(self,position,hap):
        self.haplotype1[position] = hap;

    def add_het_site2(self,position,hap):
        self.haplotype2[position] = hap;


##########################################################
def load_sv_file(sv_file,sv_list,type_sv):
    lines = [line.split() for line in open(sv_file)]
    for l in lines:
        #print(l)
        if len(l) == 2:
            chr1 = l[0].split(":")[0]; pos1 = l[0].split(":")[1];
            chr2 = l[1].split(":")[0]; pos2 = l[1].split(":")[1];
            sv_list.append(sv_junction(chr1,chr2,pos1,pos2,type_sv)); #i += 1;
        elif len(l) > 1:   #print i,l
            if int(l[0]) == 0:
                sv_list[-1].add_het_site1(int(l[2]),int(l[3]))
            elif int(l[0]) == 1:
                sv_list[-1].add_het_site2(int(l[2]),int(l[3]))
    return sv_list


##########################################################
def load_cn_phase(coverage_file,chr,chr_dict):
    lines = [line.split() for line in open(hap_solution_file)]
    chr_hap_dict = {};  #print(l)
    for l in lines:   
	#position = int(l[2]); haplotype = int(l[5]);
        position = int(l[1]); haplotype = int(l[4]);
        chr_hap_dict[position] = haplotype;   #print chr,position,haplotype
    chr_dict[chr] = chr_hap_dict;
    return chr_dict


##########################################################
chr_dict = {}
for chr in chrom_list:
    hap_solution_file = cn_phased_folder + "/" + cn_file_string + chr + ".dat"
    chr_dict = load_cn_phase(hap_solution_file,chr,chr_dict)


##########################################################
sv_list = [];
sv_list = load_sv_file(sv_file_long,sv_list,"long")
sv_list = load_sv_file(sv_file_short,sv_list,"short")

##########################################################
fout = open(sv_output_file, 'w')
for sv in sv_list:
    end1 = 0;  sum1 = 0;  end2 = 0;  sum2 = 0;
    type_sv = sv.type_sv
    if sv.chr1 in chr_dict:
        for key,val in sv.haplotype1.items():     #.iteritems():
            if key in chr_dict[sv.chr1]:
                end1 += 1;  sum1 += val*chr_dict[sv.chr1][key]
    if sv.chr2 in chr_dict:
        for key,val in sv.haplotype2.items():     #.iteritems():
            if key in chr_dict[sv.chr2]:
                end2 += 1;  sum2 += val*chr_dict[sv.chr2][key]
    frac1 = round(float(sum1)/float(end1)) if end1 != 0 else 0.0
    frac2 = round(float(sum2)/float(end2)) if end2 != 0 else 0.0
    fout.write(str(sv.chr1)+"\t"+str(sv.pos1)+"\t"+str(sv.chr2)+"\t"+str(sv.pos2)+"\t"+str(frac1)+"\t"+str(frac2)+"\t"+str(type_sv)+"\n" )






































    #print sv.chr1,sv.pos1,sv.chr2,sv.pos2,len(sv.haplotype1),len(sv.haplotype2)
                #print key,val,chr_dict[sv.chr1][key]
                #print key,val,chr_dict[sv.chr2][key]

##########################################################
#variant_dict,illumina_coverage = load_illumina_coverage(illumina_coverage_file)

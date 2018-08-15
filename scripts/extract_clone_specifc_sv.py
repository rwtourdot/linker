import os
import sys

inter = "./../sample_data/Interchromosomal.txt"
intra = "./../sample_data/IntraLongRange.txt"
fldbk = "./../sample_data/FoldBackBkps.txt"

inter_out = "./../sample_data/inter_normal.dat"

def load_sv_file(sv_file):
    lines = [line.split() for line in open(sv_file)]
    lines.pop(0)
    first_line = lines[0]
    for l in lines:
         chr1,pos1,dir1,chr2,pos2,dir2,HCC1954BL_Illumina,HCC1954BL_PCRFree,HCC1954BL_PacBio,HCC1954_Illumina,HCC1954_PCRFree,HCC1954_PacBio = l[1],int(l[2]),int(l[3]),l[4],int(l[5]),int(l[6]),int(l[14]),int(l[15]),int(l[16]),int(l[17]),int(l[18]),int(l[19]) 
         #print(chr1,pos1,chr2,pos2,HCC1954BL_Illumina,HCC1954BL_PCRFree,HCC1954BL_PacBio,HCC1954_Illumina,HCC1954_PCRFree,HCC1954_PacBio)
         #if (HCC1954BL_Illumina > 10 and HCC1954BL_PCRFree > 10 and HCC1954BL_PacBio > 5):
         #     print("normal SV")
         if (HCC1954_Illumina > 10 and HCC1954_PCRFree > 10 and HCC1954_PacBio > 5):
               print("SV")
               print(chr1,pos1,chr2,pos2,HCC1954BL_Illumina,HCC1954BL_PCRFree,HCC1954BL_PacBio,HCC1954_Illumina,HCC1954_PCRFree,HCC1954_PacBio)

load_sv_file(intra)
#load_sv_file(fldbk)

#print(first_line)




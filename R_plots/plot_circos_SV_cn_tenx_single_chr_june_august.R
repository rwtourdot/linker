# script to read and parse data files
rm(list=ls())
#library(tidyverse)
library(ggplot2)
library(readr)
library(tidyr)
library(dplyr)
library(grid)
library(ggplot2)
library(circlize)
require(gridExtra)

#############################################################
#setwd("~/2018_7_july_workdir/chosen_phase3/")
setwd("~/2018_7_july_workdir/hap_gold_plots/")
source("load_fun.R")

#############################################################
#chrom_list <- c("chr1","chr2","chr3","chr4","chr5","chr6","chr7","chr8","chr9","chr10","chr11","chr12","chr13","chr14","chr15","chr16","chr17","chr18","chr19","chr20","chr21","chr22","chrX","chrY")
chrom_list <- c("chrX")

#############################################################
hap_file <- "./hap_gold_solution_july_"
binsize <- 10000
band_width <- 10000
cutoff_cov <- 5
#binsize_str <- "10000" 

#############################################################
het_file <- "./../phasing_HCC1954_linker_v1.8_illumina_PCRFree_tenx/output/gold_filtered_coverage_may23_tumor_tenx_"
het_file_normal <- "./../phasing_HCC1954_linker_v1.8_illumina_PCRFree_tenx/output/gold_filtered_coverage_may23_normal_tenx_"
het_file_bx <- "./../phasing_HCC1954_linker_v1.8_illumina_PCRFree_tenx/output/het_coverage_may22_tumor_tenx_"
bx_coverage_file <- "./../phasing_HCC1954_linker_v1.8_illumina_PCRFree_tenx/output/bx_unique_binsize_10000_june8_tumor_"
bx_coverage_normal_file <- "./../phasing_HCC1954_linker_v1.8_illumina_PCRFree_tenx/output/bx_unique_binsize_10000_june8_normal_"
one_copy <- 0.17 # tenx

#het_file <- "./../phasing_HCC1954_linker_v1.8_illumina_PCRFree_tenx/output/gold_filtered_coverage_june18_tumor_EA_tenx_"
#het_file_normal <- "./../phasing_HCC1954_linker_v1.8_illumina_PCRFree_tenx/output/gold_filtered_coverage_june18_normal_EA_tenx_"
#het_file_bx <- "./../phasing_HCC1954_linker_v1.8_illumina_PCRFree_tenx/output/het_coverage_june18_tumor_EA_tenx_"
#bx_coverage_file <- "./../phasing_HCC1954_linker_v1.8_illumina_PCRFree_tenx/output/bx_unique_binsize_10000_june18_tumor_EA_"
#bx_coverage_normal_file <- "./../phasing_HCC1954_linker_v1.8_illumina_PCRFree_tenx/output/bx_unique_binsize_10000_june18_normal_EA_"
#one_copy <- 0.24 # tenx_EA

#############################################################
hap_tracks <- load_hap_circos_tenx(hap_file,binsize,chrom_list,one_copy)
hapA <- hap_tracks %>% select(chr,min_pos.x,max_pos.x,value1=tumorAnorm_unique)
hapB <- hap_tracks %>% select(chr,min_pos.x,max_pos.x,value1=tumorBnorm_unique)  

#############################################################
phased_file <- "./phased_sv_sites_gold_july.dat"
#sv_out <- load_phasedsv_file_total(phased_file,chrom_list,band_width)
sv_out <- load_phasedsv_file_total(phased_file,chrom_list,band_width)
anchor <- load_phasedsv_positions(phased_file,chrom_list,band_width)
anchor <- anchor %>% mutate(max_pos = as.integer(min_pos + 1))
anchor <- anchor %>% select(chr,min_pos,max_pos,value1)
anchor <- anchor %>% rename("value2"="value1") %>% mutate(value1=1) %>% mutate(value2=value2+1)

#############################################################
sv_out$m_link1_long <- sv_out$m_link1 %>% filter(call_type == "long") %>% select(chr,start,end)
sv_out$m_link2_long <- sv_out$m_link2 %>% filter(call_type == "long") %>% select(chr,start,end)
sv_out$a_link1_long <- sv_out$a_link1 %>% filter(call_type == "long") %>% select(chr,start,end)
sv_out$a_link2_long <- sv_out$a_link2 %>% filter(call_type == "long") %>% select(chr,start,end)
sv_out$b_link1_long <- sv_out$b_link1 %>% filter(call_type == "long") %>% select(chr,start,end)
sv_out$b_link2_long <- sv_out$b_link2 %>% filter(call_type == "long") %>% select(chr,start,end)
sv_out$n_link1_long <- sv_out$n_link1 %>% filter(call_type == "long") %>% select(chr,start,end)
sv_out$n_link2_long <- sv_out$n_link2 %>% filter(call_type == "long") %>% select(chr,start,end)

sv_out$m_link1_short <- sv_out$m_link1 %>% filter(call_type == "short") %>% select(chr,start,end)
sv_out$m_link2_short <- sv_out$m_link2 %>% filter(call_type == "short") %>% select(chr,start,end)
sv_out$a_link1_short <- sv_out$a_link1 %>% filter(call_type == "short") %>% select(chr,start,end)
sv_out$a_link2_short <- sv_out$a_link2 %>% filter(call_type == "short") %>% select(chr,start,end)
sv_out$b_link1_short <- sv_out$b_link1 %>% filter(call_type == "short") %>% select(chr,start,end)
sv_out$b_link2_short <- sv_out$b_link2 %>% filter(call_type == "short") %>% select(chr,start,end)
sv_out$n_link1_short <- sv_out$n_link1 %>% filter(call_type == "short") %>% select(chr,start,end)
sv_out$n_link2_short <- sv_out$n_link2 %>% filter(call_type == "short") %>% select(chr,start,end)

#############################################################
log_scale_break <- 2
addition_num <- 1
hapA <- hapA %>% mutate( ln2_mean = case_when( value1 > log_scale_break ~ log2(value1) + addition_num, value1 <= log_scale_break ~ value1 ) )
hapB <- hapB %>% mutate( ln2_mean = case_when( value1 > log_scale_break ~ log2(value1) + addition_num, value1 <= log_scale_break ~ value1 ) )
hapA <- hapA %>% select(chr,min_pos.x,max_pos.x,ln2_mean) %>% rename("value1"="ln2_mean")
hapB <- hapB %>% select(chr,min_pos.x,max_pos.x,ln2_mean) %>% rename("value1"="ln2_mean")
combined_hap <- list(hapA,hapB)

#############################################################  #FALSE
circos.par( clock.wise = FALSE, cell.padding = c(0, 0, 0, 0), gap.degree = 180)  #180   #90       #"canvas.xlim" = c(-2, 2), "canvas.ylim" = c(-1, 1),
circos.initializeWithIdeogram( species="hg38", chromosome.index = chrom_list, ideogram.height = 0.01)  #plotType=NULL

#############################################################    #, cell.padding = c(0.0, 0, 0.0, 0)   900  , 800,    ### 40
circos.genomicTrackPlotRegion(combined_hap, ylim = c(0, 9), track.margin	= c(0.0, 0.0), cell.padding = c(0.0, 0, 0.0, 0), bg.border = "white", track.height = 0.4, panel.fun = function(region, value, ...) {
  index = getI(...)
  col = ifelse(index > 1, "deepskyblue3", "brown3")
  circos.genomicPoints(region, value, col = col, cex = 0.3, pch = 16) # index    #cex is size   #16   12 , log10(value)
  cell.xlim = get.cell.meta.data("cell.xlim")
  for( h in c( 0, 1, 2, 3, 4, 5, 6, 7, 8) ) { circos.lines(cell.xlim, c(h, h), col = "#00000040") }   #, 300, 400, 500
})

#############################################################    #c(-1, 1)
foldback <- bind_rows(sv_out$n_link1_short,sv_out$n_link2_short,sv_out$m_link1_short,sv_out$m_link2_short,sv_out$a_link1_short,sv_out$a_link2_short,sv_out$b_link1_short,sv_out$b_link2_short)
foldback <- foldback %>% mutate(value=1) 
circos.genomicTrackPlotRegion(foldback, ylim = c(0, 1), bg.border = NA, track.height = 0.04, panel.fun = function(region, value, ...) {
  circos.genomicLines(region, value, type = "h", col = "orange3", pch = 4, cex = 0.1, lwd=1.5)
})

#############################################################    #c(-1, 1)
circos.genomicTrackPlotRegion(anchor, ylim = c(0, 1), bg.border = NA, track.height = 0.04, panel.fun = function(region, value, ...) {
  circos.genomicLines(region, value[2], type = "h", col = ifelse(value[1] > 1, "brown3", ifelse(value[1] < 1,"deepskyblue3","grey50")), pch = 4, cex = 0.1, lwd=1.5)   #cex = 0.3  #pch = 16
})

############################################### ##############
circos.genomicLink(sv_out$n_link1_long,sv_out$n_link2_long,col="grey50",lwd=1.5)
circos.genomicLink(sv_out$m_link1_long,sv_out$m_link2_long,col="darkorchid3",lwd=1.5)   # lwd = 1.5
circos.genomicLink(sv_out$a_link1_long,sv_out$a_link2_long,col="brown3",lwd=1.5) #, col="#3182bd"
circos.genomicLink(sv_out$b_link1_long,sv_out$b_link2_long,col="deepskyblue3",lwd=1.5)

#circos.genomicLink(sv_out$n_link1_short,sv_out$n_link2_short,col="orange3",lwd=1.5)
#circos.genomicLink(sv_out$m_link1_short,sv_out$m_link2_short,col="orange3",lwd=1.5)
#circos.genomicLink(sv_out$a_link1_short,sv_out$a_link2_short,col="orange3",lwd=1.5) #, col="#3182bd"
#circos.genomicLink(sv_out$b_link1_short,sv_out$b_link2_short,col="orange3",lwd=1.5)

#highlight.chromosome("chr21")
#highlight.chromosome("chr5")
#highlight.chromosome("chr8")
#highlight.chromosome("chr11")
circos.clear()






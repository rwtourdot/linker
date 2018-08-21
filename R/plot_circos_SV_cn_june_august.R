# script to read and parse data files
rm(list=ls())
#library(tidyverse)
library(ggplot2)
library(readr)
library(tidyr)
library(dplyr)
library(grid)
library(circlize)
require(gridExtra)


#############################################################
#setwd("~/2018_7_july_workdir/chosen_phase3/")
setwd("~/2018_7_july_workdir/hap_gold_plots/")
source("load_fun.R")

#############################################################
#chrom_list <- c("chr1","chr2","chr3","chr4","chr5","chr6","chr7","chr8","chr9","chr10","chr11","chr12","chr13","chr14","chr15","chr16","chr17","chr18","chr19","chr20","chr21","chr22","chrX","chrY")
chrom_list <- c("chr1","chr2","chr3","chr4","chr5","chr6","chr7","chr8","chr9","chr10","chr11","chr12","chr13","chr14","chr15","chr16","chr17","chr18","chr19","chr20","chr21","chr22","chrX")

#chrom_list <- c("chr4","chr5","chr6","chr7","chr8","chr9","chr10","chr11","chr12","chr14","chr15","chr16","chr17","chr18","chr19","chr20","chr21","chr22")
#chrom_list <- c("chr4","chr5","chr6","chr8","chr9","chr10","chr11","chr12","chr14","chr17","chr19","chr20","chr21","chr22")
#chrom_list <- c("chr4","chr5","chr8","chr9","chr20","chr21","chr22")  #"chr8"

#chrom_list <- c("chr5","chr8","chr11","chr21")
#chrom_list <- c("chr5","chr8","chr11","chr17","chr21")
#chrom_list <- c("chr5","chr8","chr9","chr11","chr21")
#chrom_list <- c("chr1","chr20")

#############################################################
hap_file <- "./hap_gold_solution_july_"
binsize <- 10000
band_width <- 10000
cutoff_cov <- 5
#binsize_str <- "10000" 

#############################################################
#het_file <- "./../phasing_HCC1954_linker_v1.8_illumina_PCRFree_tenx/output/gold_filtered_coverage_may23_tumor_tenx_"
#het_file_normal <- "./../phasing_HCC1954_linker_v1.8_illumina_PCRFree_tenx/output/gold_filtered_coverage_may23_normal_tenx_"
#one_copy <- 0.21  # tenx
#het_file <- "./../phasing_HCC1954_linker_v1.8_illumina_PCRFree_tenx/output/gold_filtered_coverage_june18_tumor_EA_tenx_"
#het_file_normal <- "./../phasing_HCC1954_linker_v1.8_illumina_PCRFree_tenx/output/gold_filtered_coverage_june18_normal_EA_tenx_"
#one_copy <- 0.21 # tenx_EA
het_file <- "./../phasing_HCC1954_linker_v1.8_illumina_PCRFree_tenx/output/gold_filtered_coverage_june11_PCRFree_tumor_"
het_file_normal <- "./../phasing_HCC1954_linker_v1.8_illumina_PCRFree_tenx/output/gold_filtered_coverage_june11_PCRFree_normal_"
one_copy <- 0.165 # PCRFree
#het_file <- "./../phasing_HCC1954_linker_v1.8_illumina_PCRFree_tenx/output/gold_filtered_coverage_june11_illumina_tumor_"
#het_file_normal <- "./../phasing_HCC1954_linker_v1.8_illumina_PCRFree_tenx/output/gold_filtered_coverage_june11_illumina_normal_"
#one_copy <- 0.19 # illumina

#############################################################
hap_tracks <- load_hap_circos(hap_file,binsize,chrom_list,one_copy)
hapA <- hap_tracks %>% select(chr,min_pos,max_pos,value1=tumorAnorm)
hapB <- hap_tracks %>% select(chr,min_pos,max_pos,value1=tumorBnorm) 

#############################################################
phased_file <- "./phased_sv_sites_gold_july.dat"
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
hapA <- hapA %>% select(chr,min_pos,max_pos,ln2_mean) %>% rename("value1"="ln2_mean")
hapB <- hapB %>% select(chr,min_pos,max_pos,ln2_mean) %>% rename("value1"="ln2_mean")
combined_hap <- list(hapA,hapB)

#breaks=c(0,1,2,3,4,5,6,7,8,9),label=c("0","1","2","4","8","16","32","64","128","256")

#combined_hap <- combined_hap %>% select(chr,min_pos,max_pos,ln2_mean) %>% rename("value1"=="ln2_mean")
#circos.par("cell.padding" = c(0, 0, 0, 0),"start.degree" = 90)
circos.initializeWithIdeogram(species="hg38",chromosome.index = chrom_list,ideogram.height = 0.01)

#############################################################    #, cell.padding = c(0.0, 0, 0.0, 0)   900  , 800,    ### 40
circos.genomicTrackPlotRegion(combined_hap, ylim = c(0, 9), track.margin	= c(0.0, 0.0), cell.padding = c(0.0, 0, 0.0, 0), bg.border = "white", track.height = 0.3, panel.fun = function(region, value, ...) {
  index = getI(...)
  col = ifelse(index > 1, "deepskyblue3", "brown3")
  circos.genomicPoints(region, value, col = col, cex = 0.2, pch = 16) # index    #cex is size   #16   12 , log10(value)
  cell.xlim = get.cell.meta.data("cell.xlim")
  for( h in c(0,1,2,3,4,5,6,7,8) ) { circos.lines(cell.xlim, c(h, h), col = "#00000040") }
})

#############################################################    #c(-1, 1)
foldback <- bind_rows(sv_out$n_link1_short,sv_out$n_link2_short,sv_out$m_link1_short,sv_out$m_link2_short,sv_out$a_link1_short,sv_out$a_link2_short,sv_out$b_link1_short,sv_out$b_link2_short)
foldback <- foldback %>% mutate(value=1) 
circos.genomicTrackPlotRegion(foldback, ylim = c(0.3, 0.8), bg.border = NA, track.height = 0.05, panel.fun = function(region, value, ...) {
  circos.genomicLines(region, value, type = "h", col = "orange3", pch = 4, cex = 0.1, lwd=1.5)
})

#############################################################    #c(-1, 1)
circos.genomicTrackPlotRegion(anchor, ylim = c(0.3, 0.8), bg.border = NA, track.height = 0.05, panel.fun = function(region, value, ...) {
  circos.genomicLines(region, value[2], type = "h", col = ifelse(value[1] > 1, "brown3", ifelse(value[1] < 1,"deepskyblue3","grey50")), pch = 4, cex = 0.1, lwd=1.5)   #cex = 0.3  #pch = 16
})

#############################################################
#circos.genomicLink(sv_out$n_link1_long,sv_out$n_link2_long,col="grey50",lwd=1.5)
circos.genomicLink(sv_out$m_link1_long,sv_out$m_link2_long,col="darkorchid3",lwd=1.5)   # lwd = 1.5
circos.genomicLink(sv_out$a_link1_long,sv_out$a_link2_long,col="brown3",lwd=1.5) #, col="#3182bd"
circos.genomicLink(sv_out$b_link1_long,sv_out$b_link2_long,col="deepskyblue3",lwd=1.5)

#circos.genomicLink(sv_out$n_link1_short,sv_out$n_link2_short,col="orange3",lwd=1.5)  #,transparency = 0.5
#circos.genomicLink(sv_out$m_link1_short,sv_out$m_link2_short,col="orange3",lwd=1.5)
#circos.genomicLink(sv_out$a_link1_short,sv_out$a_link2_short,col="orange3",lwd=1.5) #, col="#3182bd"
#circos.genomicLink(sv_out$b_link1_short,sv_out$b_link2_short,col="orange3",lwd=1.5)

#highlight.chromosome("chr21")
#highlight.chromosome("chr5")
#highlight.chromosome("chr8")
#highlight.chromosome("chr11")

circos.clear()









#  circos.genomicLines(region, value[2], type = "h", col = ifelse(value[1] > 1, "brown3", ifelse(value[1] < 1,"deepskyblue3","grey50")), pch = 16, cex = 0.3, lwd=1.5)




#hapA_high = hapA %>% filter(value1 > cutoff_cov)
#hapB_high = hapB %>% filter(value1 > cutoff_cov)
#hapA_low = hapA %>% filter(value1 < cutoff_cov)
#hapB_low = hapB %>% filter(value1 < cutoff_cov)
#combined_hap_high <- list(hapA_high,hapB_high)
#combined_hap_low <- list(hapA_low,hapB_low)




#############################################################    #, cell.padding = c(0.0, 0, 0.0, 0)   900  , 800,    ### 40
#circos.genomicTrackPlotRegion(combined_hap_high, ylim = c(cutoff_cov, 50), track.margin	= c(0.0, 0.0), cell.padding = c(0.0, 0, 0.0, 0), bg.border = "white", track.height = 0.1, panel.fun = function(region, value, ...) {
#  index = getI(...)
#  col = ifelse(index > 1, "deepskyblue3", "brown3")
#  circos.genomicPoints(region, value, col = col, cex = 0.2, pch = 16) # index    #cex is size   #16   12 , log10(value)
#  cell.xlim = get.cell.meta.data("cell.xlim")
#  for( h in c( 10, 20, 30, 40) ) { circos.lines(cell.xlim, c(h, h), col = "#00000040") }   #, 300, 400, 500
#})

#############################################################  #1.0
#circos.genomicTrackPlotRegion(combined_hap_low, ylim = c(0, cutoff_cov), track.margin	= c(0.0, 0.0), cell.padding = c(0.0, 0, 0.0, 0), bg.border = "white", track.height = 0.15, panel.fun = function(region, value, ...) {
#  index = getI(...)
#  col = ifelse(index > 1, "deepskyblue3", "brown3")
#  circos.genomicPoints(region, value, col = col, cex = 0.2, pch = 16) # index    #cex is size   #16   12 , log10(value)    #cex = 0.25
#  cell.xlim = get.cell.meta.data("cell.xlim")
#  for( h in c( 0, 1, 2, 3, 4, 5, 6 ) ) { circos.lines(cell.xlim, c(h, h), col = "#00000040") }   #, 300, 400, 500
#})














#log_scale_break <- 2
#addition_num <- 1
#binned <- combined_hap_high %>% mutate( ln2_meanA = case_when( value1 > log_scale_break ~ log2(value1) + addition_num, value1 <= log_scale_break ~ value1 ) ) # , ln2_meanB = case_when( tumorBnorm > log_scale_break ~ log2(tumorBnorm) + addition_num, tumorBnorm <= log_scale_break ~ tumorBnorm ) )



#chrom_list <- c("chr3","chr4","chr20")
#chrom_list <- c("chr4","chr5","chr8","chr9","chr11","chr12")
#chrom_list <- c("chr1","chr2","chr5","chr8","chr3","chr9")
#chrom_list <- c("chr1","chr4","chr11","chr17")


##############################################################
#hdata <- load_hap(hap_file,binsize,chromosome)
#blocks <- hdata %>% group_by(block) %>% summarise(min_pos = min(pos),max_pos = max(pos))
#binned <- load_het_cov(het_file,het_file_normal,one_copy,chromosome,hdata)
#binned <- binned %>% filter(n>5)

#############################################################
#hap_tracks = tibble() 
#for (chr in chrom_list) { 
#  hdata <- load_hap(hap_file,binsize,chr)
#  blocks <- hdata %>% group_by(block) %>% summarise(min_pos = min(pos),max_pos = max(pos))
#  binned <- load_het_cov(het_file,het_file_normal,one_copy,chr,hdata)
#  binned <- binned %>% filter(n>5)
#  binned <- binned %>% mutate("chr" = chr);
#  hap_tracks <- bind_rows(hap_tracks,binned)
#}



#circos.genomicLines(region, value[2], type = "h", col = ifelse(value[1] > 1, "red", "blue"), pch = 16, cex = 0.3)


#print("hello")
#print(region)
#print(value[1]$value2)
#print(value[2])
#print(value[1])
#print(data.frame(color[value[1]$value2]))
#df2 <- data.frame(color[value[1]$value2])
#col = ifelse(data.frame(value[1]$value2) > 0, "red", "blue")
#col = ifelse(value > 0, "red", "blue")
#col = data.frame(ifelse(value[1] > 1, "red", "blue"))
#ifelse(value[1] > 1, "red", ifelse(value[1] == 1,"gray","gray"))







#bed = generateRandomBed(nr = 5)
#bed = generateRandomBed(nr = 10)








#circos.genomicPoints(region, 0, col = "red", pch = 16, cex = 0.3)
#circos.genomicLines(region, value, lty = "segment")  #lwd=2,
#circos.genomicRect(region, value, ytop = 1, ybottom = -1, border = NA)  #, posTransform = posTransform.default
#circos.genomicRect(region, ytop = 1, ybottom = -1, col = "orange", border = NA)

#bed = generateRandomBed(nr = 200)  #anchor   bed


#circos.track(ylim = c(0, 1), panel.fun = function(x, y) {
#  chr = CELL_META$sector.index
#  xlim = CELL_META$xlim
#  ylim = CELL_META$ylim
#  circos.rect(xlim[1], 0, xlim[2], 1, col = rand_color(1))
#}, track.height = 0.05, bg.border = NA)

#circos.genomicRect(region[l, , drop = FALSE], ytop = n - i + 1 + 0.4,ybottom = n - i + 1 - 0.4, col = "orange", border = NA)



#############################################################
#circos.track(ylim = c(0, 1), track.height = 0.05, bg.border = NA, panel.fun = function(x, y) {
#  chr = CELL_META$sector.index
#  xlim = CELL_META$xlim
#  ylim = CELL_META$ylim
#  circos.rect(xlim[1], 0, xlim[2], 1, col = rand_color(1))
#})

#color = c("#66C2A5", "#FC8D62", "#8DA0CB", "#E78AC3")



#print(region)
#print(value)


#i = getI(...) 
#anchor_hap <- as.numeric(factor(value[1]))
#print(anchor_hap)
#print(i)

#############################################################
#circos.clear()



#circos.genomicTrackPlotRegion(anchor, stack = TRUE,track.margin=c(0,0), ylim = c(0, 1), bg.border = NA, bg.col = "yellow", track.height = 0.02, panel.fun = function(region, value, ...) {
#  circos.genomicRect(region, value, col = "red", border = NA)
#})



#col = ifelse(value > 1, "deepskyblue3", "brown3")


#color[CELL_META$sector.numeric.index]
#col = ifelse(index > 1, "deepskyblue3", "brown3")


#############################################################


#############################################################
#circos.genomicTrackPlotRegion(anchor, panel.fun = function(region, value, ...) {  #ylim = c(0.5, 5.5),
#  circos.genomicRect(region, value, col = "red", border = NA)
#}, bg.border = NA, bg.col = "yellow", track.height = 0.02)


#cutoff <- "300" 












#chrom_list <- c("chr4","chr5")
#chrom_list <- c("chr10","chr11","chr20")
#chrom_list <- c("chr9","chr4","chr5")
#chrom_list <- c("chr8","chr21")
#chrom_list <- c("chr9","chr22")
#chrom_list <- c("chr4","chr5","chr9","chr11","chr19","chr20","chr21","chr22")
#chrom_list <- c("chr1","chr12","chr17")




#hapA <- hap_tracks %>% select(chr,min_pos,max_pos,value1=meanAcov)
#hapB <- hap_tracks %>% select(chr,min_pos,max_pos,value1=meanBcov) 


#############################################################
#cn_phased_file <- "./cn_phased_may24_tumor_"
#cn_phased_file <- "./cn_phased_may14_"



#phased_data <- read_delim(phased_file,delim="\t",col_names = FALSE) #%>% filter( chr1 %in% chrom_list ) %>% filter( chr2 %in% chrom_list ) #%>% filter( TotalCount > 10 )
#phased_data <- rename(phased_data,"chr1"="X1","pos1"="X2","chr2"="X3","pos2"="X4","hap1"="X5","hap2"="X6","call_type"="X7")
#phased_data <- phased_data %>% filter( chr1 %in% chrom_list ) %>% filter( chr2 %in% chrom_list )
#phased_data <- phased_data %>% mutate(type = hap1*hap2)
#phased_data <- phased_data %>% mutate(type = replace(type,which(type == -1),3))
#phased_data <- phased_data %>% mutate(type = replace(type,which(type == 1 & hap2 == -1),1))
#phased_data <- phased_data %>% mutate(type = replace(type,which(type == 1 & hap2 == 1),2))

#phasedb_links <- phased_data %>% filter(type == 1)
#phaseda_links <- phased_data %>% filter(type == 2)
#mix_links <- phased_data %>% filter(type == 3)
#null_links <- phased_data %>% filter(type == 0)

#b_link1 <- phasedb_links %>% select(chr1,pos1) %>% mutate(end = as.integer(pos1 + band_width)) %>% rename(start=pos1,chr=chr1) %>% as.data.frame()
#b_link2 <- phasedb_links %>% select(chr2,pos2) %>% mutate(end = as.integer(pos2 + band_width)) %>% rename(start=pos2,chr=chr2) %>% as.data.frame()

#a_link1 <- phaseda_links %>% select(chr1,pos1) %>% mutate(end = as.integer(pos1 + band_width)) %>% rename(start=pos1,chr=chr1) %>% as.data.frame()
#a_link2 <- phaseda_links %>% select(chr2,pos2) %>% mutate(end = as.integer(pos2 + band_width)) %>% rename(start=pos2,chr=chr2) %>% as.data.frame()

#m_link1 <- mix_links %>% select(chr1,pos1) %>% mutate(end = as.integer(pos1 + band_width)) %>% rename(start=pos1,chr=chr1) %>% as.data.frame()
#m_link2 <- mix_links %>% select(chr2,pos2) %>% mutate(end = as.integer(pos2 + band_width)) %>% rename(start=pos2,chr=chr2) %>% as.data.frame()

#n_link1 <- null_links %>% select(chr1,pos1) %>% mutate(end = as.integer(pos1 + band_width)) %>% rename(start=pos1,chr=chr1) %>% as.data.frame()
#n_link2 <- null_links %>% select(chr2,pos2) %>% mutate(end = as.integer(pos2 + band_width)) %>% rename(start=pos2,chr=chr2) %>% as.data.frame()






#file_string <- paste(cn_phased_file,chr,".dat",sep="") 
#cn_block_data <- read_delim(file_string ,delim="\t",col_names = FALSE)
#raw_data <- rename(cn_block_data,"index"="X1","bin"="X2","pos"="X3","ref_base"="X4","var_base"="X5","hap"="X6","hapA_cov"="X7","hapB_cov"="X8","block"="X9")
#binned <- raw_data %>% group_by(bin) %>% summarise(meanAcov = mean(hapA_cov), meanBcov = mean(hapB_cov), n=n(),min_pos = min(pos),max_pos = max(pos))
#binned <- binned %>% mutate("chr" = chr);
#hap_tracks <- bind_rows(hap_tracks,binned)

#setwd("~/2018_6_june_workdir/phasing_HCC1954_linker_v1.4/output/")



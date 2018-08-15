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

setwd("./../output/")

#############################################################
chrom_list <- c("chr1","chr2","chr3","chr4","chr5","chr6","chr7","chr8","chr9","chr10","chr11","chr12","chr13","chr14","chr15","chr16","chr17","chr18","chr19","chr20","chr21","chr22","chrX")
#chrom_list <- c("chr2","chr5")

#############################################################
cn_phased_file <- "./cn_phased_may14_"

#############################################################
band_width <- 10000
cutoff_cov <- 150

#############################################################
cutoff <- "300"       
binsize <- "10000"  

hap_tracks = tibble()
for (chr in chrom_list) { 
  file_string <- paste(cn_phased_file,chr,".dat",sep="") 
  cn_block_data <- read_delim(file_string ,delim="\t",col_names = FALSE)
  raw_data <- rename(cn_block_data,"index"="X1","bin"="X2","pos"="X3","ref_base"="X4","var_base"="X5","hap"="X6","hapA_cov"="X7","hapB_cov"="X8","block"="X9")
  binned <- raw_data %>% group_by(bin) %>% summarise(meanAcov = mean(hapA_cov), meanBcov = mean(hapB_cov), n=n(),min_pos = min(pos),max_pos = max(pos))
  binned <- binned %>% mutate("chr" = chr);
  hap_tracks <- bind_rows(hap_tracks,binned)
}
hapA <- hap_tracks %>% select(chr,min_pos,max_pos,value1=meanAcov)
hapB <- hap_tracks %>% select(chr,min_pos,max_pos,value1=meanBcov) 

#############################################################
circos.initializeWithIdeogram(species="hg38",chromosome.index = chrom_list,ideogram.height = 0.01)

#############################################################         #phased_file <- "./phased_sv_sites_apr.dat"
phased_file <- "./phased_sv_sites.dat"
#phased_file <- "./phased_sv_sites_may.dat"
phased_data <- read_delim(phased_file,delim="\t",col_names = FALSE) #%>% filter( chr1 %in% chrom_list ) %>% filter( chr2 %in% chrom_list ) #%>% filter( TotalCount > 10 )
phased_data <- rename(phased_data,"chr1"="X1","pos1"="X2","chr2"="X3","pos2"="X4","hap1"="X5","hap2"="X6","call_type"="X7")
phased_data <- phased_data %>% filter( chr1 %in% chrom_list ) %>% filter( chr2 %in% chrom_list )
phased_data <- phased_data %>% mutate(type = hap1*hap2)
phased_data <- phased_data %>% mutate(type = replace(type,which(type == -1),3))
phased_data <- phased_data %>% mutate(type = replace(type,which(type == 1 & hap2 == -1),1))
phased_data <- phased_data %>% mutate(type = replace(type,which(type == 1 & hap2 == 1),2))

phasedb_links <- phased_data %>% filter(type == 1)
phaseda_links <- phased_data %>% filter(type == 2)
mix_links <- phased_data %>% filter(type == 3)
null_links <- phased_data %>% filter(type == 0)

b_link1 <- phasedb_links %>% select(chr1,pos1) %>% mutate(end = as.integer(pos1 + band_width)) %>% rename(start=pos1,chr=chr1) %>% as.data.frame()
b_link2 <- phasedb_links %>% select(chr2,pos2) %>% mutate(end = as.integer(pos2 + band_width)) %>% rename(start=pos2,chr=chr2) %>% as.data.frame()

a_link1 <- phaseda_links %>% select(chr1,pos1) %>% mutate(end = as.integer(pos1 + band_width)) %>% rename(start=pos1,chr=chr1) %>% as.data.frame()
a_link2 <- phaseda_links %>% select(chr2,pos2) %>% mutate(end = as.integer(pos2 + band_width)) %>% rename(start=pos2,chr=chr2) %>% as.data.frame()

m_link1 <- mix_links %>% select(chr1,pos1) %>% mutate(end = as.integer(pos1 + band_width)) %>% rename(start=pos1,chr=chr1) %>% as.data.frame()
m_link2 <- mix_links %>% select(chr2,pos2) %>% mutate(end = as.integer(pos2 + band_width)) %>% rename(start=pos2,chr=chr2) %>% as.data.frame()

n_link1 <- null_links %>% select(chr1,pos1) %>% mutate(end = as.integer(pos1 + band_width)) %>% rename(start=pos1,chr=chr1) %>% as.data.frame()
n_link2 <- null_links %>% select(chr2,pos2) %>% mutate(end = as.integer(pos2 + band_width)) %>% rename(start=pos2,chr=chr2) %>% as.data.frame()

#############################################################
hapA_high = hapA %>% filter(value1 > cutoff_cov)
hapB_high = hapB %>% filter(value1 > cutoff_cov)
hapA_low = hapA %>% filter(value1 < cutoff_cov)
hapB_low = hapB %>% filter(value1 < cutoff_cov)
combined_hap_high <- list(hapA_high,hapB_high)
combined_hap_low <- list(hapA_low,hapB_low)
combined_hap <- list(hapA,hapB)

#############################################################    #, cell.padding = c(0.0, 0, 0.0, 0)   900  , 800
circos.genomicTrackPlotRegion(combined_hap_high, ylim = c(cutoff_cov, 700), track.margin	= c(0.0, 0.0), cell.padding = c(0.0, 0, 0.0, 0), bg.border = "white", track.height = 0.08, panel.fun = function(region, value, ...) {
  index = getI(...)
  col = ifelse(index > 1, "deepskyblue3", "brown3")
  circos.genomicPoints(region, value, col = col, cex = 0.3, pch = 16) # index    #cex is size   #16   12 , log10(value)
  cell.xlim = get.cell.meta.data("cell.xlim")
  for( h in c( 200, 400 ) ) { circos.lines(cell.xlim, c(h, h), col = "#00000040") }   #, 300, 400, 500
})

#############################################################  #1.0
circos.genomicTrackPlotRegion(combined_hap_low, ylim = c(0, cutoff_cov), track.margin	= c(0.0, 0.0), cell.padding = c(0.0, 0, 0.0, 0), bg.border = "white", track.height = 0.3, panel.fun = function(region, value, ...) {
  index = getI(...)
  col = ifelse(index > 1, "deepskyblue3", "brown3")
  circos.genomicPoints(region, value, col = col, cex = 0.3, pch = 16) # index    #cex is size   #16   12 , log10(value)    #cex = 0.25
  cell.xlim = get.cell.meta.data("cell.xlim")
  for( h in c( 0, 20, 40, 60, 80, 100 ) ) { circos.lines(cell.xlim, c(h, h), col = "#00000040") }   #, 300, 400, 500
})

############################################### ##############
circos.genomicLink(n_link1,n_link2,col="grey50")
circos.genomicLink(m_link1,m_link2,col="darkorchid3",lwd=1.5)
circos.genomicLink(a_link1,a_link2,col="brown3",lwd=1.5) #, col="#3182bd"
circos.genomicLink(b_link1,b_link2,col="deepskyblue3",lwd=1.5)

highlight.chromosome("chr21")
highlight.chromosome("chr5")
highlight.chromosome("chr8")
highlight.chromosome("chr11")
circos.clear()














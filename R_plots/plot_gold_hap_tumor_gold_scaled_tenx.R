# script to read and parse data files
rm(list=ls())
library(ggplot2)
library(readr)
library(tidyr)
library(dplyr)
library(grid)
library(gridExtra)

#############################################################
setwd("~/2018_6_june_workdir/chosen_phase/")
source("load_fun.R")

chromosome <- "chr2" #"chr17"
hap_file <- "./hap_gold_solution_"
binsize <- 10000

#############################################################
#het_file <- "./../phasing_HCC1954_linker_v1.8_illumina_PCRFree_tenx/output/gold_filtered_coverage_may23_tumor_tenx_"
#het_file_normal <- "./../phasing_HCC1954_linker_v1.8_illumina_PCRFree_tenx/output/gold_filtered_coverage_may23_normal_tenx_"
#het_file_bx <- "./../phasing_HCC1954_linker_v1.8_illumina_PCRFree_tenx/output/het_coverage_may22_tumor_tenx_"
#bx_coverage_file <- "./../phasing_HCC1954_linker_v1.8_illumina_PCRFree_tenx/output/bx_unique_binsize_10000_june8_tumor_"
#bx_coverage_normal_file <- "./../phasing_HCC1954_linker_v1.8_illumina_PCRFree_tenx/output/bx_unique_binsize_10000_june8_normal_"
#one_copy <- 0.17 # tenx

het_file <- "./../phasing_HCC1954_linker_v1.8_illumina_PCRFree_tenx/output/gold_filtered_coverage_june18_tumor_EA_tenx_"
het_file_normal <- "./../phasing_HCC1954_linker_v1.8_illumina_PCRFree_tenx/output/gold_filtered_coverage_june18_normal_EA_tenx_"
het_file_bx <- "./../phasing_HCC1954_linker_v1.8_illumina_PCRFree_tenx/output/het_coverage_june18_tumor_EA_tenx_"
bx_coverage_file <- "./../phasing_HCC1954_linker_v1.8_illumina_PCRFree_tenx/output/bx_unique_binsize_10000_june18_tumor_EA_"
bx_coverage_normal_file <- "./../phasing_HCC1954_linker_v1.8_illumina_PCRFree_tenx/output/bx_unique_binsize_10000_june18_normal_EA_"
one_copy <- 0.24 # tenx_EA


##############################################################
hdata <- load_hap(hap_file,binsize,chromosome)
blocks <- hdata %>% group_by(block) %>% summarise(min_pos = min(pos),max_pos = max(pos))
binned <- load_het_cov(het_file,het_file_normal,one_copy,chromosome,hdata)
binned <- binned %>% filter(n>5)
#binned <- binned %>% filter(n>10)

##############################################################
un_bin <- load_het_cov_tenx(het_file_bx,chromosome,hdata)  #,het_file_normal_bx
un_bin <- un_bin %>% filter(n>5)

binned_bx_unique <- inner_join(binned,un_bin,by = c("bin","bin"))
#binned_bx_unique <- left_join(binned,un_bin,by = c("bin","bin"))
#binned_bx_unique <- binned_bx_unique %>% mutate(tumorAnorm_unique = unique_meanAcov/one_copy,tumorBnorm_unique = unique_meanBcov/one_copy)

binbx_data <- load_bx_coverage(bx_coverage_file,binsize,chromosome)
binbx_data_normal <- load_bx_coverage(bx_coverage_normal_file,binsize,chromosome)

##############################################################
binbx_data <- left_join(binbx_data,binbx_data_normal,by=c("bin"="bin"))
binbx_data <- left_join(binned_bx_unique,binbx_data,by=c("bin"="bin"))
binbx_data <- binbx_data %>% dplyr::rename("ubx"="ubx.x","normal_ubx"="ubx.y","totbx"="totbx.x","normal_totbx"="totbx.y") 
binbx_data <- binbx_data %>% mutate(ubx_newA = fracA*ubx,ubx_newB = (1-fracA)*ubx)
binbx_data <- binbx_data %>% mutate(frac_ubx = ubx/normal_ubx) %>% mutate(ubx_meanAcov = fracA*frac_ubx,ubx_meanBcov = (1-fracA)*frac_ubx)
binbx_data <- binbx_data %>% mutate(tumorAnorm_unique = ubx_meanAcov/one_copy,tumorBnorm_unique = ubx_meanBcov/one_copy)


#############################################################
phased_file <- "./phased_sv_sites_gold.dat"
sv_out <- load_phasedsv_file(phased_file,chromosome)

######################################################################
hg38_data <- load_hg38_data()

#ylimits <- c(0,500) 
#ylimits <- c(0,40)
ylimits <- c(0,9)
#ylimits <- c(0,400) 
#xlimits <- c(0,2.5E8)
xlimits <- c(0,max(hg38_data$chromEnd))


log_scale_break <- 2
addition_num <- 1
#log_scale_break <- 4
#addition_num <- 2
vline_yend <- 6


binbx_data <- binbx_data %>% mutate( ln2_meanA = case_when( tumorAnorm_unique > log_scale_break ~ log2(tumorAnorm_unique) + addition_num, tumorAnorm_unique <= log_scale_break ~ tumorAnorm_unique ) , ln2_meanB = case_when( tumorBnorm_unique > log_scale_break ~ log2(tumorBnorm_unique) + addition_num, tumorBnorm_unique <= log_scale_break ~ tumorBnorm_unique ) )

#p1 <- ggplot(data=jdata) + geom_point(aes(x = pos,y = hapA_cov),color="brown3",size=0.5) + geom_point(aes(x = pos,y = hapB_cov),color="deepskyblue3",size=0.5) + scale_y_continuous(limits = ylimits) + scale_x_continuous(limits = xlimits)
p1 <- ggplot(data=binbx_data) + ggtitle(paste("10X EA ",chromosome)) +
      geom_vline(data = sv_out$null_inter,aes(xintercept = single_chr_pos),color="grey50",alpha=0.6,size=0.8) +
      geom_vline(data = sv_out$hapB_inter,aes(xintercept = single_chr_pos),color="deepskyblue3",alpha=0.6,size=0.8) +
      geom_vline(data = sv_out$hapA_inter,aes(xintercept = single_chr_pos),color="brown3",alpha=0.6,size=0.8) +
      geom_segment(data = sv_out$null_links,aes(x = pos1, xend = pos1, y=0, yend=vline_yend),color = "grey50",size=0.8,alpha=0.6,linetype="dashed") +  #"dashed"   "longdash"
      geom_segment(data = sv_out$null_links,aes(x = pos2, xend = pos2, y=0, yend=vline_yend),color = "grey50",size=0.8,alpha=0.6,linetype="dashed") +
      geom_segment(data = sv_out$phasedb_links,aes(x = pos1, xend = pos1, y=0, yend=vline_yend),color = "deepskyblue3",size=0.8,alpha=0.6,linetype="dashed") + 
      geom_segment(data = sv_out$phasedb_links,aes(x = pos2, xend = pos2, y=0, yend=vline_yend),color = "deepskyblue3",size=0.8,alpha=0.6,linetype="dashed") +
      geom_segment(data = sv_out$phaseda_links,aes(x = pos1, xend = pos1, y=0, yend=vline_yend),color = "brown3",size=0.8,alpha=0.6,linetype="dashed") + 
      geom_segment(data = sv_out$phaseda_links,aes(x = pos2, xend = pos2, y=0, yend=vline_yend),color = "brown3",size=0.8,alpha=0.6,linetype="dashed") +
      geom_segment(data = sv_out$mix_links,aes(x = pos1, xend = pos1, y=0, yend=vline_yend),color = "darkorchid3",size=0.8,alpha=0.6,linetype="dashed") + 
      geom_segment(data = sv_out$mix_links,aes(x = pos2, xend = pos2, y=0, yend=vline_yend),color = "darkorchid3",size=0.8,alpha=0.6,linetype="dashed") +
      #geom_point(aes(x = middle_pos,y = tumorAnorm),color="brown3",size=0.7) + 
      #geom_point(aes(x = middle_pos,y = tumorBnorm),color="deepskyblue3",size=0.7) +  
      geom_point(aes(x = middle_pos,y = ln2_meanA),color="brown3",size=0.7) + 
      geom_point(aes(x = middle_pos,y = ln2_meanB),color="deepskyblue3",size=0.7) +  
      scale_x_continuous(limits = xlimits,position = "top") + 
      theme_minimal() + xlab("position") + ylab("copy number") +
      scale_y_continuous(limits = ylimits,breaks=c(0,1,2,3,4,5,6,7,8,9),label=c("0","1","2","4","8","16","32","64","128","256"),expand = c(0,0)) +
      #scale_y_continuous(limits = ylimits,breaks=c(0,1,2,3,4,5,6,8,10,12),label=c("0","1","2","3","4","5","6","8","10","12"),expand = c(0,0)) +
      theme(axis.title.x=element_blank(),axis.text.x=element_blank())


p2 <- ggplot(data=blocks) + geom_segment(aes(x=min_pos,xend=max_pos,y=-0.0,yend=0.0),color= "grey25",size=8) + scale_x_continuous(limits = xlimits) + theme_minimal() + theme(axis.line=element_blank(),axis.text.x=element_blank(),axis.text.y=element_blank(),axis.ticks=element_blank(),axis.title.x=element_blank(),axis.title.y=element_blank())   #+ theme_void()      #+ theme_minimal()


cols <- c("gneg" = "#92B3F9", "gpos25" = "#0F214D", "gpos50" = "#0F214D", "acen" = "#EAC142","gpos100" = "#0F214D","gpos75" = "#0F214D","gvar"="#FE6D39","stalk"="red")
p3 <- ggplot(data=hg38_data) + geom_segment(mapping=aes( x=chromStart, xend=chromEnd, y=1, yend=1, color=gieStain ),size=8,stat='identity') + scale_x_continuous(limits = xlimits) + theme_minimal() + theme(legend.position="none",axis.line=element_blank(),axis.text.y=element_blank(),axis.ticks=element_blank(),axis.title.x=element_blank()) + scale_color_manual(values = cols) + ylab(chromosome) 
#p3 <- ggplot(data=hg38_data) + geom_segment(mapping=aes( x=chromStart, xend=chromEnd, y=1, yend=1, color=gieStain ),size=8,stat='identity') + scale_x_continuous(limits = xlimits) + theme_minimal() + theme(legend.position="none",axis.line=element_blank(),axis.text.y=element_blank(),axis.ticks=element_blank(),axis.title.x=element_blank(),axis.title.y=element_blank()) + scale_color_manual(values = cols) #+ scale_color_manual(values = cols)  ,axis.text.x=element_blank()


p1 <- ggplot_gtable(ggplot_build(p1))   
p2 <- ggplot_gtable(ggplot_build(p2))
p3 <- ggplot_gtable(ggplot_build(p3))
maxWidth = unit.pmax(p1$widths[2:3], p2$widths[2:3], p3$widths[2:3] )  
p1$widths[2:3] <- maxWidth              
p2$widths[2:3] <- maxWidth
p3$widths[2:3] <- maxWidth

grid.arrange(p1, p2, p3,ncol = 1, heights = c(6, 0.4,0.5))  



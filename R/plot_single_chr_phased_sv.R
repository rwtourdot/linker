# script to read and parse data files
rm(list=ls())
library(ggplot2)
library(readr)
library(tidyr)
library(dplyr)
library(grid)
library(gridExtra)

#############################################################
setwd("./../output/")
cn_phased_file <- "./cn_phased_may14_"
chromosome <- "chr12"

#############################################################
file_string <- paste(cn_phased_file,chromosome,".dat",sep="")    #file_string <- paste("./cn_phased_apr24_",chromosome,".dat",sep="")
cn_block_data <- read_delim(file_string ,delim="\t",col_names = FALSE)
raw_data <- dplyr::rename(cn_block_data,"index"="X1","bin"="X2","pos"="X3","ref_base"="X4","var_base"="X5","hap"="X6","hapA_cov"="X7","hapB_cov"="X8","block"="X9","flip"="X10")
binned <- raw_data %>% group_by(bin) %>% dplyr::summarise(meanAcov = mean(hapA_cov), meanBcov = mean(hapB_cov), n=n(),min_pos = min(pos),max_pos = max(pos))
binned <- binned %>% mutate("chr" = chromosome)
blocks <- raw_data %>% group_by(block) %>% summarise(min_pos = min(pos),max_pos = max(pos))

######################################################################
file_string <- "./../hg38_info/chromosome_banding_hg38.txt"
hg38_data <- read_delim(file_string,delim="\t",col_names = TRUE) %>% dplyr::rename("chrom"=`#chrom`) %>% filter(chrom == chromosome)

#############################################################
phased_file <- "./phased_sv_sites.dat"                          #phased_file <- "./phased_sv_sites_apr.dat"
phased_data <- read_delim(phased_file,delim="\t",col_names = FALSE) #%>% filter( chr1 %in% chrom_list ) %>% filter( chr2 %in% chrom_list ) #%>% filter( TotalCount > 10 )
phased_data <- dplyr::rename(phased_data,"chr1"="X1","pos1"="X2","chr2"="X3","pos2"="X4","hap1"="X5","hap2"="X6","sv_type"="X7")

#############################################################
intra = phased_data %>% filter(chr1 == chromosome & chr2 == chromosome)
inter = phased_data %>% filter(chr1 == chromosome | chr2 == chromosome) %>% filter( chr2 != chr1 )

inter <- inter %>% mutate(single_chr_pos = case_when(chr1==chromosome ~ pos1, chr2==chromosome ~ pos2))
inter <- inter %>% mutate(sv_to = case_when(chr1==chromosome ~ chr2, chr2==chromosome ~ chr1))
inter <- inter %>% mutate(single_hap = case_when(chr1==chromosome ~ hap1, chr2==chromosome ~ hap2))

#############################################################
null_inter <- inter %>% filter(single_hap == 0)
hapA_inter <- inter %>% filter(single_hap == 1)
hapB_inter <- inter %>% filter(single_hap == -1)
null_inter <- null_inter %>% add_row(chr1 = chromosome, pos1 = 0, chr2 = chromosome, pos2 = 1, hap1=1, hap2=-1, sv_type = "inter", single_chr_pos = -10000000)
hapB_inter <- hapB_inter %>% add_row(chr1 = chromosome, pos1 = 0, chr2 = chromosome, pos2 = 1, hap1=1, hap2=-1, sv_type = "inter", single_chr_pos = -10000000)
hapA_inter <- hapA_inter %>% add_row(chr1 = chromosome, pos1 = 0, chr2 = chromosome, pos2 = 1, hap1=1, hap2=-1, sv_type = "inter", single_chr_pos = -10000000)

#############################################################
intra <- intra %>% mutate(type = hap1*hap2)
intra <- intra %>% mutate(type = replace(type,which(type == -1),3))
intra <- intra %>% mutate(type = replace(type,which(type == 1 & hap2 == -1),1))
intra <- intra %>% mutate(type = replace(type,which(type == 1 & hap2 == 1),2))

phasedb_links <- intra %>% filter( type == 1 )
phaseda_links <- intra %>% filter( type == 2 )
mix_links <- intra %>% filter( type == 3 )
null_links <- intra %>% filter( type == 0 )

phaseda_links <- phaseda_links %>% add_row(chr1 = chromosome, pos1 = -10000000, chr2 = chromosome, pos2 = -10000001, hap1=1, hap2=1, sv_type = "intra", type = 2)
phasedb_links <- phasedb_links %>% add_row(chr1 = chromosome, pos1 = -10000000, chr2 = chromosome, pos2 = -10000001, hap1=-1, hap2=-1, sv_type = "intra", type = 1)
mix_links <- mix_links %>% add_row(chr1 = chromosome, pos1 = -10000000, chr2 = chromosome, pos2 = -10000001, hap1=1, hap2=-1, sv_type = "intra", type = 3)
null_links <- null_links %>% add_row(chr1 = chromosome, pos1 = -10000000, chr2 = chromosome, pos2 = -10000001, hap1=0, hap2=-1, sv_type = "intra", type = 0)

#############################################################
xlimits <- c(0,2.5E8)
#ximits <- c(0,0.5E8)
#xlimits <- c(0.11E8,0.13E8)
ylimits <- c(0,400) 
vline_yend_inter <- 300
vline_yend <- 400

p1 <- ggplot(data=binned) + geom_segment(data = null_inter,aes(x = single_chr_pos, xend = single_chr_pos, y=0, yend=vline_yend_inter),color="grey50",alpha=0.8,size=0.6) +
                      geom_segment(data = hapB_inter,aes(x = single_chr_pos, xend = single_chr_pos, y=0, yend=vline_yend_inter),color="deepskyblue3",alpha=0.8,size=0.6) +
                      geom_segment(data = hapA_inter,aes(x = single_chr_pos, xend = single_chr_pos, y=0, yend=vline_yend_inter),color="brown3",alpha=0.8,size=0.6) +
                      geom_point(aes(x = (min_pos + max_pos)/2.0,y = meanAcov),color="brown3",size=0.8) +
                      geom_point(aes(x = (min_pos + max_pos)/2.0,y = meanBcov),color="deepskyblue3",size=0.8) +
                      geom_segment(data = null_links,aes(x = pos1, xend = pos1, y=0, yend=vline_yend),color = "grey50",size=0.6,alpha=0.8,linetype="dashed") +  #"dashed"   "longdash"
                      geom_segment(data = null_links,aes(x = pos2, xend = pos2, y=0, yend=vline_yend),color = "grey50",size=0.6,alpha=0.8,linetype="dashed") +
                      geom_segment(data = phasedb_links,aes(x = pos1, xend = pos1, y=0, yend=vline_yend),color = "deepskyblue3",size=0.6,alpha=0.8,linetype="dashed") + 
                      geom_segment(data = phasedb_links,aes(x = pos2, xend = pos2, y=0, yend=vline_yend),color = "deepskyblue3",size=0.6,alpha=0.8,linetype="dashed") +
                      geom_segment(data = phaseda_links,aes(x = pos1, xend = pos1, y=0, yend=vline_yend),color = "brown3",size=0.6,alpha=0.8,linetype="dashed") + 
                      geom_segment(data = phaseda_links,aes(x = pos2, xend = pos2, y=0, yend=vline_yend),color = "brown3",size=0.6,alpha=0.8,linetype="dashed") +
                      geom_segment(data = mix_links,aes(x = pos1, xend = pos1, y=0, yend=vline_yend),color = "darkorchid3",size=0.6,alpha=0.8,linetype="dashed") + 
                      geom_segment(data = mix_links,aes(x = pos2, xend = pos2, y=0, yend=vline_yend),color = "darkorchid3",size=0.6,alpha=0.8,linetype="dashed") +
                      theme_minimal() + xlab("position") + ylab("copy number") +
                      scale_y_continuous(limits = ylimits) +
                      theme(axis.title.x=element_blank(),axis.text.x=element_blank()) + 
                      scale_x_continuous(limits = xlimits,position = "top")

p3 <- ggplot(data=blocks) + geom_segment(aes(x=min_pos,xend=max_pos,y=-0.0,yend=0.0),color= "grey25",size=8) + scale_x_continuous(limits = xlimits) + theme_minimal() + theme(axis.line=element_blank(),axis.text.x=element_blank(),axis.text.y=element_blank(),axis.ticks=element_blank(),axis.title.x=element_blank(),axis.title.y=element_blank())   #+ theme_void()      #+ theme_minimal()

cols <- c("gneg" = "#92B3F9", "gpos25" = "#0F214D", "gpos50" = "#0F214D", "acen" = "#EAC142","gpos100" = "#0F214D","gpos75" = "#0F214D","gvar"="#FE6D39","stalk"="red")
p4 <- ggplot(data=hg38_data) + geom_segment(mapping=aes( x=chromStart, xend=chromEnd, y=1, yend=1, color=gieStain ),size=8,stat='identity') + scale_x_continuous(limits = xlimits) + theme_minimal() + theme(legend.position="none",axis.line=element_blank(),axis.text.y=element_blank(),axis.ticks=element_blank(),axis.title.x=element_blank(),axis.title.y=element_blank()) + scale_color_manual(values = cols) #+ scale_color_manual(values = cols)  ,axis.text.x=element_blank()

p1 <- ggplot_gtable(ggplot_build(p1))    #p2 <- ggplot_gtable(ggplot_build(p2))
p3 <- ggplot_gtable(ggplot_build(p3))
p4 <- ggplot_gtable(ggplot_build(p4))
maxWidth = unit.pmax(p1$widths[2:3], p3$widths[2:3], p4$widths[2:3])   #, p2$widths[2:3]
p1$widths[2:3] <- maxWidth               #p2$widths[2:3] <- maxWidth
p3$widths[2:3] <- maxWidth
p4$widths[2:3] <- maxWidth

grid.arrange(p1, p3, p4, ncol = 1, heights = c(6, 0.4, 0.4))               
                      



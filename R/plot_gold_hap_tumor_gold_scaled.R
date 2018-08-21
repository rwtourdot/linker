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

chromosome <- "chrX"   #"chr17"   "chr11"
hap_file <- "./hap_gold_solution_"
binsize <- 10000

#############################################################
######het_file <- "./../phasing_HCC1954_linker_v1.8_illumina_PCRFree_tenx/output/gold_filtered_coverage_may23_tumor_tenx_"
######het_file_normal <- "./../phasing_HCC1954_linker_v1.8_illumina_PCRFree_tenx/output/gold_filtered_coverage_may23_normal_tenx_"
######one_copy <- 0.21  # tenx
######het_file <- "./../phasing_HCC1954_linker_v1.8_illumina_PCRFree_tenx/output/gold_filtered_coverage_june18_tumor_EA_tenx_"
######het_file_normal <- "./../phasing_HCC1954_linker_v1.8_illumina_PCRFree_tenx/output/gold_filtered_coverage_june18_normal_EA_tenx_"
######one_copy <- 0.21 # tenx_EA
#het_file <- "./../phasing_HCC1954_linker_v1.8_illumina_PCRFree_tenx/output/gold_filtered_coverage_june11_PCRFree_tumor_"
#het_file_normal <- "./../phasing_HCC1954_linker_v1.8_illumina_PCRFree_tenx/output/gold_filtered_coverage_june11_PCRFree_normal_"
#one_copy <- 0.165 # PCRFree
het_file <- "./../phasing_HCC1954_linker_v1.8_illumina_PCRFree_tenx/output/gold_filtered_coverage_june11_illumina_tumor_"
het_file_normal <- "./../phasing_HCC1954_linker_v1.8_illumina_PCRFree_tenx/output/gold_filtered_coverage_june11_illumina_normal_"
one_copy <- 0.19 # illumina
#het_file <- "./../phasing_HCC1954_linker_v1.8_illumina_PCRFree_tenx/output/gold_filtered_coverage_june21_pacbio_tumor_"
#het_file_normal <- "./../phasing_HCC1954_linker_v1.8_illumina_PCRFree_tenx/output/gold_filtered_coverage_june21_pacbio_normal_"
#one_copy <- 0.25 # pacbio


##############################################################
hdata <- load_hap(hap_file,binsize,chromosome)
blocks <- hdata %>% group_by(block) %>% summarise(min_pos = min(pos),max_pos = max(pos))
binned <- load_het_cov(het_file,het_file_normal,one_copy,chromosome,hdata)
binned <- binned %>% filter(n>5)

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
#xlimits <- c(5E7,11E7)

log_scale_break <- 2
addition_num <- 1

vline_yend <- 6

binned <- binned %>% mutate( ln2_meanA = case_when( tumorAnorm > log_scale_break ~ log2(tumorAnorm) + addition_num, tumorAnorm <= log_scale_break ~ tumorAnorm ) , ln2_meanB = case_when( tumorBnorm > log_scale_break ~ log2(tumorBnorm) + addition_num, tumorBnorm <= log_scale_break ~ tumorBnorm ) )

#p1 <- ggplot(data=jdata) + geom_point(aes(x = pos,y = hapA_cov),color="brown3",size=0.5) + geom_point(aes(x = pos,y = hapB_cov),color="deepskyblue3",size=0.5) + scale_y_continuous(limits = ylimits) + scale_x_continuous(limits = xlimits)
p1 <- ggplot(data=binned) + ggtitle(paste("Pacbio ",chromosome)) +
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
p3 <- ggplot(data=hg38_data) + geom_segment(mapping=aes( x=chromStart, xend=chromEnd, y=1, yend=1, color=gieStain ),size=8,stat='identity') + scale_x_continuous(limits = xlimits) + theme_minimal() + theme(legend.position="none",axis.line=element_blank(),axis.text.y=element_blank(),axis.ticks=element_blank(),axis.title.x=element_blank()) + scale_color_manual(values = cols) + ylab(chromosome) #+ scale_color_manual(values = cols)  ,axis.text.x=element_blank()


p1 <- ggplot_gtable(ggplot_build(p1))   
p2 <- ggplot_gtable(ggplot_build(p2))
p3 <- ggplot_gtable(ggplot_build(p3))
maxWidth = unit.pmax(p1$widths[2:3], p2$widths[2:3], p3$widths[2:3] )  
p1$widths[2:3] <- maxWidth              
p2$widths[2:3] <- maxWidth
p3$widths[2:3] <- maxWidth

grid.arrange(p1, p2, p3,ncol = 1, heights = c(6, 0.4,0.5))  








#phased_data <- read_delim(phased_file,delim="\t",col_names = FALSE) #%>% filter( chr1 %in% chrom_list ) %>% filter( chr2 %in% chrom_list ) #%>% filter( TotalCount > 10 )
#phased_data <- dplyr::rename(phased_data,"chr1"="X1","pos1"="X2","chr2"="X3","pos2"="X4","hap1"="X5","hap2"="X6","sv_type"="X7")

#############################################################
#intra <- phased_data %>% filter(chr1 == chromosome & chr2 == chromosome)
#inter <- phased_data %>% filter(chr1 == chromosome | chr2 == chromosome) %>% filter( chr2 != chr1 )

#inter <- inter %>% mutate(single_chr_pos = case_when(chr1==chromosome ~ pos1, chr2==chromosome ~ pos2))
#inter <- inter %>% mutate(sv_to = case_when(chr1==chromosome ~ chr2, chr2==chromosome ~ chr1))
#inter <- inter %>% mutate(single_hap = case_when(chr1==chromosome ~ hap1, chr2==chromosome ~ hap2))

#null_inter <- inter %>% filter(single_hap == 0)
#hapA_inter <- inter %>% filter(single_hap == 1)
#hapB_inter <- inter %>% filter(single_hap == -1)

#############################################################
#intra <- intra %>% mutate(type = hap1*hap2)
#intra <- intra %>% mutate(type = replace(type,which(type == -1),3))
#intra <- intra %>% mutate(type = replace(type,which(type == 1 & hap2 == -1),1))
#intra <- intra %>% mutate(type = replace(type,which(type == 1 & hap2 == 1),2))

#phasedb_links <- intra %>% filter( type == 1 )
#phaseda_links <- intra %>% filter( type == 2 )
#mix_links <- intra %>% filter( type == 3 )
#null_links <- intra %>% filter( type == 0 )

#phaseda_links <- phaseda_links %>% add_row(chr1 = chromosome, pos1 = -1000, chr2 = chromosome, pos2 = -1001, hap1=1, hap2=1, sv_type = "intra", type = 2)
#phasedb_links <- phasedb_links %>% add_row(chr1 = chromosome, pos1 = -1000, chr2 = chromosome, pos2 = -1001, hap1=-1, hap2=-1, sv_type = "intra", type = 1)
#mix_links <- mix_links %>% add_row(chr1 = chromosome, pos1 = -1000, chr2 = chromosome, pos2 = -1001, hap1=1, hap2=-1, sv_type = "intra", type = 3)
#null_links <- null_links %>% add_row(chr1 = chromosome, pos1 = -1000, chr2 = chromosome, pos2 = -1001, hap1=0, hap2=-1, sv_type = "intra", type = 0)


#file_string <- "./../hg38_info/chromosome_banding_hg38.txt"
#hg38_data <- read_delim(file_string,delim="\t",col_names = TRUE) %>% dplyr::rename("chrom"=`#chrom`) %>% filter(chrom == chromosome)


#############################################################
#hap_data <- paste(hap_file,chromosome,".dat",sep="") %>% read_delim(delim="\t",col_names = FALSE)
#raw_data <- dplyr::rename(hap_data,"index"="X1","pos"="X2","ref"="X3","var"="X4","hap"="X5","hapA_cov"="X6","hapB_cov"="X7","block"="X10")
#hdata <- raw_data %>% separate(var,c("a","b","rbase","vbase"),"_",extra = "drop") %>% select("pos","rbase","vbase","hap","block")
#hdata <- hdata %>% mutate(bin = floor(pos/binsize))
#blocks <- raw_data %>% group_by(block) %>% summarise(min_pos = min(pos),max_pos = max(pos))



#############################################################
#het_data <- paste(het_file,chromosome,".dat",sep="") %>% read_delim(delim="\t",col_names = FALSE)
#raw_data <- dplyr::rename(het_data,"pos"="X2","ref:var"="X3","ibase"="X4","dbase"="X5","gbase"="X6","cbase"="X7","abase"="X8","tbase"="X9","tot_cov"="X10")
#cdata <- raw_data %>% separate(col = dbase, into = c("d","dbase"), sep = "\\|") %>% separate(col = ibase, into = c("i","ibase"), sep = "\\|") %>% separate(col = gbase, into = c("g","gbase"), sep = "\\|") %>% separate(col = cbase, into = c("c","cbase"), sep = "\\|") %>% separate(col = abase, into = c("a","abase"), sep = "\\|") %>% separate(col = tbase, into = c("t","tbase"), sep = "\\|")
#cdata <- cdata %>% select("pos","ibase","dbase","gbase","cbase","abase","tbase","tot_cov")

#############################################################
#het_data_normal <- paste(het_file_normal,chromosome,".dat",sep="") %>% read_delim(delim="\t",col_names = FALSE)
#raw_data_normal <- dplyr::rename(het_data_normal,"pos"="X2","ref:var"="X3","ibase"="X4","dbase"="X5","gbase"="X6","cbase"="X7","abase"="X8","tbase"="X9","tot_cov"="X10")
#cdata_normal <- raw_data_normal %>% separate(col = dbase, into = c("d","dbase"), sep = "\\|") %>% separate(col = ibase, into = c("i","ibase"), sep = "\\|") %>% separate(col = gbase, into = c("g","gbase"), sep = "\\|") %>% separate(col = cbase, into = c("c","cbase"), sep = "\\|") %>% separate(col = abase, into = c("a","abase"), sep = "\\|") %>% separate(col = tbase, into = c("t","tbase"), sep = "\\|")
#cdata_normal <- cdata_normal %>% select("pos","ibase","dbase","gbase","cbase","abase","tbase","tot_cov")

#############################################################
#jdata_temp <- left_join(hdata,cdata,by = c("pos","pos"))  #,cdata_normal ,"pos"
#jdata <- left_join(jdata_temp,cdata_normal,by = c("pos","pos"))
#jdata <- jdata %>% mutate(hapA_base = case_when(hap==1 ~ rbase, hap==-1 ~ vbase)) %>% mutate(hapB_base = case_when(hap==1 ~ vbase, hap==-1 ~ rbase))
#############################################################
#jdata <- jdata %>% mutate(hapA_cov = case_when(hapA_base=="A" ~ as.numeric(abase.x), hapA_base=="C" ~ as.numeric(cbase.x), hapA_base=="G" ~ as.numeric(gbase.x),  hapA_base=="T" ~ as.numeric(tbase.x) ))
#jdata <- jdata %>% mutate(hapB_cov = case_when(hapB_base=="A" ~ as.numeric(abase.x), hapB_base=="C" ~ as.numeric(cbase.x), hapB_base=="G" ~ as.numeric(gbase.x),  hapB_base=="T" ~ as.numeric(tbase.x) ))
#jdata <- jdata %>% mutate(hapA_cov_normal = case_when(hapA_base=="A" ~ as.numeric(abase.y), hapA_base=="C" ~ as.numeric(cbase.y), hapA_base=="G" ~ as.numeric(gbase.y),  hapA_base=="T" ~ as.numeric(tbase.y) ))
#jdata <- jdata %>% mutate(hapB_cov_normal = case_when(hapB_base=="A" ~ as.numeric(abase.y), hapB_base=="C" ~ as.numeric(cbase.y), hapB_base=="G" ~ as.numeric(gbase.y),  hapB_base=="T" ~ as.numeric(tbase.y) ))
#jdata <- drop_na(jdata)

#binned <- jdata %>% group_by(bin) %>% dplyr::summarise(meanAcov = mean(hapA_cov), meanBcov = mean(hapB_cov), meantotcov = mean(tot_cov.x),meanAcov_normal = mean(hapA_cov_normal), meanBcov_normal = mean(hapB_cov_normal), meantotcov_norm = mean(tot_cov.y), n=n(),min_pos = min(pos),max_pos = max(pos))
#binned <- binned %>% mutate( tumorA = meanAcov/meantotcov_norm, tumorB = meanBcov/meantotcov_norm ) %>% mutate( middle_pos = ((min_pos + max_pos)/2.0) )
#binned <- binned %>% mutate(tumorAnorm = tumorA/one_copy, tumorBnorm = tumorB/one_copy)
#binned <- binned %>% filter(n>10)




#vline_yend_inter <- 40

# scale_y_continuous(limits = ylimits) +

#jdata <- jdata 
#file_string <- paste(hap_file,chromosome,".dat",sep="")    #file_string <- paste("./cn_phased_apr24_",chromosome,".dat",sep="") file_string ,


####### scale by copy number
#one_copy <- 0.215  # tenx
#one_copy <- 0.16 # PCRFree
#one_copy <- 0.19 # illumina
#one_copy <- 0.205 # tenx_EA

#file_string <- paste(het_file,chromosome,".dat",sep="")    # file_string ,
#file_string_normal <- paste(het_file_normal,chromosome,".dat",sep="") %>% read_delim(delim="\t",col_names = FALSE)  #file_string_normal ,

#one_copy <- 0.215  # tenx
#one_copy <- 0.16 # PCRFree
#one_copy <- 0.19 # illumina
#one_copy <- 0.205 # tenx_EA

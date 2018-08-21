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

chromosome <- "chr4"
hap_file <- "./hap_gold_solution_"
binsize <- 10000

#############################################################
file_string <- paste(hap_file,chromosome,".dat",sep="")    #file_string <- paste("./cn_phased_apr24_",chromosome,".dat",sep="")
hap_data <- read_delim(file_string ,delim="\t",col_names = FALSE)
raw_data <- dplyr::rename(hap_data,"index"="X1","pos"="X2","ref"="X3","var"="X4","hap"="X5","hapA_cov"="X6","hapB_cov"="X7","block"="X10")
hdata <- raw_data %>% separate(var,c("a","b","rbase","vbase"),"_",extra = "drop") %>% select("pos","rbase","vbase","hap","block")
hdata <- hdata %>% mutate(bin = floor(pos/binsize))

blocks <- raw_data %>% group_by(block) %>% summarise(min_pos = min(pos),max_pos = max(pos))

#############################################################

#het_file <- "./../phasing_HCC1954_linker_v1.8_illumina_PCRFree_tenx/output/gold_filtered_coverage_may23_tumor_tenx_"
#het_file <- "./../phasing_HCC1954_linker_v1.8_illumina_PCRFree_tenx/output/gold_filtered_coverage_may23_normal_tenx_"
#het_file <- "./../phasing_HCC1954_linker_v1.8_illumina_PCRFree_tenx/output/gold_filtered_coverage_june11_PCRFree_tumor_"
#het_file <- "./../phasing_HCC1954_linker_v1.8_illumina_PCRFree_tenx/output/gold_filtered_coverage_june11_PCRFree_normal_"
het_file <- "./../phasing_HCC1954_linker_v1.8_illumina_PCRFree_tenx/output/gold_filtered_coverage_june11_illumina_tumor_"
#het_file <- "./../phasing_HCC1954_linker_v1.8_illumina_PCRFree_tenx/output/gold_filtered_coverage_june11_illumina_normal_"

file_string <- paste(het_file,chromosome,".dat",sep="")  
het_data <- read_delim(file_string ,delim="\t",col_names = FALSE)
raw_data <- dplyr::rename(het_data,"pos"="X2","ref:var"="X3","ibase"="X4","dbase"="X5","gbase"="X6","cbase"="X7","abase"="X8","tbase"="X9","tot_cov"="X10")
cdata <- raw_data %>% separate(col = dbase, into = c("d","dbase"), sep = "\\|") %>% separate(col = ibase, into = c("i","ibase"), sep = "\\|") %>% separate(col = gbase, into = c("g","gbase"), sep = "\\|") %>% separate(col = cbase, into = c("c","cbase"), sep = "\\|") %>% separate(col = abase, into = c("a","abase"), sep = "\\|") %>% separate(col = tbase, into = c("t","tbase"), sep = "\\|")
cdata <- cdata %>% select("pos","ibase","dbase","gbase","cbase","abase","tbase","tot_cov")

jdata <- left_join(hdata,cdata,by = c("pos","pos"))
jdata <- jdata %>% mutate(hapA_base = case_when(hap==1 ~ rbase, hap==-1 ~ vbase))
jdata <- jdata %>% mutate(hapB_base = case_when(hap==1 ~ vbase, hap==-1 ~ rbase))
jdata <- jdata %>% mutate(hapA_cov = case_when(hapA_base=="A" ~ as.numeric(abase), hapA_base=="C" ~ as.numeric(cbase), hapA_base=="G" ~ as.numeric(gbase),  hapA_base=="T" ~ as.numeric(tbase) ))
jdata <- jdata %>% mutate(hapB_cov = case_when(hapB_base=="A" ~ as.numeric(abase), hapB_base=="C" ~ as.numeric(cbase), hapB_base=="G" ~ as.numeric(gbase),  hapB_base=="T" ~ as.numeric(tbase) ))
jdata <- drop_na(jdata)

binned <- jdata %>% group_by(bin) %>% dplyr::summarise(meanAcov = mean(hapA_cov), meanBcov = mean(hapB_cov), meantotcov = mean(tot_cov), n=n(),min_pos = min(pos),max_pos = max(pos))
binned <- binned %>% mutate(middle_pos = ((min_pos + max_pos)/2.0))
binned <- binned %>% filter(n>10)
#binned <- binned %>% filter(n>5)

ylimits <- c(0,300) 
#ylimits <- c(0,40)
#ylimits <- c(0,400) 
#ylimits <- c(0,180) 
#ylimits <- c(0,4070)
#xlimits <- c(1.08E8,1.12E8)
#xlimits <- c(2.6E7,3.8E7)
#xlimits <- c(7.4E7,7.55E7)
#xlimits <- c(0,5E7)
#xlimits <- c(2E7,2.1E7)
#xlimits <- c(0.85E8,0.9E8)
xlimits <- c(0,2.5E8)
#xlimits <- c(0.2E8,0.3E8)
#xlimits <- c(0.24E8,0.26E8)
#xlimits <- c(0.0525E8,0.0575E8)
#xlimits <- c(0,0.8E8)
#xlimits <- c(3.1E7,3.5E7)
#xlimits <- c(8E7,9E7)

#p1 <- ggplot(data=jdata) + geom_point(aes(x = pos,y = hapA_cov),color="brown3",size=0.5) + geom_point(aes(x = pos,y = hapB_cov),color="deepskyblue3",size=0.5) + scale_y_continuous(limits = ylimits) + scale_x_continuous(limits = xlimits)
p1 <- ggplot(data=binned) + geom_point(aes(x = middle_pos,y = meanAcov),color="brown3",size=0.2) + geom_point(aes(x = middle_pos,y = meanBcov),color="deepskyblue3",size=0.2) + scale_y_continuous(limits = ylimits) + scale_x_continuous(limits = xlimits)
p2 <- ggplot(data=blocks) + geom_segment(aes(x=min_pos,xend=max_pos,y=-0.0,yend=0.0),color= "grey25",size=8) + scale_x_continuous(limits = xlimits) + theme_minimal() + theme(axis.line=element_blank(),axis.text.x=element_blank(),axis.text.y=element_blank(),axis.ticks=element_blank(),axis.title.x=element_blank(),axis.title.y=element_blank())   #+ theme_void()      #+ theme_minimal()

p1 <- ggplot_gtable(ggplot_build(p1))   
p2 <- ggplot_gtable(ggplot_build(p2))
maxWidth = unit.pmax(p1$widths[2:3], p2$widths[2:3] )  
p1$widths[2:3] <- maxWidth              
p2$widths[2:3] <- maxWidth

grid.arrange(p1, p2, ncol = 1, heights = c(6, 0.4))  

#hap_file <- "./hap_gold2_solution_"
#blocks <- raw_data %>% group_by(block) %>% summarise(min_pos = min(pos),max_pos = max(pos))


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










#binbx_data <- binbx_data %>% mutate(normAcov = meanAcov/one_copy, normBcov = meanBcov/one_copy)


#xlimits <- c(7E7,8E7)

##############################################################
##############################################################
##############################################################

#raw_het_data <- paste(het_file_bx,chromosome,".dat",sep="") %>% read_delim(delim="\t",col_names = FALSE) %>% dplyr::rename("pos"="X2","tot_cov"="X10","ref_bx"="X11","var_bx"="X12")
#phets <- intersect(raw_het_data$pos,hdata$pos)
#raw_data_phets <- hdata %>% filter(pos %in% phets) %>% arrange(desc(pos))
#raw_het_data_phets <- raw_het_data %>% filter(pos %in% phets) %>% arrange(desc(pos)) %>% select(pos,ref_bx,var_bx,tot_cov)
#combined_het_data <- bind_cols(raw_data_phets,raw_het_data_phets)
#combined_het_data$ref_bx <- substring(combined_het_data$ref_bx,3)
#combined_het_data$var_bx <- substring(combined_het_data$var_bx,3)
#combined_het_data <- combined_het_data %>% filter(tot_cov > 10)
#combined_het_data <- combined_het_data %>% mutate(bxA = case_when(hap == 1 ~ ref_bx, hap == -1 ~ var_bx)) %>% mutate(bxB = case_when(hap == -1 ~ ref_bx, hap == 1 ~ var_bx))
#bx_binned <- combined_het_data %>% group_by(bin) %>% mutate(rbxA = paste0(bxA,collapse=","),rbxB = paste0(bxB,collapse=",")) %>% dplyr::summarise(n=n(),min_pos = min(pos),max_pos = max(pos),rbxA=first(rbxA),rbxB=first(rbxB),sumcov = sum(tot_cov)/n())  #meanAcov = mean(hapA_cov), meanBcov = mean(hapB_cov),

#bxA_unique = c(); bxB_unique = c(); bx_intersection = c(); bx_intersection_set = c()
#for(i in 1:nrow(bx_binned)) {
#  list_bxA <- str_split(bx_binned$rbxA[i],",") %>% unlist()
#  list_bxB <- str_split(bx_binned$rbxB[i],",") %>% unlist()
#  list_bxA <- list_bxA[list_bxA != ""]
#  list_bxB <- list_bxB[list_bxB != ""]
#  number_common_bx <- length( intersect(unique(list_bxA),unique(list_bxB)) )
#  bx_intersection <- append(bx_intersection,length( intersect(list_bxA,list_bxB) ))
#  bx_intersection_set <- append(bx_intersection_set,length(intersect( unique(list_bxA),unique(list_bxB) )))
#  bxA_unique <- append(bxA_unique,length(unique(list_bxA))-number_common_bx)
#  bxB_unique <- append(bxB_unique,length(unique(list_bxB))-number_common_bx)
#}
#unique_bx_binned <- bind_cols(bx_binned,tibble(bxA_unique),tibble(bxB_unique),tibble(bx_intersection),tibble(bx_intersection_set))
#unique_bx_binned <- unique_bx_binned %>% mutate(fracA = bxA_unique/( bxB_unique + bxA_unique )) %>% mutate(unique_meanAcov = fracA*sumcov,unique_meanBcov = (1-fracA)*sumcov)         #,tot_cov = meanAcov+meanBcov  #unique_bx_binned <- unique_bx_binned
#un_bin <- unique_bx_binned %>% select(n,bin,min_pos,max_pos,unique_meanAcov,unique_meanBcov,fracA,bxA_unique,bxB_unique,bx_intersection,bx_intersection_set)
#un_bin <- un_bin %>% filter(n>5)
##############################################################
##############################################################
##############################################################




#het_file_normal <- "./../phasing_HCC1954_linker_v1.8_illumina_PCRFree_tenx/output/gold_filtered_coverage_may23_normal_tenx_"
#one_copy <- 0.215  # tenx
#het_file_normal_bx <- "./../phasing_HCC1954_linker_v1.8_illumina_PCRFree_tenx/output/het_coverage_june18_normal_EA_tenx_"
#one_copy <- 0.205 # tenx_EA

##############################################################
##raw_het_data <- read_delim(file_string ,delim="\t",col_names = FALSE) %>% dplyr::rename("pos"="X2","tot_cov"="X10","ref_bx"="X11","var_bx"="X12")




#p1 <- ggplot(data=binned) + geom_point(aes(x = middle_pos,y = tumorAnorm),color="brown3",size=0.7) + 
#      geom_point(aes(x = middle_pos,y = tumorBnorm),color="deepskyblue3",size=0.7) +  
#      scale_x_continuous(limits = xlimits,position = "top") + 
#      theme_minimal() + xlab("position") + ylab("copy number") +
#      scale_y_continuous(limits = ylimits,breaks=c(0,1,2,3,4,5,6,8,10,12),label=c("0","1","2","3","4","5","6","8","10","12"),expand = c(0,0)) +
#      theme(axis.title.x=element_blank(),axis.text.x=element_blank())







######################################################################
#file_string <- "./../hg38_info/chromosome_banding_hg38.txt"
#hg38_data <- read_delim(file_string,delim="\t",col_names = TRUE) %>% dplyr::rename("chrom"=`#chrom`) %>% filter(chrom == chromosome)




#hg38_data <- read_delim("./../hg38_info/chromosome_banding_hg38.txt",delim="\t",col_names = TRUE) %>% dplyr::rename("chrom"=`#chrom`) %>% filter(chrom == chromosome)


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
#binned <- binned %>% mutate( tumorAnorm = tumorA/one_copy, tumorBnorm = tumorB/one_copy )
#binned <- binned %>% filter(n>10)
#binned <- binned %>% filter(n>10)




#het_file <- "./../phasing_HCC1954_linker_v1.8_illumina_PCRFree_tenx/output/het_coverage_may22_tumor_tenx_"
#het_file_normal <- "./../phasing_HCC1954_linker_v1.8_illumina_PCRFree_tenx/output/gold_filtered_coverage_may23_normal_tenx_"
#ne_copy <- 0.215  # tenx




#het_file <- "./../phasing_HCC1954_linker_v1.8_illumina_PCRFree_tenx/output/gold_filtered_coverage_june11_PCRFree_tumor_"
#het_file_normal <- "./../phasing_HCC1954_linker_v1.8_illumina_PCRFree_tenx/output/gold_filtered_coverage_june11_PCRFree_normal_"
#one_copy <- 0.16 # PCRFree
#het_file <- "./../phasing_HCC1954_linker_v1.8_illumina_PCRFree_tenx/output/gold_filtered_coverage_june11_illumina_tumor_"
#het_file_normal <- "./../phasing_HCC1954_linker_v1.8_illumina_PCRFree_tenx/output/gold_filtered_coverage_june11_illumina_normal_"
#one_copy <- 0.18 # illumina








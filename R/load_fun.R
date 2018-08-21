
#######################################################################################################################
load_hap <- function(hap_file,binsize,chromosome) {
  #### load haplotype
  #############################################################
  hap_data <- paste(hap_file,chromosome,".dat",sep="") %>% read_delim(delim="\t",col_names = FALSE)
  raw_data <- dplyr::rename(hap_data,"index"="X1","pos"="X2","ref"="X3","var"="X4","hap"="X5","hapA_cov"="X6","hapB_cov"="X7","block"="X10")
  hdata <- raw_data %>% separate(var,c("a","b","rbase","vbase"),"_",extra = "drop") %>% select("pos","rbase","vbase","hap","block")
  hdata <- hdata %>% mutate(bin = floor(pos/binsize))
  return(hdata)
}

#######################################################################################################################
load_hap_circos <- function(hap_file,binsize,chrom_list,one_copy) {
  #### load haplotype
  #############################################################
  hap_tracks = tibble()
  for (chr in chrom_list) {
    hdata <- load_hap(hap_file,binsize,chr)
    blocks <- hdata %>% group_by(block) %>% summarise(min_pos = min(pos),max_pos = max(pos))
    binned <- load_het_cov(het_file,het_file_normal,one_copy,chr,hdata)
    binned <- binned %>% filter(n>5)
    binned <- binned %>% mutate("chr" = chr);
    hap_tracks <- bind_rows(hap_tracks,binned)
  }
  return(hap_tracks)
}

#######################################################################################################################
load_hap_circos_tenx <- function(hap_file,binsize,chrom_list,one_copy) {
  #### load haplotype
  #############################################################
  hap_tracks = tibble()
  for (chr in chrom_list) {
    ##############################################################
    hdata <- load_hap(hap_file,binsize,chr)
    blocks <- hdata %>% group_by(block) %>% summarise(min_pos = min(pos),max_pos = max(pos))
    binned <- load_het_cov(het_file,het_file_normal,one_copy,chr,hdata)
    binned <- binned %>% filter(n>5)    #binned <- binned %>% filter(n>10)

    ##############################################################
    un_bin <- load_het_cov_tenx(het_file_bx,chr,hdata)     #,het_file_normal_bx
    un_bin <- un_bin %>% filter(n>5)

    binned_bx_unique <- inner_join(binned,un_bin,by = c("bin","bin"))
    #binned_bx_unique <- left_join(binned,un_bin,by = c("bin","bin"))
    #binned_bx_unique <- binned_bx_unique %>% mutate(tumorAnorm_unique = unique_meanAcov/one_copy,tumorBnorm_unique = unique_meanBcov/one_copy)

    binbx_data <- load_bx_coverage(bx_coverage_file,binsize,chr)
    binbx_data_normal <- load_bx_coverage(bx_coverage_normal_file,binsize,chr)

    ##############################################################
    binbx_data <- left_join(binbx_data,binbx_data_normal,by=c("bin"="bin"))
    binbx_data <- left_join(binned_bx_unique,binbx_data,by=c("bin"="bin"))
    binbx_data <- binbx_data %>% dplyr::rename("ubx"="ubx.x","normal_ubx"="ubx.y","totbx"="totbx.x","normal_totbx"="totbx.y")
    binbx_data <- binbx_data %>% mutate(ubx_newA = fracA*ubx,ubx_newB = (1-fracA)*ubx)
    binbx_data <- binbx_data %>% mutate(frac_ubx = ubx/normal_ubx) %>% mutate(ubx_meanAcov = fracA*frac_ubx,ubx_meanBcov = (1-fracA)*frac_ubx)
    binbx_data <- binbx_data %>% mutate(tumorAnorm_unique = ubx_meanAcov/one_copy,tumorBnorm_unique = ubx_meanBcov/one_copy)
    binbx_data <- binbx_data %>% mutate("chr" = chr);
    hap_tracks <- bind_rows(hap_tracks,binbx_data)
  }
  return(hap_tracks)
}

#######################################################################################################################
load_het_cov <- function(het_file,het_file_normal,one_copy,chromosome,hdata) {
  #### load tumor coverage
  #############################################################
  het_data <- paste(het_file,chromosome,".dat",sep="") %>% read_delim(delim="\t",col_names = FALSE)
  raw_data <- dplyr::rename(het_data,"pos"="X2","ref:var"="X3","ibase"="X4","dbase"="X5","gbase"="X6","cbase"="X7","abase"="X8","tbase"="X9","tot_cov"="X10")
  cdata <- raw_data %>% separate(col = dbase, into = c("d","dbase"), sep = "\\|") %>% separate(col = ibase, into = c("i","ibase"), sep = "\\|") %>% separate(col = gbase, into = c("g","gbase"), sep = "\\|") %>% separate(col = cbase, into = c("c","cbase"), sep = "\\|") %>% separate(col = abase, into = c("a","abase"), sep = "\\|") %>% separate(col = tbase, into = c("t","tbase"), sep = "\\|")
  cdata <- cdata %>% select("pos","ibase","dbase","gbase","cbase","abase","tbase","tot_cov")

  #### load normal coverage
  #############################################################
  het_data_normal <- paste(het_file_normal,chromosome,".dat",sep="") %>% read_delim(delim="\t",col_names = FALSE)
  raw_data_normal <- dplyr::rename(het_data_normal,"pos"="X2","ref:var"="X3","ibase"="X4","dbase"="X5","gbase"="X6","cbase"="X7","abase"="X8","tbase"="X9","tot_cov"="X10")
  cdata_normal <- raw_data_normal %>% separate(col = dbase, into = c("d","dbase"), sep = "\\|") %>% separate(col = ibase, into = c("i","ibase"), sep = "\\|") %>% separate(col = gbase, into = c("g","gbase"), sep = "\\|") %>% separate(col = cbase, into = c("c","cbase"), sep = "\\|") %>% separate(col = abase, into = c("a","abase"), sep = "\\|") %>% separate(col = tbase, into = c("t","tbase"), sep = "\\|")
  cdata_normal <- cdata_normal %>% select("pos","ibase","dbase","gbase","cbase","abase","tbase","tot_cov")

  #### combine haplotype and coverage
  #############################################################
  jdata_temp <- left_join(hdata,cdata,by = c("pos","pos"))
  jdata <- left_join(jdata_temp,cdata_normal,by = c("pos","pos"))
  jdata <- jdata %>% mutate(hapA_base = case_when(hap==1 ~ rbase, hap==-1 ~ vbase)) %>% mutate(hapB_base = case_when(hap==1 ~ vbase, hap==-1 ~ rbase))
  #############################################################
  jdata <- jdata %>% mutate(hapA_cov = case_when(hapA_base=="A" ~ as.numeric(abase.x), hapA_base=="C" ~ as.numeric(cbase.x), hapA_base=="G" ~ as.numeric(gbase.x),  hapA_base=="T" ~ as.numeric(tbase.x) ))
  jdata <- jdata %>% mutate(hapB_cov = case_when(hapB_base=="A" ~ as.numeric(abase.x), hapB_base=="C" ~ as.numeric(cbase.x), hapB_base=="G" ~ as.numeric(gbase.x),  hapB_base=="T" ~ as.numeric(tbase.x) ))
  jdata <- jdata %>% mutate(hapA_cov_normal = case_when(hapA_base=="A" ~ as.numeric(abase.y), hapA_base=="C" ~ as.numeric(cbase.y), hapA_base=="G" ~ as.numeric(gbase.y),  hapA_base=="T" ~ as.numeric(tbase.y) ))
  jdata <- jdata %>% mutate(hapB_cov_normal = case_when(hapB_base=="A" ~ as.numeric(abase.y), hapB_base=="C" ~ as.numeric(cbase.y), hapB_base=="G" ~ as.numeric(gbase.y),  hapB_base=="T" ~ as.numeric(tbase.y) ))
  jdata <- drop_na(jdata)

  #### bin data
  #############################################################
  binned <- jdata %>% group_by(bin) %>% dplyr::summarise(meanAcov = mean(hapA_cov), meanBcov = mean(hapB_cov), meantotcov = mean(tot_cov.x),meanAcov_normal = mean(hapA_cov_normal), meanBcov_normal = mean(hapB_cov_normal), meantotcov_norm = mean(tot_cov.y), n=n(),min_pos = min(pos),max_pos = max(pos))
  binned <- binned %>% mutate( tumorA = meanAcov/meantotcov_norm, tumorB = meanBcov/meantotcov_norm ) %>% mutate( middle_pos = ((min_pos + max_pos)/2.0) )
  binned <- binned %>% mutate( tumorAnorm = tumorA/one_copy, tumorBnorm = tumorB/one_copy )    #binned <- binned %>% filter(n>10)
  return(binned)
}

#######################################################################################################################
load_het_cov_tenx <- function(het_file_bx,chromosome,hdata) {
  raw_het_data <- paste(het_file_bx,chromosome,".dat",sep="") %>% read_delim(delim="\t",col_names = FALSE) %>% dplyr::rename("pos"="X2","tot_cov"="X10","ref_bx"="X11","var_bx"="X12")
  phets <- intersect(raw_het_data$pos,hdata$pos)
  raw_data_phets <- hdata %>% filter(pos %in% phets) %>% arrange(desc(pos))
  raw_het_data_phets <- raw_het_data %>% filter(pos %in% phets) %>% arrange(desc(pos)) %>% select(pos,ref_bx,var_bx,tot_cov)
  combined_het_data <- bind_cols(raw_data_phets,raw_het_data_phets)
  combined_het_data$ref_bx <- substring(combined_het_data$ref_bx,3)
  combined_het_data$var_bx <- substring(combined_het_data$var_bx,3)
  combined_het_data <- combined_het_data %>% filter(tot_cov > 10)
  combined_het_data <- combined_het_data %>% mutate(bxA = case_when(hap == 1 ~ ref_bx, hap == -1 ~ var_bx)) %>% mutate(bxB = case_when(hap == -1 ~ ref_bx, hap == 1 ~ var_bx))
  bx_binned <- combined_het_data %>% group_by(bin) %>% mutate(rbxA = paste0(bxA,collapse=","),rbxB = paste0(bxB,collapse=",")) %>% dplyr::summarise(n=n(),min_pos = min(pos),max_pos = max(pos),rbxA=first(rbxA),rbxB=first(rbxB),sumcov = sum(tot_cov)/n())  #meanAcov = mean(hapA_cov), meanBcov = mean(hapB_cov),

  bxA_unique = c(); bxB_unique = c(); bx_intersection = c(); bx_intersection_set = c()
  for(i in 1:nrow(bx_binned)) {
    list_bxA <- strsplit(bx_binned$rbxA[i],",") %>% unlist()
    list_bxB <- strsplit(bx_binned$rbxB[i],",") %>% unlist()
    #list_bxA <- str_split(bx_binned$rbxA[i],",") %>% unlist()
    #list_bxB <- str_split(bx_binned$rbxB[i],",") %>% unlist()
    list_bxA <- list_bxA[list_bxA != ""]
    list_bxB <- list_bxB[list_bxB != ""]
    number_common_bx <- length( intersect(unique(list_bxA),unique(list_bxB)) )
    bx_intersection <- append(bx_intersection,length( intersect(list_bxA,list_bxB) ))
    bx_intersection_set <- append(bx_intersection_set,length(intersect( unique(list_bxA),unique(list_bxB) )))
    bxA_unique <- append(bxA_unique,length(unique(list_bxA))-number_common_bx)
    bxB_unique <- append(bxB_unique,length(unique(list_bxB))-number_common_bx)
  }
  unique_bx_binned <- bind_cols(bx_binned,tibble(bxA_unique),tibble(bxB_unique),tibble(bx_intersection),tibble(bx_intersection_set))
  unique_bx_binned <- unique_bx_binned %>% mutate(fracA = bxA_unique/( bxB_unique + bxA_unique )) %>% mutate(unique_meanAcov = fracA*sumcov,unique_meanBcov = (1-fracA)*sumcov)         #,tot_cov = meanAcov+meanBcov  #unique_bx_binned <- unique_bx_binned
  un_bin <- unique_bx_binned %>% select(n,bin,min_pos,max_pos,unique_meanAcov,unique_meanBcov,fracA,bxA_unique,bxB_unique,bx_intersection,bx_intersection_set)
  return(un_bin)
}

#######################################################################################################################
load_hg38_data <- function() {
  #############################################################
  file_string <- "./../hg38_info/chromosome_banding_hg38.txt"
  hg38_data <- read_delim(file_string,delim="\t",col_names = TRUE) %>% dplyr::rename("chrom"=`#chrom`) %>% filter(chrom == chromosome)
  return(hg38_data)
}

#######################################################################################################################
load_phasedsv_file <- function(phased_file,chromosome) {
  sv_out <- list()
  phased_data <- read_delim(phased_file,delim="\t",col_names = FALSE) %>% dplyr::rename("chr1"="X1","pos1"="X2","chr2"="X3","pos2"="X4","hap1"="X5","hap2"="X6","sv_type"="X7")  #phased_data,

  #############################################################
  intra <- phased_data %>% filter(chr1 == chromosome & chr2 == chromosome)
  inter <- phased_data %>% filter(chr1 == chromosome | chr2 == chromosome) %>% filter( chr2 != chr1 )

  inter <- inter %>% mutate(single_chr_pos = case_when(chr1==chromosome ~ pos1, chr2==chromosome ~ pos2))
  inter <- inter %>% mutate(sv_to = case_when(chr1==chromosome ~ chr2, chr2==chromosome ~ chr1))
  inter <- inter %>% mutate(single_hap = case_when(chr1==chromosome ~ hap1, chr2==chromosome ~ hap2))

  sv_out$null_inter <- inter %>% filter(single_hap == 0)
  sv_out$hapA_inter <- inter %>% filter(single_hap == 1)
  sv_out$hapB_inter <- inter %>% filter(single_hap == -1)

  #############################################################
  intra <- intra %>% mutate(type = hap1*hap2)
  intra <- intra %>% mutate(type = replace(type,which(type == -1),3))
  intra <- intra %>% mutate(type = replace(type,which(type == 1 & hap2 == -1),1))
  intra <- intra %>% mutate(type = replace(type,which(type == 1 & hap2 == 1),2))

  phasedb_links <- intra %>% filter( type == 1 )
  phaseda_links <- intra %>% filter( type == 2 )
  mix_links <- intra %>% filter( type == 3 )
  null_links <- intra %>% filter( type == 0 )

  sv_out$phaseda_links <- phaseda_links %>% add_row(chr1 = chromosome, pos1 = -1000, chr2 = chromosome, pos2 = -1001, hap1=1, hap2=1, sv_type = "intra", type = 2)
  sv_out$phasedb_links <- phasedb_links %>% add_row(chr1 = chromosome, pos1 = -1000, chr2 = chromosome, pos2 = -1001, hap1=-1, hap2=-1, sv_type = "intra", type = 1)
  sv_out$mix_links <- mix_links %>% add_row(chr1 = chromosome, pos1 = -1000, chr2 = chromosome, pos2 = -1001, hap1=1, hap2=-1, sv_type = "intra", type = 3)
  sv_out$null_links <- null_links %>% add_row(chr1 = chromosome, pos1 = -1000, chr2 = chromosome, pos2 = -1001, hap1=0, hap2=-1, sv_type = "intra", type = 0)
  return(sv_out)
}

#######################################################################################################################
load_bx_coverage <- function(bx_coverage_file,binsize,chromosome) {
  binbx_data <- paste(bx_coverage_file,chromosome,".dat",sep="") %>% read_delim(delim="\t",col_names = FALSE) %>% dplyr::rename("chr"="X1","pos"="X2","totbx"="X3","ubx"="X4") %>% mutate(bin = pos/binsize)
  return(binbx_data)
}

#######################################################################################################################
load_phasedsv_file_total <- function(phased_file,chrom_list,band_width) {
  sv_out <- list()
  phased_data <- read_delim(phased_file,delim="\t",col_names = FALSE) #%>% filter( chr1 %in% chrom_list ) %>% filter( chr2 %in% chrom_list ) #%>% filter( TotalCount > 10 )
  phased_data <- rename(phased_data,"chr1"="X1","pos1"="X2","chr2"="X3","pos2"="X4","hap1"="X5","hap2"="X6","call_type"="X7")
  phased_data <- phased_data %>% filter( chr1 %in% chrom_list ) %>% filter( chr2 %in% chrom_list )
  phased_data <- phased_data %>% mutate(type = hap1*hap2)
  phased_data <- phased_data %>% mutate(type = replace(type,which(type == -1),3))
  phased_data <- phased_data %>% mutate(type = replace(type,which(type == 1 & hap2 == -1),1))
  phased_data <- phased_data %>% mutate(type = replace(type,which(type == 1 & hap2 == 1),2))

  sv_out$phasedb_links <- phased_data %>% filter(type == 1)
  sv_out$phaseda_links <- phased_data %>% filter(type == 2)
  sv_out$mix_links <- phased_data %>% filter(type == 3)
  sv_out$null_links <- phased_data %>% filter(type == 0)

  sv_out$b_link1 <- sv_out$phasedb_links %>% select(chr1,pos1,call_type) %>% mutate(end = as.integer(pos1 + band_width)) %>% rename(start=pos1,chr=chr1) %>% as.data.frame()
  sv_out$b_link2 <- sv_out$phasedb_links %>% select(chr2,pos2,call_type) %>% mutate(end = as.integer(pos2 + band_width)) %>% rename(start=pos2,chr=chr2) %>% as.data.frame()
  sv_out$a_link1 <- sv_out$phaseda_links %>% select(chr1,pos1,call_type) %>% mutate(end = as.integer(pos1 + band_width)) %>% rename(start=pos1,chr=chr1) %>% as.data.frame()
  sv_out$a_link2 <- sv_out$phaseda_links %>% select(chr2,pos2,call_type) %>% mutate(end = as.integer(pos2 + band_width)) %>% rename(start=pos2,chr=chr2) %>% as.data.frame()
  sv_out$m_link1 <- sv_out$mix_links %>% select(chr1,pos1,call_type) %>% mutate(end = as.integer(pos1 + band_width)) %>% rename(start=pos1,chr=chr1) %>% as.data.frame()
  sv_out$m_link2 <- sv_out$mix_links %>% select(chr2,pos2,call_type) %>% mutate(end = as.integer(pos2 + band_width)) %>% rename(start=pos2,chr=chr2) %>% as.data.frame()
  sv_out$n_link1 <- sv_out$null_links %>% select(chr1,pos1,call_type) %>% mutate(end = as.integer(pos1 + band_width)) %>% rename(start=pos1,chr=chr1) %>% as.data.frame()
  sv_out$n_link2 <- sv_out$null_links %>% select(chr2,pos2,call_type) %>% mutate(end = as.integer(pos2 + band_width)) %>% rename(start=pos2,chr=chr2) %>% as.data.frame()
  return(sv_out)
}

#######################################################################################################################
load_phasedsv_positions <- function(phased_file,chrom_list,band_width) {
  sv_out <- list()
  phased_data <- read_delim(phased_file,delim="\t",col_names = FALSE) #%>% filter( chr1 %in% chrom_list ) %>% filter( chr2 %in% chrom_list ) #%>% filter( TotalCount > 10 )
  phased_data <- rename(phased_data,"chr1"="X1","pos1"="X2","chr2"="X3","pos2"="X4","hap1"="X5","hap2"="X6","call_type"="X7")
  #phased_data <- phased_data %>% filter( chr1 %in% chrom_list ) %>% filter( chr2 %in% chrom_list )
  anchor1 <- phased_data %>% select(chr1,pos1,hap1) %>% rename(chr=chr1,min_pos=pos1,value1=hap1)
  anchor2 <- phased_data %>% select(chr2,pos2,hap2) %>% rename(chr=chr2,min_pos=pos2,value1=hap2)
  anchor <- bind_rows(anchor1,anchor2)
  anchor <- anchor %>% filter( chr %in% chrom_list )
  return(anchor)
}







  #file_string <- file_string ,
  #phased_data <- dplyr::rename(phased_data,"chr1"="X1","pos1"="X2","chr2"="X3","pos2"="X4","hap1"="X5","hap2"="X6","sv_type"="X7")


# run_full_report <- function(ptrial) {
#   if(!exists('workingDir')) {workingDir <- "/home/ubuntu/data/CART19/CART19_from_git2"}
#   inputFile <- file.path(workingDir,'ONLY_ALL',"04_report.Rmd")
#   outFile <- file.path(workingDir,'ONLY_ALL','goi_reports',
#                        paste(ptrial,'_goi_', Sys.Date(), '.pdf', sep=''))
#   rmarkdown::render( 
#     input       = inputFile, 
#     #encoding    = encoding, 
#     params      = list(pTRIAL=ptrial,
#                        tPART=ifelse(ptrial=='ALL','ALL',
#                              ifelse(ptrial=='CLL','CLL',
#                              ifelse(ptrial=='CALL','ALL and CLL','none')))),      
#     output_file = outFile)
#   return(TRUE)
# }
# run_full_report('ALL')
# run_full_report('CLL')


#### 
library(rlang)
library(tidyverse)
if(!exists('workingDir')) {workingDir <- "/home/ubuntu/data/CART19/CART19_from_git2"}

to_run <- expand_grid(x= c('ALL','CLL','CALL'),y= c('RE','na')) 

tt<- to_run %>%   
purrr::pmap(~source(file.path(workingDir,'01_condensed_intsites.R'),local=env(TRIAL=.x,RESP=.y)))
  


# cluster reports -------------------
run_cluster_report <- function(ptrial,group_c) {
  if(!exists('workingDir')) {workingDir <- "/home/ubuntu/data/CART19/CART19_from_git2"}
  inputFile <- file.path(workingDir,'ONLY_ALL',"06_write_cluster.Rmd")
  outFile <- file.path(workingDir,'ONLY_ALL','cluster_reports',
                       paste(ptrial,'_Scan_stats_',group_c,'_', Sys.Date(), '.pdf', sep=''))
  rmarkdown::render( 
    input       = inputFile, 
    #encoding    = encoding, 
    params      = list(pTRIAL=ptrial,
                       count=group_c,
                       tPART=ifelse(ptrial=='ALL','ALL',
                                         ifelse(ptrial=='CLL','CLL',
                                         ifelse(ptrial=='CALL','ALL and CLL','none')))),      
    output_file = outFile)
  return(TRUE)
}

run_cluster_report('ALL' ,'single')
run_cluster_report('CLL' ,'single')
run_cluster_report('CALL','single')

#run_cluster_report('ALL' ,'multi')
#run_cluster_report('CLL' ,'multi')
#run_cluster_report('CALL','multi')



### generate bed files ----------------------------
if(!exists('workingDir')) {workingDir <- "/home/ubuntu/data/CART19/CART19_from_git2"}
library(tidyverse)
library(GenomicRanges)
library(rtracklayer)
#function to format metadata for bigbed export
fix_to_bigbed <- function(in_bed) {
  hhh3 <- in_bed
  hhh3$timepoint <- NULL
  
  temp_meta <- list(score=hhh3$score,
                    name=hhh3$name,
                    itemRgb=hhh3$itemRgb)
  
  hhh3$score <- NULL
  hhh3$name <- NULL
  hhh3$itemRgb <- NULL
  
  hhh3$name <- temp_meta$name
  hhh3$score <- temp_meta$score
  hhh3$thick<- ranges(hhh3)
  hhh3$itemRgb <- temp_meta$itemRgb
  return(hhh3)
}



produce_bed <- function(TRIAL) {
  cond_uniq_sites_tmp <- readRDS(file.path(workingDir,paste0('ONLY_',TRIAL),paste0("only_",TRIAL,"_condensed_intsites.rds")))
  save_seqinfo <- seqinfo(cond_uniq_sites_tmp)
  cond_uniq_sites <- cond_uniq_sites_tmp %>% 
    as.data.frame() %>%
  #  head() %>% 
    select(c(seqnames,start,end,strand,estAbund,timepoint)) %>%
    mutate(name=paste0('item_',dplyr::row_number())) %>% 
    dplyr::rename('score'=estAbund) %>%
    mutate(itemRgb=ifelse(strand=='+',c('#ff0000'),c('#0000ff'))) %>% 
  #  mutate(estAbund=NULL) %>% 
    GenomicRanges::makeGRangesFromDataFrame(keep.extra.columns = T,seqinfo = save_seqinfo)

#cond_uniq_sites$score <- cond_uniq_sites$estAbund
  tdn_sites <- cond_uniq_sites[cond_uniq_sites$timepoint == "d0"]
  timepoint_sites <- cond_uniq_sites[cond_uniq_sites$timepoint != "d0"]
  
  export(fix_to_bigbed(timepoint_sites),file.path(workingDir,'BED','hg38',paste0(TRIAL,'_tp.bb')),'bb')
  export(fix_to_bigbed(tdn_sites),file.path(workingDir,'BED','hg38',paste0(TRIAL,'_tdn.bb')),'bb')
  
}

produce_bed('ALL')
produce_bed('CLL')
  



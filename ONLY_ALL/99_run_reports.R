

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

run_cluster_report('ALL' ,'multi')
run_cluster_report('CLL' ,'multi')
run_cluster_report('CALL','multi')



### generate bed files
if(!exists('workingDir')) {workingDir <- "/home/ubuntu/data/CART19/CART19_from_git2"}
library(tidyverse)
library(GenomicRanges)
TRIAL <-'ALL'
cond_uniq_sites_tmp <- readRDS(file.path(workingDir,paste0('ONLY_',TRIAL),paste0("only_",TRIAL,"_condensed_intsites.rds")))
save_seqinfo <- seqinfo(cond_uniq_sites_tmp)
cond_uniq_sites <- cond_uniq_sites_tmp %>% 
  as.data.frame() %>%
  head() %>% 
  select(c(seqnames,start,end,strand,estAbund)) %>%
  dplyr::rename('score'=estAbund) %>% 
  GenomicRanges::makeGRangesFromDataFrame(keep.extra.columns = T,seqinfo = save_seqinfo)

#cond_uniq_sites$score <- cond_uniq_sites$estAbund
tdn_sites <- cond_uniq_sites[cond_uniq_sites$timepoint == "d0"]
timepoint_sites <- cond_uniq_sites[cond_uniq_sites$timepoint != "d0"]

library(rtracklayer)
export(timepoint_sites,file.path(workingDir,'BED','hg38','ALL_tp.bb'),'bb')
export(tdn_sites,file.path(workingDir,'BED','hg38','ALL_tdn.bb'),'bb')



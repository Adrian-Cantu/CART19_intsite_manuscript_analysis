

run_cluster_report <- function(ptrial) {
  if(!exists('workingDir')) {workingDir <- "/home/ubuntu/data/CART19/CART19_from_git2"}
  inputFile <- file.path(workingDir,'ONLY_ALL',"06_write_cluster.Rmd")
  outFile <- file.path(workingDir,'ONLY_ALL','cluster_reports',paste(ptrial,'_Scan_stats_', Sys.Date(), '.pdf', sep=''))
  rmarkdown::render( 
    input       = inputFile, 
    #encoding    = encoding, 
    params      = list(pTRIAL=ptrial,tPART=ifelse(ptrial=='ALL','ALL',
                                         ifelse(ptrial=='CLL','CLL',
                                         ifelse(ptrial=='CALL','ALL and CLL','none')))),      
    output_file = outFile)
  return(TRUE)
}

run_cluster_report('ALL')
run_cluster_report('CLL')
run_cluster_report('CALL')



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

library(tidyverse)
library(RMySQL)
dbConn  <- dbConnect(MySQL(), group='specimen_management')
samples_tmp <- dbGetQuery(dbConn, 'select * from gtsp where Trial="CART19_ALL"') 
#samples_tmp %>% dplyr::select(c(SpecimenAccNum,Timepoint)) %>% filter(SpecimenAccNum %in% isnatime)
(samples_tmp %>% filter(SpecimenAccNum=='GTSP4276'))$Timepoint
(samples_tmp %>% filter(SpecimenAccNum=='GTSP4286'))$Timepoint

ll <- rep('NA',nrow(samples_tmp))
samples_good <- data.frame(
  GTSP=samples_tmp$SpecimenAccNum,
  Trial=samples_tmp$Trial,
  Sampe_Source=ll,
  Abv_Cell_Type=samples_tmp$CellType,
  Sorting_Parameters = samples_tmp$SampleCellType,
  Additional_specimen_information=samples_tmp$SpecimenInfo,
  Days_Post_Treatment=samples_tmp$DaysPostTrtmnt,
  Timepoint=samples_tmp$Timepoint,
  Sample_Patient_Code=samples_tmp$SamplePatientCode,
  Patient_ID=samples_tmp$Patient,
  Disease=rep('ALL',nrow(samples_tmp)),
  Response_to_Treatment=ll,
  TransGene_Expression=ll,
  Date_Received=samples_tmp$DateRcvd,
  Completed=rep('yes',nrow(samples_tmp)),
  Set=ll,
  clin_trial=ll
)
write_csv(samples_good,'cart19_intsite_sample_list.csv')


bsession <- makeUCSCsession('hg18')
tableNames(ucscTableQuery(bsession, track = "CpG Islands"))
ucscTables('hg19',"CpG Islands")

ucscTableQuery(bsession, track = "CpG Islands", table = "cpgIslandExt")
getTable(ucscTableQuery(bsession, track = "CpG Islands", table = "cpgIslandExt"))

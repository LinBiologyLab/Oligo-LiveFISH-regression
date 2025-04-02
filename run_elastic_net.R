
selY <- "SNR"
cvseed <- 94
source('elastic_net.R')

selY <- "SBN"
cvseed <- 900
source('elastic_net.R')

selY <- "detection_SNR_greater_than_4"
cvseed <- 66
source('elastic_net.R')

selY <- "detection_SNR_greater_than_5"
cvseed <- 94
source('elastic_net.R')

selY <- "Detection_Ratios_SBN_greater_than_5"
cvseed <- 7307
source('elastic_net.R')


selY <- "SNR"
cvseed <- 4533
source('elastic_net_allvar.R')

selY <- "SBN"
cvseed <- 567
source('elastic_net_allvar.R')

selY <- "detection_SNR_greater_than_4"
cvseed <- 9081
source('elastic_net_allvar.R')

selY <- "detection_SNR_greater_than_5"
cvseed <- 871
source('elastic_net_allvar.R')

selY <- "Detection_Ratios_SBN_greater_than_5"
cvseed <- 4896
source('elastic_net_allvar.R')


source("elastic_net_allvar_cv.R")

sessionInfo() 


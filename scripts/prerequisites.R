#! R
##install prerequisite libraries
libs <- c("RnBeads.hg19","RnBeads","RPMM","tidyverse","magrittr",	"pheatmap","mgcv","readxl","GenomicRanges","reshape2","ggplot2","VennDiagram", "ggdendro")

##install the great installer
install.packages("BiocManager",
                  repos="https://cloud.r-project.org",
                  dependencies=TRUE)
library("BiocManager")

##install libs through BiocManager
lapply(libs, function(f){
  BiocManager::install(f,
                       site_repository="https://cloud.r-project.org",
                       ask=FALSE,
                       update=TRUE,
                       dependencies=TRUE)
  })

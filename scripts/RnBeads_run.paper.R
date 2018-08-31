#! /usr/local/bin/R

#source -> "http://rnbeads.mpi-inf.mpg.de/install.R"
libs <- c("RnBeads", "RnBeads.hg19", "RPMM", "tidyverse")
print("Loading libraries...")
libsLoaded <- lapply(libs,function(l){print(l);suppressWarnings(suppressMessages(library(l, character.only = TRUE)))})

#input arguments
argsin <- commandArgs(trailingOnly=TRUE)
print(argsin)

##data is held in [1]; same dir for output, but each gets a dir
BASEDIR <- argsin[1]
IDATDIR <- argsin[2]
TAG <- argsin[3]
ANNOFILE <- argsin[4]

##prepared earlier
sample_anno <- read_csv(ANNOFILE)

#set options
rnb.options(analysis.name=TAG,
  logging=TRUE,
  min.group.size=3,
  qc=TRUE,
  qc.boxplots=TRUE,
  qc.barplots=TRUE,
  qc.snp.boxplot=TRUE,
  filtering.greedycut=TRUE,
  differential=TRUE,
  differential.enrichment=TRUE,
  differential.site.test.method="limma",
  differential.comparison.columns="Group",
  filtering.sex.chromosomes.removal=TRUE,
  filtering.snp="yes",
  identifiers.column="Sample_ID",
  import.table.separator=",",
  disk.dump.big.matrices=FALSE,
  export.to.csv=TRUE,
  import=TRUE,
  exploratory = TRUE,
  exploratory.columns = "Group",
  exploratory.principal.components = 4,
  exploratory.deviation.plots = TRUE,
  exploratory.clustering.top.sites = 100,
  exploratory.clustering.heatmaps.pdf = TRUE,
  normalization.method="bmiq",
  normalization.background.method="none")

reports.dir <- paste0(BASEDIR, "/analysis/", TAG)
rnb.run.analysis(dir.reports=reports.dir,
                 sample.sheet=as.data.frame(sample_anno),
                 data.dir=IDATDIR,
                 data.type="infinium.idat.dir")

##create tibble of all samples per site:
load(paste0(reports.dir, "/rnbSet_preprocessed/rnb.set.RData"))
rnb.set <- object
grps <- factor(pheno(rnb.set)[,"Sample_ID"])
metht <- as_tibble(meth(rnb.set, "sites", row.names=TRUE)) %>%
         round(.,3)

##add sample IDs, probe IDs, round to 3 sig. digits, arrange, save
colnames(metht) <- grps
metht <- metht %>% mutate(rowj = seq(from=1, to=dim(metht)[1], by=1))
probes <- tibble(rowj = seq(from=1, to=length(rownames(rnb.set@sites)), by=1),
                 probe = rownames(rnb.set@sites))
betasites <- left_join(probes, metht, by=c("rowj", "rowj")) %>%
             dplyr::select(-rowj) %>%
             dplyr::arrange(probe)
write_tsv(x=as.data.frame(betasites),
          path=paste0(reports.dir, "/", TAG, ".methylation_per_sample.sites.tsv"),
          col_names=TRUE)

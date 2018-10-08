#! /usr/local/bin/R

libs <- c("RnBeads", "tidyverse", "VennDiagram")
print("Loading libraries...")
libsLoaded <- lapply(libs,function(l){suppressWarnings(suppressMessages(library(l, character.only = TRUE)))})

#input arguments
argsin <- commandArgs(trailingOnly=TRUE)
print(argsin)

##inputs
BASEDIR <- argsin[1]
SCRIPTDIR <- argsin[2]
TAG <- argsin[3]
EPICBED <- argsin[4]
BETAFILE <- argsin[5]

#source functions
source(paste0(SCRIPTDIR, "/set_patient_EMR.func.R"))

##identifying EMRs within each patient’s tumor tissue
##proportions of those detected in the patient’s plasma and urine sediment
# data.dir <- "/Volumes/GoogleDrive/My Drive/Romina/infiniumEPIC_0916/1359_Perry_EPIC_data/ANALYSIS_RnBeads/paper/HyperHypo_20pc/split_hyper-hypo"
# epic.file <- "/Volumes/GoogleDrive/My Drive/Romina/infiniumEPIC_0916/1359_Perry_EPIC_data/ANALYSIS_RnBeads/EPIC_annotation/MethylationEPIC_v-1-0_B4.Legacy_B2.bed"
epic.bed <- read_tsv(EPICBED, col_names=F, col_types="ciicc")
colnames(epic.bed) <- c("seqnames", "start", "end", "strand","anno")
epic.bed$probe <- suppressWarnings(separate(epic.bed, anno, "probe", ";")$probe)

##revision: load betasites directly
betasitet <- read_tsv(BETAFILE)

##set output directories
OUTDIR <- paste0(BASEDIR, "/analysis/probe_EMR")
probe.dir <- paste0(OUTDIR, "/probes")
emr.dir <- paste0(OUTDIR, "/EMRs")
pat.probe.dir <-  paste0(OUTDIR, "/probes/patient")
set.probe.dir <- paste0(OUTDIR, "/probes/set")
pat.emr.dir <- paste0(OUTDIR, "/EMRs/patient")
set.emr.dir <- paste0(OUTDIR, "/EMRs/set")

lapply(c(OUTDIR, probe.dir, emr.dir, pat.probe.dir, set.probe.dir, pat.emr.dir, set.emr.dir), function(f){
  dir.create(f, showWarnings=FALSE, recursive=TRUE)
})

if(file.exists(paste0(OUTDIR, "/set_patient-probe_EMRs.GRanges.RData"))){
  print(paste0("Found: ", OUTDIR, "/set_patient-probe_EMRs.GRanges.RData"))
  print("Loading...")
  load(paste0(OUTDIR, "/set_patient-probe_EMRs.GRanges.RData"))
  print("Load complete")
}

if(!file.exists(paste0(OUTDIR, "/set_patient-probe_EMRs.GRanges.RData"))){

  print(paste0("No file: ", OUTDIR, "/set_patient-probe_EMRs.GRanges.RData"))
  print("Making...")

  ##define levels on which to work: set, patient
  sampleSet <- c("Benign", "Tumour", "Plasma", "Urine")
  patients <- colnames(betasitet)[grep("probe", colnames(betasitet), invert=TRUE)]
  patGroups <- unique(unlist(lapply(patients,function(f){
    strsplit(f,"_")[[1]][[3]]
  })))

  ############
  ## PROBES ##
  ############

  ##per patient
  patProbeList <- probeHyperHypoLister(betavalues=betasitet, samples=patients, filters=c(">0.8","<0.2"))
  names(patProbeList) <- patients

  #write output, make venn4 and store that data
  patGroupVennProbeList <- as.list(patGroups)

  for(pg in 1:length(patGroups)){

    #individual dir for each patient
    pat.out.dir <- paste0(pat.probe.dir, "/patient_", patGroups[pg])
    dir.create(pat.out.dir, showWarnings=F)

    #setup
    patSet <- grep(patGroups[pg], patients)
    patGroupVennProbeList[[pg]] <- list("hyper", "hypo")

    #write list of patient hyper, hypo probes
    for(pt in 1:length(grep(patGroups[pg], patients))){
      print(paste0("Working on: ", patients[patSet[pt]]))
      write.table(patProbeList[[patients[patSet[pt]]]]$hyper$probe,
                  file=as.character(paste0(pat.out.dir,"/", patients[patSet[pt]], ".hyper.probes.txt")),
                  quote=F, sep="\t", col=F, row=F)
      write.table(patProbeList[[patients[patSet[pt]]]]$hypo$probe,
                  file=as.character(paste0(pat.out.dir,"/", patients[patSet[pt]], ".hypo.probes.txt")),
                  quote=F, sep="\t", col=F, row=F)
    }

    #individual probe sets
    benign_hyper_probes <- patProbeList[[patients[patSet[1]]]]$hyper$probe
    benign_hypo_probes <- patProbeList[[patients[patSet[1]]]]$hypo$probe
    tumour_hyper_probes <- patProbeList[[patients[patSet[2]]]]$hyper$probe
    tumour_hypo_probes <- patProbeList[[patients[patSet[2]]]]$hypo$probe
    plasma_hyper_probes <- patProbeList[[patients[patSet[3]]]]$hyper$probe
    plasma_hypo_probes <- patProbeList[[patients[patSet[3]]]]$hypo$probe
    urine_hyper_probes <- patProbeList[[patients[patSet[4]]]]$hyper$probe
    urine_hypo_probes <- patProbeList[[patients[patSet[4]]]]$hypo$probe

    setwd(pat.out.dir)
    patGroupVennProbeList[[pg]][[1]] <- venn4ProbeList(A=benign_hyper_probes,
          B=tumour_hyper_probes,
          C=plasma_hyper_probes,
          D=urine_hyper_probes,
          sampVec=c("Matched Normal", "Tumor", "Plasma", "Urine"),
          tag="Hyper Probes",
          nam=paste0(pat.out.dir, "/patient_", patGroups[pg], ".hyper"),
          fil=c("forestgreen","black","darkred","yellow"))
    patGroupVennProbeList[[pg]][[2]] <- venn4ProbeList(A=benign_hypo_probes,
          B=tumour_hypo_probes,
          C=plasma_hypo_probes,
          D=urine_hypo_probes,
          sampVec=c("Matched Normal", "Tumor", "Plasma", "Urine"),
          tag="Hypo Probes",
          nam=paste0(pat.out.dir, "/patient_", patGroups[pg], ".hypo"),
          fil=c("forestgreen","black","darkred","yellow"))
   names(patGroupVennProbeList[[pg]]) <- c("hyper", "hypo")

   ##distribution of T, U, S probes, not in benign
   #per tumour, plasma, urine
   selpg <- c("probe", grep("Benign", grep(paste0("_", pg), names(betasitet), value=TRUE), invert=T, value=T))
   selpg <- match(selpg,names(betasitet))
   tprobens <- match(c("n2","n23","n24","n234"), names(patGroupVennProbeList[[pg]][["hyper"]]))
   tprobeu <- unique(unlist(lapply(tprobens,function(f){
     hrs <- patGroupVennProbeList[[pg]][["hyper"]][[f]]
     hos <- patGroupVennProbeList[[pg]][["hypo"]][[f]]
     return(list(hrs,hos))
   })))
   tbetasitet <- betasitet %>% dplyr::filter(probe %in% tprobeu) %>%
                               dplyr::select(selpg)

   pprobens <- match(c("n3","n23","n34","n234"), names(patGroupVennProbeList[[pg]][["hyper"]]))
   pprobeu <- unique(unlist(lapply(pprobens,function(f){
     hrs <- patGroupVennProbeList[[pg]][["hyper"]][[f]]
     hos <- patGroupVennProbeList[[pg]][["hypo"]][[f]]
     return(list(hrs,hos))
   })))
   pbetasitet <- betasitet %>% dplyr::filter(probe %in% pprobeu) %>%
                               dplyr::select(selpg)

   uprobens <- match(c("n4","n24","n34","n234"), names(patGroupVennProbeList[[pg]][["hyper"]]))
   uprobeu <- unique(unlist(lapply(uprobens,function(f){
     hrs <- patGroupVennProbeList[[pg]][["hyper"]][[f]]
     hos <- patGroupVennProbeList[[pg]][["hypo"]][[f]]
     return(list(hrs,hos))
   })))
   ubetasitet <- betasitet %>% dplyr::filter(probe %in% uprobeu) %>%
                               dplyr::select(selpg)

   kdensPlot(meltdf=melt(tbetasitet[2:4]), pat=paste0("patient_", pg), nam="Tumor 20%", colz=c("black","red","yellow"))
   kdensPlot(meltdf=melt(pbetasitet[2:4]), pat=paste0("patient_", pg), nam="Plasma 20%", colz=c("black","red","yellow"))
   kdensPlot(meltdf=melt(ubetasitet[2:4]), pat=paste0("patient_", pg), nam="Urine 20%", colz=c("black","red","yellow"))
  }

  #set level
  setProbeList <- probeHyperHypoLister(betavalues=betasitet, samples=sampleSet, exclude="_1", filters=c(">0.8","<0.2"))
  names(setProbeList) <- sampleSet
  setProbeVennList <- list("hyper", "hypo")

  #write output
  for(st in 1:length(sampleSet)){
    print(paste0("Working on: ", sampleSet[st]))

    write.table(setProbeList[[sampleSet[st]]]$hyper$probe,
                file=as.character(paste0(set.probe.dir,"/", sampleSet[st], ".hyper.probes.txt")),
                quote=F, sep="\t", col=F, row=F)
    write.table(setProbeList[[sampleSet[st]]]$hypo$probe,
                file=as.character(paste0(set.probe.dir,"/", sampleSet[st], ".hypo.probes.txt")),
                quote=F, sep="\t", col=F, row=F)
  }

  benign_hyper_probes <- setProbeList[["Benign"]]$hyper$probe
  benign_hypo_probes <- setProbeList[["Benign"]]$hypo$probe
  tumour_hyper_probes <- setProbeList[["Tumour"]]$hyper$probe
  tumour_hypo_probes <- setProbeList[["Tumour"]]$hypo$probe
  plasma_hyper_probes <- setProbeList[["Plasma"]]$hyper$probe
  plasma_hypo_probes <- setProbeList[["Plasma"]]$hypo$probe
  urine_hyper_probes <- setProbeList[["Urine"]]$hyper$probe
  urine_hypo_probes <- setProbeList[["Urine"]]$hypo$probe

  setwd(set.probe.dir)

  setProbeVennList[[1]] <- venn4ProbeList(A=benign_hyper_probes,
        B=tumour_hyper_probes,
        C=plasma_hyper_probes,
        D=urine_hyper_probes,
        nam=paste0(set.probe.dir,"/set.hyper"),
        sampVec=c("Matched Normal", "Tumor", "Plasma", "Urine"),
        tag="Hyper Probes",
        fil=c("forestgreen","black","darkred","yellow"))
  setProbeVennList[[2]] <- venn4ProbeList(A=benign_hypo_probes,
        B=tumour_hypo_probes,
        C=plasma_hypo_probes,
        D=urine_hypo_probes,
        nam=paste0(set.probe.dir,"/set.hypo"),
        sampVec=c("Matched Normal", "Tumor", "Plasma", "Urine"),
        tag="Hypo Probes",
        fil=c("forestgreen","black","darkred","yellow"))
  names(setProbeVennList) <- c("hyper", "hypo")

  #samples and schema for means
  sampleSetList <- lapply(sampleSet,function(f){grep(f, colnames(betasitet), value=T)})
  names(sampleSetList) <- sampleSet
  tss <- paste0("(c(",paste(sampleSetList[["Tumour"]],collapse=","),"))")
  tsetmean <- paste0("mean",tss)
  pss <- paste0("(c(",paste(sampleSetList[["Plasma"]],collapse=","),"))")
  psetmean <- paste0("mean",pss)
  uss <- paste0("(c(",paste(sampleSetList[["Urine"]],collapse=","),"))")
  usetmean <- paste0("mean",uss)

  #per tumour, plasma, urine
  tprobens <- match(c("n2","n23","n24","n234"), names(setProbeVennList[[1]]))
  tprobeu <- unique(unlist(lapply(tprobens,function(f){
    hrs <- setProbeVennList[["hyper"]][[f]]
    hos <- setProbeVennList[["hypo"]][[f]]
    return(list(hrs,hos))
  })))
  tprobes <- tibble(probe=unlist(betasitet$probe)[unlist(betasitet$probe) %in% tprobeu])
  tbetasitet <- dplyr::bind_cols(tprobes, betasitet %>% dplyr::filter(probe %in% unlist(tprobes$probe)) %>%
                              rowwise() %>%
                              dplyr::mutate_(Tumor_Mean = tsetmean, Plasma_Mean = psetmean, Urine_Mean = usetmean) %>%
                              dplyr::select(Tumor_Mean, Plasma_Mean, Urine_Mean) %>%
                              round(.,3) %>%
                              ungroup())

  pprobens <- match(c("n3","n23","n34","n234"), names(setProbeVennList[[1]]))
  pprobeu <- unique(unlist(lapply(pprobens,function(f){
    hrs <- setProbeVennList[["hyper"]][[f]]
    hos <- setProbeVennList[["hypo"]][[f]]
    return(list(hrs,hos))
  })))
  pprobes <- tibble(probe=unlist(betasitet$probe)[unlist(betasitet$probe) %in% pprobeu])
  pbetasitet <- dplyr::bind_cols(pprobes, betasitet %>% dplyr::filter(probe %in% unlist(pprobes$probe)) %>%
                              rowwise() %>%
                              dplyr::mutate_(Tumor_Mean = tsetmean, Plasma_Mean = psetmean, Urine_Mean = usetmean) %>%
                              dplyr::select(Tumor_Mean, Plasma_Mean, Urine_Mean) %>%
                              round(.,3) %>%
                              ungroup())

  uprobens <- match(c("n4","n24","n34","n234"), names(setProbeVennList[[1]]))
  uprobeu <- unique(unlist(lapply(uprobens,function(f){
    hrs <- setProbeVennList[["hyper"]][[f]]
    hos <- setProbeVennList[["hypo"]][[f]]
    return(list(hrs,hos))
  })))
  uprobes <- tibble(probe=unlist(betasitet$probe)[unlist(betasitet$probe) %in% uprobeu])
  ubetasitet <- dplyr::bind_cols(uprobes, betasitet %>% dplyr::filter(probe %in% unlist(uprobes$probe)) %>%
                              rowwise() %>%
                              dplyr::mutate_(Tumor_Mean = tsetmean, Plasma_Mean = psetmean, Urine_Mean = usetmean) %>%
                              dplyr::select(Tumor_Mean, Plasma_Mean, Urine_Mean) %>%
                              round(.,3) %>%
                              ungroup())

  kdensPlot(meltdf=melt(tbetasitet[2:4]), pat="Tumor", nam="Mean Tumor 20%", colz=c("black","red","yellow"))
  kdensPlot(meltdf=melt(pbetasitet[2:4]), pat="Plasma", nam="Mean Plasma 20%", colz=c("black","red","yellow"))
  kdensPlot(meltdf=melt(ubetasitet[2:4]), pat="Urine", nam="Mean Urine 20%", colz=c("black","red","yellow"))

  ##########################
  ## EMRs FROM PROBE SETS ##
  ##########################

  ##EMRs from probes per patient, and split to hyper or hypo
  patEMRList <- eMRListFunc(inVec=patients, bed=epic.bed, betavalues=betasitet, filterVec=c(">0.8","<0.2"))
  names(patEMRList) <- patients

  ##per patient EMR overlap
  patEMROutList <- patEMRInList <- list()

  for(pg in 1:length(patGroups)){

    print(paste0("Working on: patient_", pg))
    pat.emr.out <- paste0(pat.emr.dir, "/patient_", pg)
    dir.create(pat.emr.out, showWarnings=F)

    p1 <- patEMRList[[grep(pg,patients)[1]]]
    p2 <- patEMRList[[grep(pg,patients)[2]]]
    p3 <- patEMRList[[grep(pg,patients)[3]]]
    p4 <- patEMRList[[grep(pg,patients)[4]]]

    patEMRInList[[pg]] <- list(p1,p2,p3,p4)
    names(patEMRInList[[pg]]) <- grep(pg, patients, value=TRUE)
    fullEMROut <- fullEMRFunc(EMRList=patEMRInList[[pg]], outdir=pat.emr.out, nam=paste0("patient_", pg))
    patEMROutList[[pg]] <- list("hyper", "hypo")
    patEMROutList[[pg]][[1]] <- fullEMROut[[1]]
    patEMROutList[[pg]][[2]] <- fullEMROut[[2]]
    names(patEMROutList[[pg]]) <- c("hyper", "hypo")
    names(patEMROutList[[pg]]$hyper) <- fullEMROut[[3]]
    names(patEMROutList[[pg]]$hypo) <- fullEMROut[[4]]
    for(x in 1:length(names(patEMROutList[[pg]]$hyper))){
      hyp <- c("hyper", "hypo")
      for(hy in 1:2){
        hypc <- hyp[hy]
        nmo <- paste0(pat.emr.out, "/", gsub(paste0("_", pg, "_", hypc), "", gsub(",", ".", names(patEMROutList[[pg]][[hypc]])[x])),".", hypc, ".patient_", pg, ".EMRs.GRanges.tsv")
        nmo <- gsub("FFPE_Benign", "MN", nmo)
        nmo <- gsub("FFPE_Tumour", "T", nmo)
        nmo <- gsub("Urine_Sediment", "U", nmo)
        nmo <- gsub("Plasma_cfDNA", "P", nmo)
        print(nmo)
        write.table(patEMROutList[[pg]][[hypc]][[names(patEMROutList[[pg]][[hypc]])[x]]], file=nmo, sep="\t", quote=F, row=F, col=T)
      }
    }
  }

  ##make from set of all EMRs in groups samples (i.e. not per sample, but in any set (not inc. 1))
  ##output dir
  setEMRList <- eMRListFunc(inVec=sampleSet, exVec="_1", bed=epic.bed, betavalues=betasitet, filterVec=c(">0.8","<0.2"))
  names(setEMRList) <- sampleSet

  ##make full setEMR outputs
  setEMROutList <- fullEMRFunc(EMRList=setEMRList, outdir=set.emr.dir, nam="set")
  names(setEMROutList) <- c("hyper", "hypo")
  names(setEMROutList$hyper) <- setEMROutList[[3]]
  names(setEMROutList$hypo) <- setEMROutList[[4]]
  hyp <- c("hyper", "hypo")
  for(x in 1:length(setEMROutList[[3]])){
  for(hy in 1:2){
    hypc <- hyp[hy]
    nmo <- paste0(set.emr.dir, "/", gsub(paste0("_", hypc), "", gsub(",", ".", names(setEMROutList[[hypc]])[x])), ".", hypc, ".EMRs.GRanges.tsv")
    nmo <- gsub("Benign", "MN", nmo)
    nmo <- gsub("Tumour", "T", nmo)
    nmo <- gsub("Urine", "U", nmo)
    nmo <- gsub("Plasma", "P", nmo)
    print(nmo)
    write.table(setEMROutList[[hypc]][x][[names(setEMROutList[[hypc]])[x]]], file=nmo, sep="\t", quote=F, row=F, col=T)
  }}

  #########################
  ##Save what we now have##
  #########################
  save(patProbeList,
       setProbeList,
       patGroupVennProbeList,
       setProbeVennList,
       patEMRList,
       setEMRList,
       setEMROutList,
       patEMROutList,
       sampleSet,
       patients,
       patGroups,
       betasitet,
       epic.bed,
       file=paste0(OUTDIR, "/set_patient-probe_EMRs.GRanges.RData"))
}

#! /usr/local/bin/R
libs <- c("GGally", "magrittr", "tidyverse", "GenomicRanges", "reshape2", "ggplot2", "VennDiagram")
libsLoaded <- lapply(libs,function(l){;suppressWarnings(suppressMessages(library(l, character.only = TRUE)))})

eMRListFunc <- function(betavalues, bed, inVec=NULL, exVec=NULL, filterVec=NULL){
  ##betavalues is a tibble of beta values per sample, and a 'probe' column
  ##bed is a bedfile from EPIC annotation, in bed format
  ##inVec is set of identifiers that group the samples in betavalues by colnames
  ##exVec is grep for exclusion of samples
  ##filterVec is vector of conditions e.g. c(">0.8","<0.2")
  if(length(grep("probe", colnames(betavalues)))==0){
    print("Betavalues requires probe column")
    break
  }
  if(is.null(inVec)){
    print("Specify input vector of sample/group IDs grep'd from colnames of betavalues")
    break
  }
  if(is.null(filterVec)){
    print("Specify filters required for sample/group IDs of betavalues")
    break
  }
  if(!is.null(exVec)){
    exCols <- unlist(lapply(exVec, function(xx){
      grep(exVec,colnames(betavalues))
    }))
    betavalues <- betavalues[,-exCols]
  }
  outEMRList <- lapply(seq_along(inVec), function(x){
    sampl <- inVec[x]
    print(paste0("Working on: ", sampl))
    colnamesT <- grep(sampl,colnames(betavalues))
    colnamesSet <- colnames(betavalues)[colnamesT]

    listo <- lapply(filterVec, function(flt){
      filterQuery <- paste(unlist(lapply(colnamesSet,function(f){
        paste(f,flt)})), collapse="&")
      outEMR <- betavalues %>% dplyr::select(colnamesSet, probe) %>% dplyr::filter_(filterQuery)
      bed <- bed %>% dplyr::filter(probe %in% as.vector(outEMR$probe))
      gR <- GRanges(seqnames=bed$seqnames,
                        ranges=IRanges(start=bed$start,
                                       end=bed$end),
                        strand=bed$strand)
      ##this defines 2kb between probes which we specify as the EMR
      gRR <- GenomicRanges::reduce(gR, min.gapwidth=2000L, with.revmap=TRUE)
      eMR <- as_tibble(gRR$revmap) %>% group_by(group) %>% filter(length(group)>2)
      groupLengths <- as.vector(table(eMR$group))
      eMR <- gRR[unique(as.vector(eMR$group))]
      mcolseMRList <- as.list(1:length(eMR$revmap))
      bedanno <- as.vector(bed$anno)
      print("Annotating...")
      pb <- txtProgressBar(min = 0, max = length(eMR$revmap), initial = 0, char = "=",
                     width = NA, style = 3)
      for(y in 1:length(eMR$revmap)){
        setTxtProgressBar(pb, y)
        strsp <- strsplit(as.vector(bedanno[eMR$revmap[[y]]]),";")
        mcolseMRList[[y]] <- apply(as.data.frame(lapply(strsp,function(fs){
          if(length(fs)!=5){return(c(fs,""))}else{return(fs)}})),
                                   1,
                                  function(ff){paste(ff,collapse=",")})
      }
      close(pb)
      values(eMR) <- cbind(values(eMR), DataFrame(groupLengths),
                                        DataFrame(sapply(mcolseMRList, function(x){x[1]})),
                                        DataFrame(sapply(mcolseMRList, function(x){x[2]})),
                                        DataFrame(sapply(mcolseMRList, function(x){x[3]})),
                                        DataFrame(sapply(mcolseMRList, function(x){x[4]})),
                                        DataFrame(sapply(mcolseMRList, function(x){x[5]}))
                                      )
      names(values(eMR)) <- c("revmap", "groupLengths", "probe", "UCSC_RefGene_Name", "UCSC_RefGene_Accession", "UCSC_RefGene_Group", "Relation_to_UCSC_CpG_Island")
      return(eMR)
    })
    names(listo) <- c("hyper","hypo")
    return(listo)
  })
  return(outEMRList)
}

##specify all Tumour EMRs (patients 2:4) not in MN
allEMR <- function(queryTerm,exterm,setVec,listData){
  unique(do.call(c,lapply(grep(queryTerm,setVec[grep(exterm,setVec,invert=T)]),function(f){
    listData[[f]]})))
}

##full set of EMRs, and Venns of contribution of each
fullEMRFunc <- function(EMRList, outdir, nam){

  if(class(EMRList[[1]][[1]])[[1]] != "GRanges"){
    print("Requires EMRList[[x]][[y]] as GRanges object")
    break
  }
  fullEMRList <- list(EMRList[[1]][[1]][1,], EMRList[[1]][[1]][1,])
  names(fullEMRList) <- c("hyper", "hypo")

  hyp <- c("hyper", "hypo")
  ooList <- list(1,2,3)

  for (h in 1:2){
    hypc <- hyp[h]
    print(paste0("Working on: ", hypc))

    fullEMRList[[hypc]]$includes <- "FAKE"
    for(x in 1:length(names(EMRList))){
      emr <- EMRList[[x]][[hypc]]
      emr$includes <- paste0(names(EMRList)[[x]],"_",hypc)
      fullEMRList[[hypc]] <- c(fullEMRList[[hypc]], emr)
    }

    ##merge GRanges
    fullEMR <- fullEMRList[[hypc]][-1,-1]
    rfullEMR <- GenomicRanges::reduce(fullEMR, with.revmap=TRUE)

    ##specify which samples are found in which merged regions
    rfullEMR$includes <- unlist(lapply(rfullEMR$revmap, function(f){
      paste(unique(sort(fullEMR[f]$includes)), collapse=",")
    }))

    ##tabulate number of regions found to include which samples and venn4
    rfullEMRvenn4 <- table(sort(rfullEMR$includes))
    listNms <- names(EMRList)

    nms2 <- apply(combn(x=4,m=2), 2, function(f){paste(sort(paste0(listNms[f],"_",hypc)), collapse=",")})
    nms3 <- apply(combn(x=4,m=3), 2, function(f){paste(sort(paste0(listNms[f],"_",hypc)), collapse=",")})
    nms4 <- apply(combn(x=4,m=4), 2, function(f){paste(sort(paste0(listNms[f],"_",hypc)), collapse=",")})

    area01 <- as.numeric(sum(rfullEMRvenn4[grep(listNms[1], names(rfullEMRvenn4))]))
    area02 <- as.numeric(sum(rfullEMRvenn4[grep(listNms[2], names(rfullEMRvenn4))]))
    area03 <- as.numeric(sum(rfullEMRvenn4[grep(listNms[3], names(rfullEMRvenn4))]))
    area04 <- as.numeric(sum(rfullEMRvenn4[grep(listNms[4], names(rfullEMRvenn4))]))
    areas <- c(area01, area02, area03, area04)

    n012 <- as.numeric(rfullEMRvenn4[names(rfullEMRvenn4) %in% nms2[1]][[1]])
    n013 <- as.numeric(rfullEMRvenn4[names(rfullEMRvenn4) %in% nms2[2]][[1]])
    n014 <- as.numeric(rfullEMRvenn4[names(rfullEMRvenn4) %in% nms2[3]][[1]])
    n023 <- as.numeric(rfullEMRvenn4[names(rfullEMRvenn4) %in% nms2[4]][[1]])
    n024 <- as.numeric(rfullEMRvenn4[names(rfullEMRvenn4) %in% nms2[5]][[1]])
    n034 <- as.numeric(rfullEMRvenn4[names(rfullEMRvenn4) %in% nms2[6]][[1]])
    n0123 <- as.numeric(rfullEMRvenn4[names(rfullEMRvenn4) %in% nms3[1]][[1]])
    n0124 <- as.numeric(rfullEMRvenn4[names(rfullEMRvenn4) %in% nms3[2]][[1]])
    n0134 <- as.numeric(rfullEMRvenn4[names(rfullEMRvenn4) %in% nms3[3]][[1]])
    n0234 <- as.numeric(rfullEMRvenn4[names(rfullEMRvenn4) %in% nms3[4]][[1]])
    n01234 <- as.numeric(rfullEMRvenn4[names(rfullEMRvenn4) %in% nms4[1]][[1]])

    n01 <- area01 - (n012 + n013 + n014 + n0123 + n0134 + n0124 + n01234)
    n02 <- area02 - (n012 + n023 + n024 + n0123 + n0124 + n0234 + n01234)
    n03 <- area03 - (n013 + n023 + n034 + n0123 + n0134 + n0234 + n01234)
    n04 <- area04 - (n014 + n024 + n034 + n0124 + n0134 + n0234 + n01234)
    areaVec <- c(n03, n034, n04, n013, n0134, n01234, n0234, n024, n01, n014, n0124, n0123, n023, n02, n012)

    listNms[grep("Tumour", listNms)] <- gsub("Tumour", "Tumor", listNms[grep("Tumour", listNms)])
    category <- paste0(c("Matched Normal", "Tumor", "Plasma", "Urine"),"\n", hypc, " EMRs","\n","n = ", prettyNum(areas, big.mark=",", trim=TRUE))
    pdf(paste0(outdir,"/", nam, ".", hypc, "_EMRs.venn4.pdf"), width=14,height=10)
      draw.quad.venn(area.vector=prettyNum(areaVec, big.mark=",", trim=TRUE), direct.area=TRUE, category=category, cat.default.pos=c("outer"), cex=1.8, cat.cex=1.8, fill=c("forestgreen","black","darkred","yellow"), cat.dist=c(0.25, 0.25, 0.135, 0.135))
    dev.off()

    ##write each individual set out
    inc <- rfullEMR$includes
    uinc <- unique(inc)
    ooList[[2+h]] <- uinc
    oList <- as.list(uinc)

    for (th in 1:length(uinc)){
      print(paste0("Working on: ", uinc[th]))

      uincs <- strsplit(uinc[th],",")[[1]]
      lenguincs <- length(uincs)

      ##if only single sample included
      if(lenguincs==1){
        oList[[th]] <- fullEMR[unlist(lapply(rfullEMR[inc==uinc[th]]$revmap, function(f){return(f)}))]
      }

      else{

        ##combine each two/three/four revmap ID lines from fullEMR into single line
        catGR <- do.call(c,lapply(rfullEMR[inc==uinc[th]]$revmap, function(f){

          # print(paste0("TH -> ", th, ";F -> "))
          # print(f)
          ##rows to work on
          ff <- fullEMR[unlist(f)]
          ##possible that all rows are exactly same
          if(length(unique(ff))==1){
            rowOut <- unique(ff)
            samp <- uincs[uincs %in% rowOut$includes]
            for(ii in grep(samp, uincs, invert=T)){

              rowOut$includes <- paste0(rowOut$includes, ";", ff[ii]$includes, "[all]")
            }
          }

          ##max groupLengths
          else{
            #ff <- unique(ff)
            sampleng <- 1:length(ff$includes)
            rowUse <- which(ff$groupLengths %in% max(ff$groupLengths), arr.ind=TRUE)[1]
            rowOut <- ff[rowUse]
            ##attemtp to collapse rest
            for(ii in grep(rowUse, sampleng, invert=T)){
              matchpr <- match(strsplit(ff[ii]$probe,",")[[1]], strsplit(rowOut$probe,",")[[1]])
              rowOut$includes <- paste0(rowOut$includes, ";", ff[ii]$includes, "[", paste(matchpr,collapse=","), "]")
            }
          }
          return(rowOut)
        }))
        oList[[th]] <- catGR
      }
    }
    ooList[[h]] <- oList
  }
  return(ooList)
}


##function to compare overlap of GRanges, output list of uniq to [1] v. [2]
uniqEMR <- function(gr1,gr2,namegrs){
  ##take overlaps between sets
  over12 <- findOverlaps(gr1, gr2)
  ##all overlap in 1, then 2
  qh1 <- queryHits(over12)
  qh2 <- subjectHits(over12)
  ##unique regions to 1, 2 when removing 'hit' regions
  ##therefore (eg 1):not 1:x, when 1:x found in 2 (queryHits)
  oList <- list(unique(gr1[c(1:length(gr1))[! 1:length(gr1) %in% unique(sort(qh1))]]),
          unique(gr2[c(1:length(gr2))[! 1:length(gr2) %in% unique(sort(qh2))]]))
  names(oList) <- c(paste0(namegrs[1],".",namegrs[2]),
                    paste0(namegrs[2],".",namegrs[1]))
  return(oList)
}
##function to compare overlap of GRanges, output list of uniq to [1] v. [2]
uniq3EMR <- function(gr1,gr2,gr3,namegrs){
  ##take overlaps between sets
  over12 <- findOverlaps(gr1, gr2)
  ##all overlap in 1, then 2
  qh1 <- queryHits(over12)
  qh2 <- subjectHits(over12)
  ##unique regions to 1, 2 when removing 'hit' regions
  ##therefore (eg 1):not 1:x, when 1:x found in 2 (queryHits)
  gr1o <- unique(gr1[c(1:length(gr1))[! 1:length(gr1) %in% unique(sort(qh1))]])
  gr2o <- unique(gr2[c(1:length(gr2))[! 1:length(gr2) %in% unique(sort(qh2))]])
  ##queryhits over 3
  qh1o3 <- queryHits(findOverlaps(gr1o, gr3))
  qh2o3 <- queryHits(findOverlaps(gr2o, gr3))
  oList <- list(unique(gr1o[c(1:length(gr1o))[! 1:length(gr1o) %in% unique(sort(qh1o3))]]),
                unique(gr2o[c(1:length(gr2o))[! 1:length(gr2o) %in% unique(sort(qh2o3))]]))
  names(oList) <- c(paste0(namegrs[1],".",namegrs[2],"-",namegrs[3]),
                    paste0(namegrs[2],".",namegrs[1],"-",namegrs[3]))
  return(oList)
}
##inverse of above, doing the same not the unique
sameEMR <- function(gr1,gr2,namegrs){
  ##take overlaps between sets
  over12 <- findOverlaps(gr1, gr2)
  ##all overlap in 1, then 2
  qh1 <- queryHits(over12)
  qh2 <- subjectHits(over12)
  ##unique regions to 1, 2 when removing 'hit' regions
  ##therefore (eg 1): 1:x, when 1:x found in 2 (queryHits)
  oList <- list(unique(gr1[c(1:length(gr1))[1:length(gr1) %in% unique(sort(qh1))]]),
          unique(gr2[c(1:length(gr2))[1:length(gr2) %in% unique(sort(qh2))]]))
  names(oList) <- c(paste0(namegrs[1],"_",namegrs[2]),
                    paste0(namegrs[2],"_",namegrs[1]))
  return(oList)
}

##probe sets as list per patient [[patientx]][[hyper, hypo]]
probeHyperHypoLister <- function(betavalues, samples, exclude=NULL, filters){
  options(pillar.sigfig=7)
  if(!is.null(exclude)){
    colsOut <- grep(exclude, colnames(betavalues))
    betavalues <- betavalues[,-colsOut]
  }
  lapply(seq_along(samples), function(f){
    print(paste0("Working on: ", samples[f]))
    colnamesT <- grep(samples[f],colnames(betavalues))
    colnamesSet <- colnames(betavalues)[colnamesT]
    listo <- lapply(filters, function(flt){
      filterQuery <- paste(unlist(lapply(colnamesSet,function(ff){
        paste(ff,flt)})), collapse="&")
      betavalues %>% dplyr::select(colnamesSet, probe) %>%
                     dplyr::filter_(filterQuery)
    })
    names(listo) <- c("hyper","hypo")
    return(listo)
  })
}

kdensPlot <- function(meltdf, pat, nam, colz){

  colnames(meltdf) <- c("Sample", "Beta")
  colevs <- levels(meltdf$Sample)

  ##colour in meltdf
  colVec <- c()
  for(ff in 1:length(colevs)){
      colVec <- c(colVec, rep(colz[ff],length=table(meltdf$Sample==colevs[ff])[[2]]))
  }
  meltdf$Colour <- colVec

  ggplot(meltdf, aes(x=Beta, group=Sample, fill=Colour)) +
  geom_density(alpha=0.5, show.legend=TRUE) +
  scale_fill_manual(values=colz, labels=colevs) +
  labs(x="Beta Value",
       y="Density",
       title=paste0("Beta Value Distribution, ", nam, " probes (n = ", prettyNum(dim(meltdf)[1]/3, trim=T, big.mark=","), ")"))
  ggsave(paste0(pat, ".probes.kdensity.", strsplit(nam," ")[[1]][1], ".pdf"))

}

venn4ProbeList <- function(A, B, C, D, sampVec, tag, nam, listo, fil=NULL){
  if(is.null(fil)){
    print("No fill colours provided, defaults will be used!")
    fil <- c("red", "green", "dodgerblue", "yellow")
  }
  sampleVec <- unlist(lapply(sampVec, function(f){gsub(" ","-",f)}))
  sampVec <- paste0(sampVec,paste0("\n", tag))

  area1 <- length(A)
  area2 <- length(B)
  area3 <- length(C)
  area4 <- length(D)

  Anam <- sampVec[1]
  Bnam <- sampVec[2]
  Cnam <- sampVec[3]
  Dnam <- sampVec[4]

  #found in all
  n1234 <- A[A %in% B & A %in% C & A %in% D]

  #threes
  n123 <- A[A %in% B & A %in% C &! A %in% n1234]
  n124 <- A[A %in% B & A %in% D &! A %in% n1234]
  n134 <- A[A %in% C & A %in% D &! A %in% n1234]
  n234 <- B[B %in% C & B %in% D &! B %in% n1234]

  #twos
  n12 <- A[A %in% B &! A %in% n123 &! A %in% n124 &! A %in% n1234]
  n13 <- A[A %in% C &! A %in% n123 &! A %in% n134 &! A %in% n1234]
  n14 <- A[A %in% D &! A %in% n124 &! A %in% n134 &! A %in% n1234]
  n23 <- B[B %in% C &! B %in% n123 &! B %in% n234 &! B %in% n1234]
  n24 <- B[B %in% D &! B %in% n124 &! B %in% n234 &! B %in% n1234]
  n34 <- C[C %in% D &! C %in% n134 &! C %in% n234 &! C %in% n1234]

  #singles
  n1 <- A[! A %in% n12 &! A %in% n13 &! A %in% n14 &! A %in% n123 &! A %in% n134 &! A %in% n124 &! A %in% n1234]
  n2 <- B[! B %in% n12 &! B %in% n23 &! B %in% n24 &! B %in% n123 &! B %in% n124 &! B %in% n234 &! B %in% n1234]
  n3 <- C[! C %in% n13 &! C %in% n23 &! C %in% n34 &! C %in% n123 &! C %in% n134 &! C %in% n234 &! C %in% n1234]
  n4 <- D[! D %in% n14 &! D %in% n24 &! D %in% n34 &! D %in% n124 &! D %in% n134 &! D %in% n234 &! D %in% n1234]

  #area.vector
  areaVec <- c(length(n3), length(n34), length(n4), length(n13), length(n134), length(n1234), length(n234), length(n24), length(n1), length(n14), length(n124), length(n123), length(n23), length(n2), length(n12))
  pareaVec <- prettyNum(areaVec, trim=TRUE, big.mark=",")

  category <- c(paste0(Anam, "\n", prettyNum(area1, trim=TRUE, big.mark=",")),
                paste0(Bnam, "\n", prettyNum(area2, trim=TRUE, big.mark=",")),
                paste0(Cnam, "\n", prettyNum(area3, trim=TRUE, big.mark=",")),
                paste0(Dnam, "\n", prettyNum(area4, trim=TRUE, big.mark=",")))

  pdf(paste0(nam,".venn4.pdf"),width=20,height=14)
  draw.quad.venn(area.vector=pareaVec, direct.area=TRUE,
    category=category,
    cat.default.pos=c("outer"),
    cex=1.8,cat.cex=1.8,
    fill=fil)
  dev.off()

  listo <- list(n1,n2,n3,n4,n12,n13,n14,n23,n24,n34,n123,n124,n134,n234,n1234)
  nlisto <- c("n1","n2","n3","n4","n12","n13","n14","n23","n24","n34","n123","n124","n134","n234","n1234")
  unlist(lapply(seq_along(nlisto), function(f){
    nlist <- paste(sampleVec[as.numeric(grep("n",strsplit(nlisto[f], "")[[1]], invert=T, value=T))], collapse="_")
    write.table(listo[[f]], file=paste0(nam, ".", nlist , ".probes.txt"), row=F, col=F, quote=F)
  }))
  names(listo) <- nlisto
  return(listo)
}

heatmapCorPair <- function(hyper=NULL, hypo=NULL, sampleset=NULL, samplerenames=NULL){

    if(is.null(sampleset)){print("Sample set identifier (found in colnames) required"); break}
    if(is.null(hyper) & is.null(hypo)){
      print("Input of at least hyper or hypo data required (matrix, data.frame, tibble of beta values)")
      break
    }
    else{
      ##hyper above diagonal, hypo below
      used <- list(hypo, hyper)
      names(used) <- c("hypo", "hyper")
      if(is.null(hyper)){used <- list(1,hypo); names(used) <- c(1,"hypo")}
      if(is.null(hypo)){used <- list(hyper,1); names(used) <- c("hyper",1)}
      ##means per 'used' entries per sampleset IDs
      meansetList <- lapply(seq_along(used), function(f){
        fsets <- used[[f]]
        meanfset <- lapply(sampleset,function(ff){
          rowMeans(fsets[,grep(ff,colnames(fsets))])
        })
      })
      corpso <- lapply(seq_along(meansetList), function(f){
        return(apply(combn(x=length(sampleset), m=2), 2, function(ff){
          corp <- cor.test(meansetList[[f]][[ff[1]]],
                           meansetList[[f]][[ff[2]]],
                           method="spearman")
          return(corp$estimate[[1]])
        }))
      })
      corpsi <- matrix(ncol=length(sampleset), nrow=length(sampleset))
      combni <- combn(x=length(sampleset), m=2)
      for(f in 1:2){
        if(f==1){
          for(x in 1:length(combni[1,])){
            corpsi[combni[1,x],combni[2,x]] <- corpso[[f]][x]
          }
        }
        if(f==2){
          for(x in 1:length(combni[1,])){
            corpsi[combni[2,x],combni[1,x]] <- -1*corpso[[f]][x]
          }
        }
      }
      if(is.null(samplerenames)){
        colnames(corpsi) <- sampleset -> rownames(corpsi)
      }
      if(!is.null(samplerenames)){
        colnames(corpsi) <- samplerenames -> rownames(corpsi)
      }
    return(corpsi)
  }
}

ggplotcor <- function(plotdf, names){
  smgam <- summary(gam(plotdf$x ~ plotdf$y))
  ggpo <- ggplot(plotdf, aes(x=x, y=y)) +
    geom_hex(bins=100) +
    scale_fill_gradient(low="grey90", high="black") +
    geom_smooth(method="auto",
                formula="y ~ x",
                linetype=2, size=0.5, colour="black",
                se=T) +
    annotate("text",
             label = paste0("R-sq. = ",round(smgam$r.sq,digits=3)),
             x = 0.4, y = 0.625, size = 4, colour = "black") +
    ggtitle(label = paste0(names[1],"-vs-", names[2],"_3x3 (n_probes = ", dim(plotdf)[1], ")")) +
    labs(x=names[1], y=names[2]) +
    theme(plot.title = element_text(hjust = 0.5),
          panel.background = element_rect(fill = "white"),
          panel.grid.major = element_line(colour = "black"),
          panel.grid.minor = element_line(colour = "black"))
  return(ggpo)
}

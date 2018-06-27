#! /usr/local/bin/R

libs <- c("GGally", "tidyverse", "magrittr", "pheatmap", "RnBeads", "mgcv", "ggplot2", "readxl")
print("Loading libraries...")
libsLoaded <- lapply(libs,function(l){suppressWarnings(suppressMessages(library(l, character.only = TRUE)))})

##script takes beta values across samples, set of probes (here in each of liquid)
##takes probeset across samples, changes those >0.2, < 0.8 to NA
##produces correlation matrix to visulaise similarity per sample types

#input arguments
argsin <- commandArgs(trailingOnly=TRUE)
print(argsin)

##inputs
BASEDIR <- argsin[1]
SCRIPTDIR <- argsin[2]
TAG <- argsin[3]

#source functions
source(paste0(SCRIPTDIR, "/set_patient_EMR.func.R"))

#set dirs, load data
OUTDIR <- paste0(BASEDIR, "/analysis/probe_EMR")
load(paste0(OUTDIR,"/set_patient-probe_EMRs.GRanges.RData"))
probe.dir <- paste0(OUTDIR, "/probes/set")
roadmap.dir <- paste0(BASEDIR,"/Roadmap_450K")

##filter patient_1
not1 <- grep("_1", colnames(betasitet), invert=T, value=T)
filterQuery <- paste(paste0(grep("probe", not1, invert=T, value=T),
                      " > 0.8 | ",
                      not1,
                      " < 0.2"), collapse=" | ")
beta3set <- betasitet %>% dplyr::select(not1) %>%
                          dplyr::filter_(filterQuery)

#scan vectors of probes
urprobes <- scan(paste0(probe.dir,"/Urine.hyper.probes.txt"), what="character")
uoprobes <- scan(paste0(probe.dir,"/Urine.hypo.probes.txt"), what="character")
prprobes <- scan(paste0(probe.dir,"/Plasma.hyper.probes.txt"), what="character")
poprobes <- scan(paste0(probe.dir,"/Plasma.hypo.probes.txt"), what="character")

uprobes <- c(urprobes, uoprobes)
pprobes <- c(prprobes, poprobes)

#bvalues
ubeta3set <- beta3set %>% dplyr::filter(probe %in% uprobes)
pbeta3set <- beta3set %>% dplyr::filter(probe %in% pprobes)

#means
umean3set <- as_tibble(data.frame(MN_mean=rowMeans(ubeta3set[,grep("Benign", colnames(ubeta3set))]),
                                  Tumor_mean=rowMeans(ubeta3set[,grep("Tumour", colnames(ubeta3set))]),
                                  Plasma_mean=rowMeans(ubeta3set[,grep("Plasma", colnames(ubeta3set))]),
                                  Urine_mean=rowMeans(ubeta3set[,grep("Urine", colnames(ubeta3set))]))) %>%
             round(.,3)

pmean3set <- as_tibble(data.frame(MN_mean=rowMeans(pbeta3set[,grep("Benign", colnames(ubeta3set))]),
                                  Tumor_mean=rowMeans(pbeta3set[,grep("Tumour", colnames(ubeta3set))]),
                                  Plasma_mean=rowMeans(pbeta3set[,grep("Plasma", colnames(ubeta3set))]),
                                  Urine_mean=rowMeans(pbeta3set[,grep("Urine", colnames(ubeta3set))]))) %>%
             round(.,3)

##set all above 0.2, below 0.8 as NA?
##don't! nicer plots if we show intermediate
# ubeta3set %<>% transmute_at(.vars = colnames(ubeta3set), .funs = funs(if_else(.>0.8 | .<0.2,.,-1)))
# ubeta3set[as.data.frame(ubeta3set) == -1] <- NA
#
# pbeta3set %<>% transmute_at(.vars = colnames(pbeta3set), .funs = funs(if_else(.>0.8 | .<0.2,.,-1)))
# pbeta3set[as.data.frame(pbeta3set) == -1] <- NA

##correlation matrix; no longer wanted
# ggu <- ggpairs(umean3set, upper = list(continuous = wrap("cor", size = 6)))
# ggp <- ggpairs(pmean3set, upper = list(continuous = wrap("cor", size = 6)))
#
# ggu <- ggu + theme(axis.text.x = element_text(size = 5), axis.text.y = element_text(size = 5))
# ggp <- ggp + theme(axis.text.x = element_text(size = 5), axis.text.y = element_text(size = 5))
#
# ggsave(ggu,file=paste0(probedir,"/Urine.hyper-hypo.methylation_per_sample.annotate-sites.means.corrMat.png"), dpi=1200)
# ggsave(ggp,file=paste0(probedir,"/Plasma.hyper-hypo.methylation_per_sample.annotate-sites.means.corrMat.png"), dpi=1200)
#
####again, with feeling
##heatmap hyper, hypo over, under diagonal for urine, plasma
urbeta <- beta3set %>% dplyr::filter(probe %in% urprobes)
uobeta <- beta3set %>% dplyr::filter(probe %in% uoprobes)
prbeta <- beta3set %>% dplyr::filter(probe %in% prprobes)
pobeta <- beta3set %>% dplyr::filter(probe %in% poprobes)

#lapply to get mean per sample set, then cor each set against others, for each hyper, o, and then plot
sampleSet <- c("Benign", "Tumour", "Plasma", "Urine")
print("Plotting correlation heatmaps...")
ucor <- heatmapCorPair(hyper=urbeta,hypo=uobeta,sampleset=sampleSet,samplerenames=c("MN","Tumor","Plasma","Urine"))
pcor <- heatmapCorPair(hyper=prbeta,hypo=pobeta,sampleset=sampleSet,samplerenames=c("MN","Tumor","Plasma","Urine"))

pdf(paste0(probe.dir, "/Urine_hyper-hypo_probes.correlation.heatmap.pdf"), onefile=F)
  pheatmap(ucor, color=colorRampPalette(c("darkblue","blue","white","white","white","white","white","white","yellow","darkgoldenrod"))(100), display_numbers=TRUE, number_color="white", fontsize_number=12, cluster_col = FALSE, cluster_row = FALSE)
  dev.off()
pdf(paste0(probe.dir, "/Plasma_hyper-hypo_probes.correlation.heatmap.pdf"), onefile=F)
  pheatmap(pcor, color=colorRampPalette(c("darkblue","blue","white","white","white","white","white","white","yellow","darkgoldenrod"))(100), display_numbers=TRUE, number_color="white", fontsize_number=12, cluster_col = FALSE, cluster_row = FALSE)
dev.off()


##############
## 450K GAM ##
##############
#load 450K data
betatab <- readxl::read_excel(paste0(roadmap.dir, "/Roadmap_BMIQ_Beta_Values_Complete.xlsx"))
benID <- c("ID",grep("B",colnames(betatab), value=TRUE))
benbeta <- betatab %>% dplyr::select(benID) %>%
                       dplyr::rename(probe = ID) %>%
                       dplyr::arrange(probe)

##filter < 0.2, > 0.8
benIDn <- benID[2:length(benID)]
filterQuery <- paste(benIDn, " > 0.8 | ", benIDn, " < 0.2", collapse=" | ")
benbetahyp <- benbeta %>% dplyr::filter_(filterQuery)

#write output
write_tsv(as.data.frame(benbetahyp), paste0(roadmap.dir, "/Roadmap_BMIQ_Beta_Values_Benign_EMP.txt"))

#compare to EPIC
K450inEPIC <- sort(unlist(benbeta$probe)[unlist(benbeta$probe) %in% unlist(beta3set$probe)])
EPICinK450 <- sort(unlist(beta3set$probe)[unlist(beta3set$probe) %in% unlist(benbeta$probe)])

##Matched Normal
MNss <- paste0("(c(",paste(grep("Benign",colnames(beta3set), value=T) ,collapse=","),"))")
MNsetmean <- paste0("mean",MNss)
MNrmeans <- beta3set %>% dplyr::filter(probe %in% EPICinK450) %>%
                        dplyr::select(grep("Benign",colnames(beta3set))) %>%
                        rowwise() %>%
                        dplyr::mutate_(MN_mean = MNsetmean) %>%
                        dplyr::select(MN_mean)

##Tumour
Tuss <- paste0("(c(",paste(grep("Tumour",colnames(beta3set), value=T) ,collapse=","),"))")
Tusetmean <- paste0("mean",Tuss)
Turmeans <- beta3set %>% dplyr::filter(probe %in% EPICinK450) %>%
                        dplyr::select(grep("Tumour",colnames(beta3set))) %>%
                        rowwise() %>%
                        dplyr::mutate_(Tu_mean = Tusetmean) %>%
                        dplyr::select(Tu_mean)
##True Benign
TBss <- paste0("(c(",paste(grep("B",colnames(benbeta), value=T) ,collapse=","),"))")
TBsetmean <- paste0("mean",TBss)
TBrmeans <- benbeta %>% dplyr::filter(probe %in% EPICinK450) %>%
                        dplyr::select(grep("B",colnames(benbeta))) %>%
                        rowwise() %>%
                        dplyr::mutate_(TB_mean = TBsetmean) %>%
                        dplyr::select(TB_mean)

##named x, y for ggplot; x->MN/Tumour, y->TN
plotMNTB <- data.frame(x=c(MNrmeans),
                     y=c(TBrmeans))
plotMNTu <- data.frame(x=c(MNrmeans),
                       y=c(Turmeans))
colnames(plotMNTB) <- colnames(plotMNTu) <- c("x", "y")

##no longer using
# ggplotMNTB <- ggplotcor(plotMNTB, c("MN","TB"))
# ggplotMNTu <- ggplotcor(plotMNTu, c("MN","Tumor"))

##plot all into single plot
smgamTB <- summary(gam(plotMNTB$x ~ plotMNTB$y))
smgamTu <- summary(gam(plotMNTu$x ~ plotMNTu$y))

##melt
mdfTB <- melt(data.frame(probe=EPICinK450,
                  x=plotMNTB$x,
                  y=plotMNTB$y),
            id.vars = c("x", "y"))
mdfTB$variable <- "MN vs. TB"
mdfTu <- melt(data.frame(probe=EPICinK450,
                  x=plotMNTu$x,
                  y=plotMNTu$y),
            id.vars = c("x", "y"))
mdfTu$variable <- "MN vs. T"
mdf <- rbind(mdfTB, mdfTu)

print("Plotting EPIC-450K GAM plots...")
ggpo <- ggplot() +
  geom_hex(data=mdf, aes(x=x, y=y, fill=variable), alpha=0.5, bins=100) +
  scale_fill_manual(values = c("black", "dodgerblue")) +
  geom_smooth(data=mdfTB, aes(x=x, y=y),
              method="auto",
              formula="y ~ x",
              linetype=4, size=1, colour="white",
              se=T) +
  geom_smooth(data=mdfTu, aes(x=x, y=y),
              method="auto",
              formula="y ~ x",
              linetype=5, size=1, colour="white",
              se=T) +
  annotate("text",
           label = paste0("TB vs. MN R-sq. = ", round(smgamTB$r.sq,digits=3)),
           x = 0.63, y = 0.425, size = 4, colour = "white") +
  annotate("text",
           label = paste0("T vs. MN R-sq. = ", round(smgamTu$r.sq,digits=3)),
           x = 0.35, y = 0.65, size = 4, colour = "white") +
  ggtitle(label = "Matched Normal (MN) vs. True Benign (TB), Tumor (T)",
          subtitle = paste0("Probes = ", prettyNum(dim(mdfTB)[1], big.mark=","))) +
  labs(x=expression(paste("MN ", beta,"-value")),
       y=expression(paste("TB, T ", beta,"-value")))

ggsave(filename=paste0(probe.dir, "/MN-TN_MN-Tu_3x4.pdf"), plot=ggpo, dpi=1200)

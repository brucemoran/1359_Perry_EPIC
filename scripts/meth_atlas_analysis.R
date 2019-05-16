#! R

##clustering analysis using meth_atlas reference data
library(tidyverse)
library(ggdendro)

#inputs
argsin <- commandArgs(trailingOnly=TRUE)

##data is held in [1]; same dir for output, but each gets a dir
OUTDIR <- paste0(argsin[1], "/analysis/meth_atlas")
dir.create(OUTDIR, showWarnings=FALSE)
BETASITES <- argsin[2]
REFATLAS <- argsin[3]

betasitet <- read_tsv(BETASITES)

atlas <- read_csv(REFATLAS) %>% rename(probe="CpGs") %>% arrange(probe) %>% unique()

##join, hclust
atlas.betasitet <- inner_join(atlas,betasitet) %>%
                   na.omit() %>%
                   column_to_rownames("probe") %>%
                   as.data.frame()
d <- dist(t(atlas.betasitet), method = "euclidean")
hc <- hclust(d, method="average")
dend <- as.dendrogram(hc)
ddata <- dendro_data(dend)

##colour by Group from our data
ddata$labels$Group <- "atlas"
ddata$labels$Group[grep("FFPE_Benign", ddata$labels$label)] <- "Matched Normal"
ddata$labels$Group[grep("FFPE_Tumour", ddata$labels$label)] <- "Tumor"
ddata$labels$Group[grep("Plasma_cfDNA", ddata$labels$label)] <- "Plasma"
ddata$labels$Group[grep("Urine_Sediment", ddata$labels$label)] <- "Urine Sediment"

ddata$labels$label <- gsub("FFPE_Benign_", "Matched Normal ", ddata$labels$label)
ddata$labels$label <- gsub("FFPE_Tumour_", "Tumor ", ddata$labels$label)
ddata$labels$label <- gsub("Plasma_cfDNA_", "Plasma ", ddata$labels$label)
ddata$labels$label <- gsub("Urine_Sediment_", "Urine Sediment ", ddata$labels$label)

##ggdendroogram
# ggdendrogram(hc, rotate=TRUE, labels=FALSE) +
#   geom_text(data=ddata$labels,
#             aes(label=label, x=x, y=y-0.5, colour=factor(group)),
#             hjust = 0) +
#   scale_colour_manual(values=c("grey", "forestgreen", "black", "darkred", "yellow"))  +
#   scale_y_reverse(expand=c(0, 20))

#https://stackoverflow.com/questions/12630565/adding-labels-to-a-dendrogram-in-ggplot-using-ggdendro-in-r
ggplot() +
geom_segment(data=segment(ddata), aes(x=x, y=y, xend=xend, yend=yend)) +
geom_text(data=ddata$labels, aes(x=x, y=y, label=label, hjust=0, color = Group), size=3, fontface = "bold", show.legend = FALSE) +
coord_flip() +
scale_y_reverse(expand=c(0.2, 8)) +
scale_colour_manual(values=c("grey70", "#41bf1b", "#6d6d6d", "#f42a35", "#f4f920")) +
theme_bw() +
theme(panel.border = element_blank(),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      axis.line = element_blank(),
      axis.title = element_blank(),
      axis.text = element_blank(),
      axis.ticks = element_blank(),
      legend.title = element_blank())
ggsave("methylation_atlas.betasites.hcluster.dendrogram.pdf")

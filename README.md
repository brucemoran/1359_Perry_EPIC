# 1359_Perry_EPIC Array Analysis

Run the EPIC array analysis for Silva et al. (in submission)

## Getting started

In a Linux environment simply issue the commands:

```
git clone https://github.com/brucemoran/1359_Perry_EPIC
cd 1359_Perry_EPIC
sh scripts/run_RnBeads-probeEMR-corrMatrix.sh \
   $(readlink -e ./1359_Perry_EPIC) \
   $(readlink -e ./1359_Perry_EPIC/scripts)
   "1359_Perry_EPIC"
   $(readlink -e ./1359_Perry_EPIC/data/sample_annotation.csv)
```

This downloads the raw IDAT files (~440MB), EPIC array annotation (~570MB), and runs RnBeads, and our full analysis.

### Prerequisites

You will need a copy of R installed, as well as required packages. These can be installed thus:

```
install.packages("BiocManager")
library("BiocManager")
BiocManager::install("GenomicRanges",
                     "RnBeads",
                     "FDb.InfiniumMethylation.hg19",
                     "tidyverse",
                     "GGally",
                     "magrittr",
                     "pheatmap",
                     "mgcv",
                     "ggplot2",
                     "readxl",
                     "VennDiagram",
                     "reshape2")

```

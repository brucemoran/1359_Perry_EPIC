# 1359_Perry_EPIC Array Analysis

Run the EPIC array analysis for Silva et al. (in submission)

## Getting started

In a Linux environment simply issue the commands (N.B. for OSX you can use 'greadlink' via Homebrew):

```
git clone https://github.com/brucemoran/1359_Perry_EPIC
sh scripts/run_RnBeads-probeEMR-corrMatrix.sh \
   "$(readlink -e ./1359_Perry_EPIC)" \
   "$(readlink -e ./1359_Perry_EPIC/scripts)" \
   "1359_Perry_EPIC" \
   "$(readlink -e ./1359_Perry_EPIC/data/sample_annotation.csv)"
```

This downloads the raw IDAT files (~440MB), EPIC array annotation (~570MB), and runs RnBeads, and our full analysis. To run the R analysis, and not the full RnBeads analysis from IDATs, include a fifth input to the shell script

### Prerequisites

You will need a copy of R installed (3.5.0 working), as well as required packages. These can be installed using the the script specified below. NB that this uses BiocManager, which may not be available to older versions of R.

```
Rscript scripts/prerequisites.R
```

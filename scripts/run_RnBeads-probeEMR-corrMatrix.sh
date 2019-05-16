#! /bin/bash

##this script runs the component parts of the analysis of data in GSEXYZ
##first part downloads IDAT files, second runs RnBeads, third runs probe, EMR ananlysis, fourth other plots

BASEDIR=$1
SCRIPTDIR=$2
RUNID=$3
ANNOFILE=$4
RONLY=$5

OUTDIR="$BASEDIR/analysis"
DATADIR="$BASEDIR/data"
IDATDIR="$DATADIR/IDAT_FILES"
ANNODIR="$DATADIR/EPIC_ANNO"
ANNOFILE="$DATADIR/sample_annotation.csv"

#mkdirs, test if made
if [[ ! -d "$OUTDIR" ]];then
  mkdir -p "$OUTDIR"
fi

######################
## 0: Download IDAT ##
######################
if [[ $RONLY == "" ]];then
  if [[ ! -d "$IDATDIR" ]];then
    echo "Making $IDATDIR and downloading data..."
    mkdir -p "$IDATDIR"
    for x in {390..405};do
      wget -P "$IDATDIR" \
        ftp://ftp.ncbi.nlm.nih.gov/geo/samples/GSM3362nnn/GSM33624${x}/suppl/*
    done
  else
    echo -e "Found $IDATDIR, delete to rerun download of data\nRunning analysis..."
  fi
else
  echo "IDAT downloads skipped\nDownloading BMIQ data..."
  wget -P "$OUTDIR/RnBeads" \
    ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE119nnn/GSE119260/suppl/*BMIQ*
  mv "$OUTDIR/RnBeads/*BMIQ*" "$OUTDIR/RnBeads/$RUNID.methylation_per_sample.sites.tsv"
fi

##########################
## 0.1: EPIC Annotation ##
##########################
if [[ -e "$ANNOFILE" ]]; then
  echo "Generating EPIC annotation BED file..."
  sh "$SCRIPTDIR/EPIC_annotation_maker.sh" "$ANNODIR"
else
  echo "Found EPIC annotation BED file $ANNOFILE, delete to rerun, proceeding..."
fi

####################
## 2: RnBeads Run ##
####################
if [[ $RONLY == "" ]];then
echo "Running RnBeads..."
  Rscript --vanilla "$SCRIPTDIR/RnBeads_run.paper.R" \
    "$BASEDIR" \
    "$IDATDIR" \
    "$RUNID" \
    "$ANNOFILE"
fi

########################
## 3: Probe, EMR sets ##
########################
echo "Generating probe, EMR sets per sample-type and patient..."
Rscript --vanilla "$SCRIPTDIR/set_patient_EMR.callerv2.R" \
  "$BASEDIR" \
  "$SCRIPTDIR" \
  "$RUNID" \
  "$DATADIR/EPIC_ANNO/MethylationEPIC_v-1-0_B4.Legacy_B2.bed" \
  "$OUTDIR/RnBeads/$RUNID.methylation_per_sample.sites.tsv"

####################
## 4: Other plots ##
####################
echo "Plotting correlation matrix, GAM with 450K"
Rscript --vanilla "$SCRIPTDIR/correlation_matrix_liquids-450K_GAM.R" \
  "$BASEDIR" \
  "$SCRIPTDIR" \
  "$RUNID"

###################
## 5: meth_atlas ##
###################
echo "Running meth_atlas"
sh "$SCRIPTDIR/meth_atlas_analysis.sh" "$BASEDIR" "$OUTDIR/RnBeads/$RUNID.methylation_per_sample.sites.tsv"
Rscript --vanilla "$SCRIPTDIR/meth_atlas_analysis.R" \
  "$BASEDIR" \
  "$OUTDIR/RnBeads/$RUNID.methylation_per_sample.sites.tsv" \
  "$DATADIR/meth_atlas/reference_atlas.csv"

#! /bin/bash

##this script runs the component parts of the analysis of data in GSEXYZ
##first part downloads IDAT files, second runs RnBeads, third runs probe, EMR ananlysis, fourth other plots

BASEDIR=$1
SCRIPTDIR=$2
RUNID=$3
ANNOFILE=$4

ANALDIR="$BASEDIR""/analysis"
IDATDIR="$BASEDIR""/IDAT_FILES"
ANNODIR="$BASEDIR""/EPIC_ANNO"
ANNOFILE="$IDATDIR""/sample_annotations/sample_annotation.tsv"

#mkdirs, test if made
if [[ ! -d "$IDATDIR" ]];then
  echo "$IDATDIR was not found"
  exit 127;
fi
if [[ ! -d "$ANALDIR" ]];then
  mkdir -p "$ANALDIR"
fi

######################
## 0: Download IDAT ##
######################

##########################
## 0.1: EPIC Annotation ##
##########################
echo "Generating EPIC annotation BED file..."
sh $SCRIPTDIR"/Epic_annotation_maker.sh" "$ANNODIR"

####################
## 2: RnBeads Run ##
####################
echo "Running RnBeads..."
Rscript --vanilla $SCRIPTDIR"/RnBeads_run.paper.R" \
  "$BASEDIR" \
  "$IDATDIR/ALL" \
  $RUNID \
  "$ANNOFILE"

########################
## 3: Probe, EMR sets ##
########################
echo "Generating probe, EMR sets per sample-type and patient"
Rscript --vanilla $SCRIPTDIR"/set_patient_EMR.callerv2.R" \
  "$BASEDIR" \
  "$SCRIPTDIR" \
  "$RUNID" \
  "$BASEDIR/EPIC_ANNO/MethylationEPIC_v-1-0_B4.Legacy_B2.bed" \
  "$BASEDIR/analysis/$RUNID/$RUNID.methylation_per_sample.sites.tsv"

####################
## 4: Other plots ##
####################
echo "Plotting correlation matrix, GAM with 450K"
Rscript --vanilla $SCRIPTDIR"/correlation_matrix_liquids-450K_GAM.R" \
  "$BASEDIR" \
  "$SCRIPTDIR" \
  "$RUNID"

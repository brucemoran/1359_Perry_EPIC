#! /bin/bash

##analysis using meth_atlas
BASEDIR="$1"
BETASITES="$2"
WORKDIR="$BASEDIR/analysis/meth_atlas"

##clone
cd "$BASEDIR/data"
git clone https://github.com/nloyfer/meth_atlas
git reset --hard "a1036496bafccdfcc001032460e3c18788ee653b"

##make pools of probes for input, parse from main table
cd $WORKDIR
SETHDIR="$BASEDIR/analysis/probe_EMR/probes/set"
echo "Running: meth atlas analysis"
cut -f 1 "$BETASITES" | tail -n+2 > all.probes
for TYPE in Plasma Urine Tumour Benign; do
  echo -e "\t"$TYPE
  cat "$SETHDIR/"${TYPE}.hy*.probes.txt | sort > "${TYPE}.tmp.txt"

  perl "$BASEDIR/scripts/perl_grep.pl" "${TYPE}.tmp.txt" "$BETASITES" "$TYPE" > "${TYPE}.EMPs.csv"

  ./deconvolve.py \
    --atlas_path reference_atlas.csv \
    "${TYPE}.EMPs.csv"

  perl "$BASEDIR/scripts/perl_grep.pl" all.probes "$BETASITES" "$TYPE" > "${TYPE}.all.csv"

  ./deconvolve.py \
    --atlas_path reference_atlas.csv \
    "${TYPE}.all.csv"
done

#!/bin/sh

set -e

SD=$(dirname $0)

ORIGINAL_GFA=$1
CONSENSUS_GAF=$2
SUPP=$3
WD=$4
THREADS=$5

mkdir -p $WD

>&2 echo "[$(date)] Augmenting"
vg augment --include-paths --min-coverage 1 --gaf $ORIGINAL_GFA $CONSENSUS_GAF | vg view - > $WD/augmented-pass0.gfa

>&2 echo "[$(date)] Cleaning augmentation"
python3 $SD/remove_low_supported_vertices.py $WD/augmented-pass0.gfa $CONSENSUS_GAF $SUPP 2> $WD/cleaning-pass1.log | vg mod --unchop -

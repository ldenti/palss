#!/bin/sh

set -e

SD=$(dirname $0)

ORIGINAL_GFA=$1
CONSENSUS_GAF=$2
SUPP=$3
OGFA=$4
OGAF=$5

>&2 echo "[$(date)] Augmenting graph"
vg augment --include-paths --min-coverage 1 --gaf $ORIGINAL_GFA $CONSENSUS_GAF | vg view - > $OGFA.unclean

>&2 echo "[$(date)] Cleaning graph"
python3 $SD/clean.py --supp $SUPP $OGFA.unclean $CONSENSUS_GAF --gaf $OGAF | vg mod --unchop - > $OGFA


# >&2 echo "[$(date)] GraphAligner to both graphs (original and augmented)"
# cut -f1,17 $CONSENSUS_GAF | sed "s/^/>/" | sed "s/\tqs:Z:/\n/g" > $WD/consensus.fa
# GraphAligner --graph $ORIGINAL_GFA --reads $WD/consensus.fa --alignments-out $WD/consensus_to_original.gaf --preset vg --threads $(($THREADS/2)) &> $WD/consensus_to_original.log &
# GraphAligner --graph $WD/augmented.gfa --reads $WD/consensus.fa --alignments-out $WD/consensus_to_augmented.gaf --preset vg --threads $(($THREADS/2)) &> $WD/consensus_to_augmented.log

# wait

# >&2 echo "[$(date)] Cleaning consensus"
# python3 $SD/filter_consensus.py $CONSENSUS_GAF $WD/consensus_to_original.gaf $WD/consensus_to_augmented.gaf > $WD/filtered_consensus.gaf

# >&2 echo "[$(date)] Augmenting graph (II)"
# vg augment --include-paths --min-coverage 1 --gaf $ORIGINAL_GFA $WD/filtered_consensus.gaf | vg view - > $WD/augmented.pass2.gfa

# >&2 echo "[$(date)] Cleaning augmentation"
# python3 $SD/clean_augment.py $WD/augmented.pass2.gfa $CONSENSUS_GAF $SUPP | vg mod --unchop -

>&2 echo "[$(date)] Done"

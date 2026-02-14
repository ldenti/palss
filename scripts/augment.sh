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
vg augment --include-paths --min-coverage 1 --gaf $ORIGINAL_GFA $CONSENSUS_GAF | vg view - > $WD/augmented.gfa

>&2 echo "[$(date)] Cleaning augmentation (1)"
python3 $SD/clean_augment.py 1 $WD/augmented.gfa $CONSENSUS_GAF $SUPP > $WD/augmented.clean.gfa 2> $WD/cleaning-pass1.log

>&2 echo "[$(date)] Unchopping (1)"
vg mod --unchop $WD/augmented.clean.gfa > $WD/augmented.clean.unchopped.gfa

>&2 echo "[$(date)] Converting GAF to FASTA"
cut -f1,17 $CONSENSUS_GAF | sed "s/^/>/" | sed "s/\tqs:Z:/\n/g" > $WD/consensus.fa

>&2 echo "[$(date)] GraphAligner to original GFA"
GraphAligner --graph $ORIGINAL_GFA --reads $WD/consensus.fa --alignments-out $WD/consensus_to_original.gaf --preset vg --threads $THREADS &> $WD/graphaligner_to_original.log

>&2 echo "[$(date)] GraphAligner to augmented GFA"
GraphAligner --graph $WD/augmented.clean.unchopped.gfa --reads $WD/consensus.fa --alignments-out $WD/consensus_to_augmented.gaf --preset vg --threads $THREADS &> $WD/graphaligner_to_augmented.log

>&2 echo "[$(date)] Cleaning augmentation (2)"
python3 $SD/clean_augment.py 2 $WD/augmented.clean.gfa $WD/augmented.clean.unchopped.gfa $CONSENSUS_GAF $WD/consensus_to_original.gaf $WD/consensus_to_augmented.gaf 2 $WD/resulting_consensus.gaf > $WD/augmented.clean.pass2.gfa 2> $WD/cleaning-pass2.log

>&2 echo "[$(date)] Unchopping (2)"
vg mod --unchop $WD/augmented.clean.pass2.gfa

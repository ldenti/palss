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
python3 $SD/remove_low_supported_vertices.py $WD/augmented-pass0.gfa $CONSENSUS_GAF $SUPP 2> $WD/cleaning-pass1.log | vg mod --unchop - > $WD/augmented-pass1.gfa

>&2 echo "[$(date)] Converting GAF to FASTA"
cut -f1,17 $CONSENSUS_GAF | sed "s/^/>/" | sed "s/\tqs:Z:/\n/g" > $CONSENSUS_GAF.fa

>&2 echo "[$(date)] GraphAligner to original GFA"
GraphAligner --graph $ORIGINAL_GFA --reads $CONSENSUS_GAF.fa --alignments-out $WD/consensus_to_original.gaf --preset vg --threads $THREADS &> $WD/graphaligner_to_original.log

>&2 echo "[$(date)] GraphAligner to augmented GFA"
GraphAligner --graph $WD/augmented-pass1.gfa --reads $CONSENSUS_GAF.fa --alignments-out $WD/consensus_to_augmented.gaf --preset vg --threads $THREADS &> $WD/graphaligner_to_augmented.log

>&2 echo "[$(date)] Selecting consensus"
python3 $SD/select_consensus.py $WD/consensus_to_original.gaf $WD/consensus_to_augmented.gaf $CONSENSUS_GAF > $WD/resulting_consensus.gaf 2> $WD/consensus_selection.log

>&2 echo "[$(date)] Augmenting (2)"
vg augment --include-paths --min-coverage 1 --gaf $ORIGINAL_GFA $WD/resulting_consensus.gaf | vg view - > $WD/augmented-pass2.gfa

>&2 echo "[$(date)] Cleaning augmentation"
python3 $SD/remove_low_supported_vertices.py $WD/augmented-pass2.gfa $CONSENSUS_GAF $SUPP 2> $WD/cleaning-pass2.log | vg mod --unchop -

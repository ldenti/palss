#!/bin/sh

set -e

SD=$(dirname $0)

ORIGINAL_GFA=$1
CONSENSUS_GAF=$2
SUPP=$3
WD=$4
THREADS=$5
MODE=$6

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

# Thanks to --include-paths, we can have a mapping between vertices and consensus. This allows to compute the "real" read support of each novel vertex. However, by doing this, starting/ending vertices of alignments are split to have "paths" and this creates issues with graphaligner since it cannot place seeds on short vertices. We cannot unchop the graph here as done during pass1 since we might lose paths and not being able to use support again. Therefore, we might: (1) compare new alignments to older one. If are worse, we keep vertices along that consensus (taken from paths in graph) anyway or (2) unchop the graph and not filter again by support (but we may end up needing (1) anyway). Going (1) in clean_unused.py for now

# XXX: however, it seems that doing (1) is not bringing much

>&2 echo "[$(date)] Augmenting (2)"
vg augment --include-paths --min-coverage 1 --gaf $ORIGINAL_GFA $WD/resulting_consensus.gaf | vg view - > $WD/augmented-pass2.gfa

>&2 echo "[$(date)] Cleaning (2)"
python3 $SD/clean_gfa.py $WD/resulting_consensus.gaf $WD/augmented-pass2.gfa | vg mod --unchop - > $WD/augmented-pass2.clean.gfa

>&2 echo "[$(date)] Converting GAF to FASTA"
cut -f1,17 $WD/resulting_consensus.gaf | sed "s/^/>/" | sed "s/\tqs:Z:/\n/g" > $WD/resulting_consensus.gaf.fa

>&2 echo "[$(date)] GraphAligner to augmented GFA (2)"
GraphAligner --graph $WD/augmented-pass2.clean.gfa --reads $WD/resulting_consensus.gaf.fa --alignments-out $WD/resulting_consensus_to_augmented.gaf --preset vg --threads $THREADS &> $WD/graphaligner_to_augmented.log

>&2 echo "[$(date)] Cleaning augmentation (2)"
python3 $SD/clean_ununsed.py $WD/augmented-pass2.clean.gfa $WD/resulting_consensus.gaf $WD/consensus_to_augmented.gaf $WD/resulting_consensus_to_augmented.gaf $SUPP $MODE 2> $WD/cleaning-pass2.log | vg mod --unchop -

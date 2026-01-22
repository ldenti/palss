#!/bin/sh

set -e
eval "$(~/miniforge3/bin/conda shell.zsh hook)"

SD=$(dirname $0)

WD=$1
THREADS=$2

reads=$WD/HG04199.pb-10x.fq.gz # XXX: corrected or not?

conda activate graphaligner

for n in 8
do
    mkdir -p $WD/graph-alignments/$n/
    # FULL
    GraphAligner --graph $WD/$n/pangenome-full.gfa --reads $reads --alignments-out $WD/graph-alignments/$n/full.gaf --preset vg --threads $THREADS
    # 1OUT
    GraphAligner --graph $WD/$n/pangenome.gfa --reads $reads --alignments-out $WD/graph-alignments/$n/1out.gaf --preset vg --threads $THREADS
    # PALSS
    for gfa in $WD/$n/palss.*/pangenome-augmented.*.w*.gfa
    do
	parameters=$(basename $gfa .gfa | cut -d"." -f2-)
	GraphAligner --graph $WD/$n/pangenome.gfa --reads $reads --alignments-out $WD/graph-alignments/$n/1out-augmented.$parameters.gaf --preset vg --threads $THREADS
    done
    for gfa in $WD/$n/palss-full.*/pangenome-augmented.*.w*.gfa
    do
	parameters=$(basename $gfa .gfa | cut -d"." -f2-)
	GraphAligner --graph $WD/$n/pangenome.gfa --reads $reads --alignments-out $WD/graph-alignments/$n/full-augmented.$parameters.gaf --preset vg --threads $THREADS
    done
    # MG-CACTUS
    # GraphAligner --graph $WD/$n/pangenome.gfa --reads $reads --alignments-out $WD/graph-alignments/$n/1out.gaf --preset vg --threads $THREADS
done

conda deactivate

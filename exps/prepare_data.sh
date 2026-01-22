#!/bin/sh

set -xe

SD=$(dirname $0)

GBZ=$1
WD=$2

Ns=(1 8 32 64)

mkdir -p $WD

N=${Ns[${#Ns[@]} - 1]}

vg gbwt --gbz-input --samples --list-names $GBZ | grep -v -P "CHM13|GRCh38" | head -n $((N+1)) > $WD/all_samples.txt

# === NEW SAMPLE ===
sample=$(head -n 1 $WD/all_samples.txt)
echo $sample > $WD/$sample.txt
vg paths --paths-by "$sample#1" --extract-fasta --xg $GBZ > $WD/$sample.hap1.fa
vg paths --paths-by "$sample#2" --extract-fasta --xg $GBZ > $WD/$sample.hap2.fa

# === OTHER SAMPLES ===
for n in "${Ns[@]}"
do
    echo "=== $n ==="
    mkdir -p $WD/$n/
    echo "CHM13" > $WD/$n/known_sample.txt
    tail -n +2 $WD/all_samples.txt | head -n $n >> $WD/$n/known_sample.txt
    $SD/utils/extract_subgraph $GBZ $WD/$n/known_sample.txt > $WD/$n/pangenome.gfa
    vg gbwt --gbz-format -g $WD/$n/pangenome.gbz --xg-name $WD/$n/pangenome.gfa --index-paths
done

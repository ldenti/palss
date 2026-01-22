#!/bin/sh

set -e

SD=$(dirname $0)

GBZ=$1
WD=$2

Ns=(1 8) # 32 64)

mkdir -p $WD

N=${Ns[${#Ns[@]} - 1]}

echo "[$(date)] Getting samples"
vg gbwt --gbz-input --samples --list-names $GBZ | grep -v -P "CHM13|GRCh38" | head -n $((N+1)) > $WD/all_samples.txt

# === NEW SAMPLE ===
sample=$(head -n 1 $WD/all_samples.txt)
echo $sample > $WD/$sample.txt
echo "[$(date)] Extracting new sample (haplotype 1)"
vg paths --paths-by "$sample#1" --extract-fasta --xg $GBZ > $WD/$sample.hap1.fa
echo "[$(date)] Extracting new sample (haplotype 2)"
vg paths --paths-by "$sample#2" --extract-fasta --xg $GBZ > $WD/$sample.hap2.fa

# === OTHER SAMPLES ===
for n in "${Ns[@]}"
do
    echo "[$(date)] Extracting pangenome ($n)"
    mkdir -p $WD/$n/
    echo "CHM13" > $WD/$n/known_samples.txt
    tail -n +2 $WD/all_samples.txt | head -n $n >> $WD/$n/known_samples.txt
    $SD/utils/extract_subgraph $GBZ $WD/$n/known_samples.txt > $WD/$n/pangenome.gfa
    vg gbwt --gbz-format -g $WD/$n/pangenome.gbz --xg-name $WD/$n/pangenome.gfa --index-paths

    echo "[$(date)] Extracting full pangenome ($n)"
    echo "CHM13" > $WD/$n/all_samples.txt
    echo $sample >> $WD/$n/all_samples.txt
    tail -n +2 $WD/all_samples.txt | head -n $n >> $WD/$n/all_samples.txt
    $SD/utils/extract_subgraph $GBZ $WD/$n/all_samples.txt > $WD/$n/pangenome-full.gfa
    vg gbwt --gbz-format -g $WD/$n/pangenome-full.gbz --xg-name $WD/$n/pangenome-full.gfa --index-paths
done

echo "[$(date)] Done"

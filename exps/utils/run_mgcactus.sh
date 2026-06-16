#!/bin/sh

set -xe

FA=$1 # reference
GBZ=$2 # known pangenome
name=$3
HAP1=$4 # first assembled haplotype of new individual
HAP2=$5 # second assembled haplotype of new individual
WD=$6
threads=$7
mcgtime=$8

mkdir -p $WD

echo -e "# Haploid sample (reference):" > $WD/seqfile
echo -e "CHM13\t$FA" >> $WD/seqfile

echo "# Diploid samples:" >> $WD/seqfile
vg gbwt --gbz-input --samples --list-names $GBZ | grep -Pv "CHM13|GRCh38" | while read idx
do
    vg paths --paths-by "$idx#1" --extract-fasta --xg $GBZ | sed "s/#/_/g" > $WD/$idx.h1.fa
    echo -e "$idx.1\t$WD/$idx.h1.fa" >> $WD/seqfile
    vg paths --paths-by "$idx#2" --extract-fasta --xg $GBZ | sed "s/#/_/g" > $WD/$idx.h2.fa
    echo -e "$idx.2\t$WD/$idx.h2.fa" >> $WD/seqfile
done

echo -e "$name.1\t$HAP1" >> $WD/seqfile
echo -e "$name.2\t$HAP2" >> $WD/seqfile

>&2 echo "[$(date)] Running minigraph-cactus"
rm -rf $WD/JOBSTORE
mkdir -p $WD/tmptmptmp
/usr/bin/time -vo ${mcgtime} cactus-pangenome $WD/JOBSTORE $WD/seqfile --outDir $WD --workDir $WD/tmptmptmp --outName pangenome --reference CHM13 --maxCores $threads --logColors False --logFile $WD/cactus-pangenome.log --binariesMode singularity

zcat $WD/pangenome.gfa.gz | head -1
zgrep "^S" $WD/pangenome.gfa.gz | cut -f1-3
zgrep -P "^[L|W|P]\t" $WD/pangenome.gfa.gz

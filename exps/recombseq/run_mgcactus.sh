#!/bin/sh

set -xe

FA=$1 # reference
VCF=$2 # VCF with known individuals
name=$3
HAP1=$4 # first assembled haplotype of new individual
HAP2=$5 # second assembled haplotype of new individual
WD=$6
threads=$7

mkdir -p $WD

echo -e "# Haploid sample (reference):" > $WD/seqfile
echo -e "CHM13\t$FA" >> $WD/seqfile

echo "# Diploid samples:" >> $WD/seqfile
bcftools view -h $VCF | tail -1 | cut -f10- | tr '\t' '\n' | while read idx
do
    bcftools consensus -s $idx -H1 --fasta-ref $FA $VCF > $WD/$idx.h1.fa
    echo -e "$idx.1\t$WD/$idx.h1.fa" >> $WD/seqfile
    bcftools consensus -s $idx -H2 --fasta-ref $FA $VCF > $WD/$idx.h2.fa
    echo -e "$idx.2\t$WD/$idx.h2.fa" >> $WD/seqfile
done

echo -e "$name.1\t$HAP1" >> $WD/seqfile
echo -e "$name.2\t$HAP2" >> $WD/seqfile

rm -rf $WD/JOBSTORE
/usr/bin/time -v cactus-pangenome $WD/JOBSTORE $WD/seqfile --outDir $WD --outName pangenome --reference CHM13

zcat $WD/pangenome.gfa.gz | head -1
zgrep "^S" $WD/pangenome.gfa.gz | cut -f1-3
zgrep -P "^[L|W|P]\t" $WD/pangenome.gfa.gz

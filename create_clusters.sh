#!/bin/sh

set -e

fq=$1
vg=$2
f=$3
odir=$4

mkdir -p $odir

i=0
while read line
do
    rnames=$(echo $line | cut -f4 -d' ')

    echo "Getting fastq ($i)..."
    echo $rnames | tr "," "\n" > $odir/$i
    samtools faidx --region-file $odir/$i $fq > $odir/$i.fa

    echo "Getting graph ($i)..."
    nodes=$(echo $line | cut -f5 -d' ')
    cl=$(echo $nodes | sed "s/,/ --subgraph /g")
    vg mod --subgraph $cl $vg | vg view - > $odir/$i.gfa

    i=$((i+1))
done < $f
#!/bin/sh

set -xe

fa=$1
bed=$2
wd=$3
truth1=$4
truth2=$5
vcfs="${@:5}"

mkdir -p $wd

for truth in $truth1 $truth2
do
    tname=$(basename $truth | cut -d"-" -f1)

    # Force GT to 1|1
    bcftools +setGT $truth -- -n c:"0|1" -t a | bgzip -c > $wd/truth.vcf.gz
    tabix -p vcf $wd/truth.vcf.gz

    for vcf in $vcfs
    do
	cname=$(basename $vcf | cut -d"-" -f1)
	subwd=$wd/$cname-against-$tname/
	mkdir -p $subwd

	# Force GT to 1|1
	bcftools +setGT $vcf -- -n c:"0|1" -t a | bgzip -c > $wd/call.vcf.gz
	tabix -fp vcf $wd/call.vcf.gz
	
	hap.py --leftshift $wd/truth.vcf.gz $wd/call.vcf.gz -o $subwd/result-nobed -r $fa
	hap.py --leftshift $wd/truth.vcf.gz $wd/call.vcf.gz -o $subwd/result-bed -r $fa -f $bed
    done
done

rm $wd/truth.vcf.gz $wd/truth.vcf.gz.tbi
rm $wd/call.vcf.gz $wd/call.vcf.gz.tbi

#!/bin/sh

set -xe

fa=$1
bed=$2
wd=$3
threads=$4
truth1=$5
truth2=$6
vcfs="${@:7}"

mkdir -p $wd

for truth in $truth1 $truth2
do
    tname=$(basename $truth | cut -d"-" -f1)

    for vcf in $vcfs
    do
	cname=$(basename $vcf | cut -d"-" -f1)

	mkdir -p $wd/$cname-against-$tname

	subwd=$wd/$cname-against-$tname/conf
        truvari bench --passonly --pick ac --dup-to-ins --includebed $bed --reference $fa --base $truth --comp $vcf --output $subwd
        # truvari refine --reference $fa --regions $subwd/candidate.refine.bed --coords R --use-original-vcfs --threads $threads --align mafft $subwd
        # truvari ga4gh --input $subwd --output $subwd/ga4gh_with_refine # --with-refine


	subwd=$wd/$cname-against-$tname/full
        truvari bench --passonly --pick ac --dup-to-ins --reference $fa --base $truth --comp $vcf --output $subwd
        # truvari refine --reference $fa --regions $subwd/candidate.refine.bed --coords R --use-original-vcfs --threads $threads --align mafft $subwd
        # truvari ga4gh --input $subwd --output $subwd/ga4gh_with_refine # --with-refine
    done
done

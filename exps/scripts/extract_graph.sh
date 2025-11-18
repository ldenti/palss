#!/bin/sh

echo $(dirname $0)

gbz=$1
# vg view $gbz > $gfa
gfa=$2
n=$3
wd=$4

vg gbwt --gbz-input --samples --list-names $gbz | grep -Pv "CHM13|GRCh38" > $wd/all_samples.list

echo "CHM13" > $wd/samples.list
head -n $n all_samples.list >> $wd/samples.list

new_sample=$(head -n $((n+1)) $wd/all_samples.list | tail -1)
echo $new_sample > $wd/new_sample.txt

# python script does not remove edges not consistent with paths
python3 $(dirname $0)/remove_samples.py $gfa $wd/samples.list | vg mod --remove-non-path - > $wd/pangenome.vg

vg gbwt --gbz-format --graph-name $wd/pangenome.gbz --xg-name $wd/pangenome.vg --index-paths

vg paths --paths-by "$sample#1" --extract-fasta --xg $gbz > $wd/new_sample.hap1.fa
vg paths --paths-by "$sample#2" --extract-fasta --xg $gbz > $wd/new_sample.hap2.fa

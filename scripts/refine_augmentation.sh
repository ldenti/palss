#!/bin/sh

set -e

SD=$(dirname $0)
HIFIASM_BIN="/home/ld/code/hifiasm/hifiasm"

GFA=$1
FA=$2
SFS=$3
WD=$4
NTH=$5
ID=$6

mkdir -p $WD

>&2 echo "[$(date)] Extracting unanchored reads"
grep -P "^1|^2|^3" $SFS | cut -f2 | sort -u > $WD/reads_with_unanchored.list
python3 $SD/subfa.py $FA $WD/reads_with_unanchored.list > $WD/reads_with_unanchored.fa

>&2 echo "[$(date)] Assembling reads"
$HIFIASM_BIN -r0 -o $WD/reads_with_unanchored -t$NTH $WD/reads_with_unanchored.fa 2> $WD/hifiasm.log
awk '/^S/{print ">"$2;print $3}' $WD/reads_with_unanchored.bp.p_ctg.gfa > $WD/reads_with_unanchored.bp.p_ctg.fa
# python3 $SD/cut_contigs.py $WD/reads_with_unanchored.bp.p_ctg.gfa $CUT > $WD/reads_with_unanchored.bp.p_ctg.fa

>&2 echo "[$(date)] Aligning contigs"
GraphAligner --graph $GFA --reads $WD/reads_with_unanchored.bp.p_ctg.fa --alignments-out $WD/reads_with_unanchored.bp.p_ctg.gaf --preset vg --threads $NTH &> $WD/graphaligner.log
$SD/add_cs $GFA $WD/reads_with_unanchored.bp.p_ctg.fa $WD/reads_with_unanchored.bp.p_ctg.gaf > $WD/reads_with_unanchored.bp.p_ctg.wcs.gaf
python3 $SD/filter_gaf.py $WD/reads_with_unanchored.bp.p_ctg.wcs.gaf $ID > $WD/reads_with_unanchored.bp.p_ctg.wcs.selected.gaf

>&2 echo "[$(date)] Augmenting"
vg augment --min-coverage 1 --gaf $GFA $WD/reads_with_unanchored.bp.p_ctg.wcs.selected.gaf | vg mod --unchop - | vg view -

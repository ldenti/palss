#!/bin/sh

set -xe

SD=$(dirname $0)

FA1=$1
FA2=$2
SAMPLE_FQ=$3
OUTPUT_FQ=$4
COV=$5

COV=$((COV / 2.0))

mkdir -p $OUTPUT_FQ.tmp
pbsim --id-prefix "S1_" --prefix $OUTPUT_FQ.tmp/hap1 --strategy wgs --method sample --sample $SAMPLE_FQ --depth $COV --genome $FA1
pbsim --id-prefix "S2_" --prefix $OUTPUT_FQ.tmp/hap2 --strategy wgs --method sample --sample $SAMPLE_FQ --depth $COV --genome $FA2

# XXX: this might crash if we have too many .fq.gz
# XXX: we could parallelize this
zcat $OUTPUT_FQ.tmp/*.fq.gz | python3 $SD/utils/remove_n.py | gzip -c > $OUTPUT_FQ
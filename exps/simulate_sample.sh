#!/bin/sh

set -e
eval "$(~/miniforge3/bin/conda shell.zsh hook)"

SD=$(dirname $0)

REF=$1
FA1=$2
FA2=$3
SAMPLE_FQ=$4
COV=$5
THREADS=$6
OUTPUT_PREFIX=$7

COV=$(echo "scale=1; $COV/2" | bc)

mkdir -p $OUTPUT_PREFIX-pbsim3.tmp

conda activate pbsim3
echo "[$(date)] Simulating from haplotype 1"
pbsim --id-prefix "S1_" --prefix $OUTPUT_PREFIX-pbsim3.tmp/hap1 --strategy wgs --method sample --sample $SAMPLE_FQ --depth $COV --genome $FA1 2> $OUTPUT_PREFIX-pbsim3.hap1.log
echo "[$(date)] Simulating from haplotype 2"
pbsim --id-prefix "S2_" --prefix $OUTPUT_PREFIX-pbsim3.tmp/hap2 --strategy wgs --method sample --sample $SAMPLE_FQ --depth $COV --genome $FA2 2> $OUTPUT_PREFIX-pbsim3.hap1.log
conda deactivate

echo "[$(date)] Cleaning reads using $THREADS threads"
i=0
for fq in $OUTPUT_PREFIX-pbsim3.tmp/*.fq.gz
do
    python3 $SD/utils/remove_n.py $fq | gzip -c > $fq.clean &
    i=$((i+1))
    if (( i % $THREADS == 0 )); then
        wait
    fi
done

echo "[$(date)] Merging reads"
# XXX: this might crash if we have too many .fq.gz
cat $OUTPUT_PREFIX-pbsim3.tmp/*.clean > $OUTPUT_PREFIX.fq.gz
rm $OUTPUT_PREFIX-pbsim3.tmp/*.clean

echo "[$(date)] Correcting reads"
conda activate hifiasm
/usr/bin/time -vo $OUTPUT_PREFIX.hifiasm.time hifiasm -t$THREADS --write-ec --bin-only $OUTPUT_PREFIX.fq.gz -o $OUTPUT_PREFIX 2> $OUTPUT_PREFIX.hifiasm.log
conda deactivate

conda activate minimap2
echo "[$(date)] Aligning reads"
minimap2 -t$THREADS --MD -ax map-hifi --eqx $REF $OUTPUT_PREFIX.fq.gz | samtools view -bS | samtools sort > $OUTPUT_PREFIX.bam
samtools index $OUTPUT_PREFIX.bam
echo "[$(date)] Aligning corrected reads"
minimap2 -t$THREADS --MD -ax map-hifi --eqx $REF $OUTPUT_PREFIX.ec.fa | samtools view -bS | samtools sort > $OUTPUT_PREFIX.ec.bam
samtools index $OUTPUT_PREFIX.ec.bam
conda deactivate

echo "Sample: $OUTPUT_PREFIX.fq.gz"
echo "Corrected sample: $OUTPUT_PREFIX.ec.fa"
echo "[$(date)] Done."

### Prerequisites
Software/libraries:
* conda
* snakemake
* vg
* bedtools
* python3-pandas, python3-matplotlib, python3-seaborn, python3-pysam


### Data
```
# Get reference
wget https://s3-us-west-2.amazonaws.com/human-pangenomics/T2T/CHM13/assemblies/analysis_set/chm13v2.0.fa.gz
gunzip chm13v2.0.fa.gz

# Get pangenome graph
wget https://s3-us-west-2.amazonaws.com/human-pangenomics/pangenomes/freeze/release2/minigraph-cactus/hprc-v2.0-mc-chm13.gbz
vg view hprc-v2.0-mc-chm13.gbz > hprc-v2.0-mc-chm13.gfa

# Get NA12878 reads
wget https://42basepairs.com/download/s3/platinum-pedigree-data/data/hifi/mapped/CHM13/NA12878.CHM13.haplotagged.bam
wget https://42basepairs.com/download/s3/platinum-pedigree-data/data/hifi/mapped/CHM13/NA12878.CHM13.haplotagged.bam.bai
samtools view -b -s 42.1 NA12878.CHM13.haplotagged.bam | samtools fastq > NA12878.10x.fq

# Get NA12878 assembly
wget https://42basepairs.com/download/s3/platinum-pedigree-data/assemblies/NA12878/verkko/1.3.1/assembly.haplotype1.fasta
wget https://42basepairs.com/download/s3/platinum-pedigree-data/assemblies/NA12878/verkko/1.3.1/assembly.haplotype2.fasta
cat assembly.haplotype1.fasta assembly.haplotype2.fasta > assembly.haplotypes.fasta

# Get complex regions
wget https://hgdownload.soe.ucsc.edu/gbdb/hs1/bbi/simpleRepeat.bb
wget https://hgdownload.gi.ucsc.edu/admin/exe/linux.x86_64/bigBedToBed
chmod +x bigBedToBed
./bigBedToBed simpleRepeat.bb simpleRepeat.bed
python3 utils/bed_fix.py data/genbank2chr.info simpleRepeat.bed | sort -k1,1 -k2,2n > simpleRepeat.fixed.bed
wget https://hgdownload.soe.ucsc.edu/gbdb/hs1/sedefSegDups/sedefSegDups.bb
bigBedToBed sedefSegDups.bb sedefSegDups.bed
cat sedefSegDups.bed <(awk -v T=2500 '$3-$2 > T' simpleRepeat.fixed.bed | cut -f1-3) | sort -k1,1 -k2,2n | bedtools merge -i stdin > complex.bed
bedtools genomecov -i complex.bed -g chm13v2.0.fa.fai

# Compile some code we will need later
cd utils
g++ -I$PWD/../../include/ -L$PWD/../../build/zstd-prefix/src/zstd/lib -L$PWD/../../lib -Wl,-rpath,$PWD/../../lib -o extract_subgraph extract_subgraph.cpp -lgbwtgraph -lgbwt -lsdsl -fopenmp -lhandlegraph -lzstd -lcrypto
cd ..
```

### Experiment 1 - Sketching
```
snakemake -c16 -p -s sketch_analysis.smk --config fa=chm13v2.0.fa gbz=hprc-v2.0-mc-chm13.gbz fq=NA12878.10x.fq wd=palss-sketch_analysis [-n]
python3 ./plots/plot_sketch_analysis.py fa palss-sketch_analysis/d0.5/missed_regions.g15000.bed chm13v2.0.fa.fai complex.bed
python3 ./plots/plot_sketch_analysis.py fq palss-sketch_analysis/d0.5/reads.txt
```

### Experiment 2 - Simulated data
```
# Get single (or multi) chromosome(s) data
samtools faidx chm13v2.0.fa chr1 chr20 > chr1-20.fa
samtools faidx chr1-20.fa
wget https://s3-us-west-2.amazonaws.com/human-pangenomics/pangenomes/freeze/release2/minigraph-cactus/v2.0/hprc-v2.0-mc-chm13/hprc-v2.0-mc-chm13.chroms/chr1.vg
wget https://s3-us-west-2.amazonaws.com/human-pangenomics/pangenomes/freeze/release2/minigraph-cactus/v2.0/hprc-v2.0-mc-chm13/hprc-v2.0-mc-chm13.chroms/chr20.vg
vg combine chr1.vg chr20.vg > chr1-20.vg
vg gbwt --gbz-format -g chr1-20.gbz --xg-name chr1-20.vg --index-paths
egrep '^chr(1|20)[[:space:]]' complex.bed > chr1-20.bed

# edit config/exp2.yaml

snakemake -c32 -p --use-conda -s augment.smk --configfile config/exp2.yaml --keep-going [-n]
```

```
python3 ./plots/plot_exp2.py ${SMK_WD}/support.csv ${SMK_WD}/nm.csv
python3 ./plots/plot_exp2.supp.py ${SMK_WD}/support.csv ${SMK_WD}/nm.csv
python3 ./plots/plot_recall.py --nm 25 --steps 1 ${SMK_WD}/nm.csv
python3 ./plots/plot_exp2.wsample.py -x HG03742 -y HG01993 ${HG03742_WD}/support.csv ${HG03742_WD}/nm.csv ${HG01993_WD}/support.csv ${HG01993_WD}/nm.csv
python3 ./plots/build_pop_table.py ${SMK_WD}/n64/samples-full.list ./data/hprc_release2_sample_metadata.csv
```

### Experiment 3 - Real data
```
### PALSS ###
mkdir palss
\time -v ../palss sketch -@32 -w palss/sketching-wd hprc-v2.0-mc-chm13.gbz > palss/pangenome.skt
\time -v sh -c "LD_LIBRARY_PATH="$PWD/../lib" ../build/gbwtgraph-prefix/src/gbwtgraph/bin/gbz_extract --threads 32 --progress hprc-v2.0-mc-chm13.gbz | ../build/rb3-prefix/src/rb3/ropebwt3 build -t32 -Ld - > palss/pangenome.fmd"
\time -v ../palss sfs -@32 hprc-v2.0-mc-chm13.gbz palss/pangenome.skt palss/pangenome.fmd palss/reads.ec.fa > palss/specific_strings.new.txt
\time -v ../palss align -@32 hprc-v2.0-mc-chm13.gbz palss/specific_strings.txt > palss/consensus.gaf
\time -v ../palss augment -@4 -c -w palss/exp3-wd -g consensus-refined.gaf hprc-v2.0-mc-chm13.gbz palss/consensus.gaf > palss/pangenome-augmented.gfa

### Reads alignment ###
# to HPRCv2 graph
GraphAligner --graph hprc-v2.0-mc-chm13.gfa --reads NA12878.10x.fq --alignments-out reads_to_original.gaf --preset vg --threads 32

# to reference genome
minimap2 -ax map-hifi --MD --eqx -Y -t32 chm13v2.0.fa NA12878.10x.fq | samtools view -bS | samtools sort > reads_to_reference.bam
samtools index reads_to_reference.bam

# to assembly
minimap2 -ax map-hifi --MD --eqx -Y -t32 assembly.haplotypes.fasta NA12878.10x.fq | samtools view -bS | samtools sort > reads_to_assembly.bam 
samtools index reads_to_reference.bam

# to augmented graph
GraphAligner --graph palss/pangenome-augmented.gfa --reads NA12878.10x.fq --alignments-out reads_to_augmented.gaf --preset vg --threads 32

### Summarize NM ###
python3 utils/get_nm_table.py --gaf HPRCv2 reads_to_original.gaf palss reads_to_augmented.gaf --bam T2T-CHM13 reads_to_reference.bam Assembly reads_to_assembly.bam > nm.csv
python3 ./plots/plot_exp3.py nm.csv
```

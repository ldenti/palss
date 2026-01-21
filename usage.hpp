#ifndef PALSS_USAGE_H
#define PALSS_USAGE_H

static const char *const VERSION = "PALSS, v0.0.1";

static const char *const MAIN_USAGE =
    "Usage: palss [sketch|sfs|align] -h\n";

static const char *const SKETCH_USAGE_MESSAGE =
    "Usage: palss sketch [options] <graph.gfa> <paths.fa.fmd>\n"
    "Options:\n"
    "        -k <INT>     kmer size (default: 31, maximum: 32)\n"
    "        -d <FLOAT>   kmer density in [0,1] (default: 1.0, maximum: 32)\n"
    // "        -@ <INT>   set threads (default: 4)\n"
    "        -h           display this help and exit\n"
    "\n";

static const char *const SFS_USAGE_MESSAGE =
    "Usage: palss sfs [options] <graph.gbz> <graph.skt> <paths.fmd> <reads.fx>\n"
    "Options:\n"
    "        -a <INT>   number of anchors to check on each side of a specific "
    "string (default: 20)\n"
    "        -b <INT>   batch size (default: 10000)\n"
    "        -@ <INT>   threads (default: 4)\n"
    "        -h         display this help and exit\n"
    "\n";

static const char *const ALIGN_USAGE_MESSAGE =
    "Usage: palss align [options] <graph.gbz> <specific_strings.txt>\n"
    "Options:\n"
    // "        -@ <INT>   set threads (default: 4)\n"
    "        -h         display this help and exit\n"
    "\n";

// static const char *const KAN_USAGE_MESSAGE =
//     "Usage: palss kan [options] <graph.gfa.skt> <file.fx>\n"
//     "Options:\n"
//     "        -r         use this if input is FASTQ, so that we won't output a "
//     "BED file (default: "
//     "FASTA, BED output)\n"
//     "        -h         display this help and exit\n"
//     "\n";

#endif

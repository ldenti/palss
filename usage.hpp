#ifndef PALSS_USAGE_H
#define PALSS_USAGE_H

static const char *const VERSION = "PALSS, v0.0.1";

static const char *const MAIN_USAGE =
    "Usage: palss [sketch|search|augment] -h\n";

static const char *const SKETCH_USAGE_MESSAGE =
    "Usage: palss sketch [options] <graph.gfa> <paths.fa.fmd>\n"
    "Options:\n"
    "        -k <INT>   kmer size (default: 31, maximum: 32)\n"
    // "        -m <INT>   mmer size (default: 9, maximum: k-1)\n"
    // "        -@ <INT>   set threads (default: 4)\n"
    "        -b         big sketch for big graph (require 64GB+overhead, "
    "realistically <40GB)\n"
    "        -h         display this help and exit\n"
    "\n";

// static const char *const SEARCH_USAGE_MESSAGE =
//     "Usage: palss search [options] <graph.gfa> <graph.gfa.skt> <paths.fa.fmd> "
//     "<reads.fx>\n"
//     "Options:\n"
//     "        -p <STR>   load only paths containing this (default: CHM13)\n"
//     "        -s <INT>   number of anchors to check on each side of a specific "
//     "string (default: 20)\n"
//     "        -b <INT>   batch size (default: 10000)\n"
//     "        -@ <INT>   set threads (default: 4)\n"
//     "        -a         use all solid anchors (default: only anchors on "
//     "reference paths)\n"
//     "        -h         display this help and exit\n"
//     "\n";

// static const char *const AUGMENT_USAGE_MESSAGE =
//     "Usage: palss augment [options] <graph.gfa> <specific_strings.txt>\n"
//     "Options:\n"
//     "        -k <INT>   kmer size (default: 31, maximum: 32)\n"
//     "        -p <STR>   load only paths containing this (default: CHM13)\n"
//     "        -g <PATH>  output GAF (default: \"\", ie, no GAF output)\n"
//     /* "        -s <INT>   minimum alignment score (default: 0)\n" */
//     "        -w <INT>   minimum support per cluster (default: 2)\n"
//     /* "        -@ <INT>   set threads (default: 4)\n" */
//     "        -h         display this help and exit\n"
//     "\n";

static const char *const KAN_USAGE_MESSAGE =
    "Usage: palss kan [options] <graph.gfa.skt> <file.fx>\n"
    "Options:\n"
    "        -r         use this if input is FASTQ, so that we won't output a "
    "BED file (default: "
    "FASTA, BED output)\n"
    "        -h         display this help and exit\n"
    "\n";

#endif

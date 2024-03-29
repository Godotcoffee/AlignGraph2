### Please visit this [Github site][aligngraph2] for the latest AlignGraph2 version.

## Introduction

AlignGraph2 is the second version of [AlignGraph][aligngraph] for PacBio long reads. It extends and refines contigs assembled from  the long reads with a published genome similar to the sequencing genome.

## How to cite AlignGraph2?
If you use AlignGraph2, please cite the following paper:  
Huang, S., He, X., Wang, G., & Bao, E. (2021). AlignGraph2: similar genome-assisted reassembly pipeline for PacBio long reads. Briefings in Bioinformatics: [epub](https://doi.org/10.1093/bib/bbab022).

## Requirements

Linux OS 64-bit

GCC >= 4.8.5

CMake >= 3.1.0

Python >= 3.6 with Biopython

## Installation

```sh
git clone --recursive https://github.com/Godotcoffee/AlignGraph2.git
cd AlignGraph2 && python ./install.py
```

## Usage

### Quick Start

```sh
python ./AlignGraph2.py -r /path/to/read -c /path/to/contig -g /path/to/genome -o /path/to/output
```

### Mandatory

`-r, --read [fastq]`, long read file

`-c, --contig [fasta]`, contig file

`-g, --genome [fasta]`, similar genome file

`-o, --output [dir]`, output directory

### Options

`-m`, customized alignment algorithm mecat2ref+ (default: none)

`-b [int]`, size of similar genome blocks for mecat2ref+ (default: 200)

`--alpha [real]`, lower bound of k-mer scoring function for mecat2ref+ (default: 0.5)

`--beta [real]`, upper bound of k-mer scoring function for mecat2ref+ (default: 2.0)

`--delta [real]`, threshold for alignment scoring (default: 0.9)

`-v [int]`, coverage to filter alignments (default: 2)

`-k [int]`, size of k-mers in A-Bruijn graph (default: 14)

`--epsilon [int]`, distance to join two vertices in A-Bruijn graph (default: 10)

`-l [int]`, minimum path length for graph traversal (default: 50)

`-a [int]`, size of long read blocks for consensus (default: 10,000)

`-t [int]`, number of threads (default: 16)

[aligngraph]: https://github.com/baoe/AlignGraph
[aligngraph2]: https://github.com/huangs001/AlignGraph2

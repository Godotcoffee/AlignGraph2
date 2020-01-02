## Introduction

AlignGraph2 is a software improved from [AlignGraph][aligngraph] that 
joins contigs generated from long reads by reassembling them with 
guidance of a reference genome of a closely related species.

## Requirements

GCC >= 4.8.5

CMake >= 3.1.0

Python >= 3.6

## Installation

```sh
git clone https://github.com/Godotcoffee/AlignGraph2
cd AlignGraph2 && python ./install.py
```

## Usage

### Quick Start

```sh
python ./AlignGraph2.py -r /path/to/read -c /path/to/contig -R /path/to/reference -o /path/to/output
```

### Options
`-r, --read [fastq/fasta]`, reads file

`-c, --contig [fasta]`, contigs file

`-R, --ref [fasta]`, reference genome file

`-o, --output [dir]`, output directory

`-k [size]`, size of k-mer (default: 14)

`--ratio [threshold]`, threshold of solid k-mer set (default: 0.2)

`-t --thread`, number of thread (default: 16)

[aligngraph]: https://github.com/baoe/AlignGraph
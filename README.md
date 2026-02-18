# Kaspar et al., 2026

[![Python 3.9+](https://img.shields.io/badge/python-3.9+-blue.svg)](https://www.python.org/downloads/release/python-390/)
[![DOI](https://img.shields.io/badge/DOI-10.xxxx%2Fxxxxx-blue)](http://doi.org)

This repository contains the official implementation of the paper:

**"[Paper Title]"**
*Authors: Kašpar Tomáš, Vojtěch Čermák, Kateřina Adamusová, and Lukáš Fischer*
**Kaspar et al., 2026**
Published in: *Journal*

# NAPALM
NAPALM (**N**anopore **APA** **L**ocus **M**etrics) allows analysis of APA sites from Nanopore RNA sequencing data.
It uses annotation file in GTF format and one or multiple Nanopore aligned reads in BAM format.


## Installation

To run this program, make sure you have Python installed. Then, install the required dependencies.
## Usage
```
python napalm.py -b ./bam1.bam ./bam2.bam ./bam3.bam -g ./annotation.gtf -t 12 -min 5 -w 10
```

### Required
- `(--genomefile GENOMEFILE `
- `--bamfile [BAMFILE ...]` : Path(s) to one or more files of aligned sRNA-seq data in BAM format. Multiple files are separated by spaces. BAM files must match the reference genome given in `--genomefile`.

### Recommended
- `-m MERGE_WINDOW`


### Other options
- `-h` : Print a help message and then quit.
- `--version` : Print the version and then quit.

## Outputs

# komplexity

A command-line tool built in Rust to quickly calculate and/or mask low-complexity sequences from a FAST[A/Q] file. This uses the number of unique _k_-mers over a sequence divided by the length to assess complexity.

We've noticed that, for _k_=4, a normalized complexity score of < 0.55 suggests across a 64-120bp region strongly that the sequence is a low-complexity repeat sequence. 

## Installation

Assuming you have Rust and Cargo installed:

```sh
git clone https://github.com/eclarke/komplexity
cd komplexity
cargo install
```

If Cargo's default bin directory is in your path (usually `$HOME/.cargo/bin`), you can run this by just typing `kz`.

## Usage

This tool has two modes: measuring (default) and masking. 

### Measuring

Measuring mode reports the length of a sequence, the number of unique kmers in that sequence, and the normalized complexity (unique kmers / length). As this looks across the entire sequence, it is most appropriate for short reads (e.g. Illumina).

```sh
# For fastq files:
kz < seqs.fq
# For fasta files:
kz --fasta < seqs.fa
# Vary length of k:
kz -k5 < seqs.fq
```

The output is a list of sequence IDs from the fastq file with the length of the sequence, the number of unique _k_-mers in the sequence, and the number of unique _k_-mers/length (normalized complexity) in a tab-delimited format to stout.

### Masking

Masking mode (`--mask`) outputs the input sequences with low-complexity regions masked by Ns. The low-complexity regions are found using a sliding window (of size `--window_size`) over the sequence; the complexity of the subsequence is assessed as in the "measuring" mode. The threshold below which the sequence is masked is configurable through the `--threshold` parameter and should be a value between 0-1 (least to most masking). 

```sh
# For fastq files:
kz --mask < seqs.fq
# For fasta files:
kz --mask --fasta < seqs.fa
# Vary threshold:
kz --mask --threshold 0.7 < seqs.fq
# Vary window size:
kz --mask --window_size 64 < seqs.fq
```

## So what?


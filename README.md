## komplexity

A command-line tool built in Rust to quickly calculate and filter low-complexity sequences from a FASTQ file.



## Installation

Assuming you have Rust and Cargo installed:

```sh
git clone https://github.com/eclarke/komplexity
cd komplexity
cargo install
```

If Cargo's default bin directory is in your path (usually `$HOME/.cargo/bin`), you can run this by just typing `kz`.

## Usage

```sh
# For fastq files:
kz < seqs.fq
# For fasta files:
kz --fasta < seqs.fa
# Vary length of k:
kz -k5 < seqs.fq
```

For a given _k_, default 4, this tool calculates the number of unique _k_-mers in each sequence, normalized by the sequence length.

The output is a list of sequence IDs from the fastq file with the length of the sequence, the number of unique _k_-mers in the sequence, and the number of unique _k_-mers/length (normalized complexity) in a tab-delimited format to stout.

## So what?

We've noticed that, for _k_=4, a normalized complexity score of < 0.55 suggests strongly that the sequence is a low-complexity repeat sequence. This is true across many different samples, predominantly tested on Illumina reads with a length <= 126bp. 
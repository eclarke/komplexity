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
kz calculate < seqs.fq | kz filter > ids_passing_filter.txt
```

### kz-calculate
For a given _k_, this tool calculates the number of unique _k_-mers in each sequence, normalized by the sequence length.

The output is a list of sequence IDs from the fastq file with the normalized frequency.

### kz-filter
This subcommand takes the input of `kz calculate`, calculates a z-score for each sequence, and outputs only those IDs 
whose z-scores are higher than the threshold specified.
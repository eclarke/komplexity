extern crate bio;
extern crate clap;
extern crate fnv;
extern crate itertools;

use std::vec::Vec;
use std::fmt::{self, Display};
use std::io;

use bio::alphabets;
use bio::alphabets::RankTransform;
use bio::io::fastq;

use clap::{Arg, App};

use itertools::Itertools;

use fnv::FnvHashSet;

fn main() {

    let cli = App::new("Unique Kmers")
        .arg(Arg::with_name("input")
                 .short("i")
                 .long("input")
                 .default_value("-")
                 .takes_value(true)
                 .help("Path to FASTQ file (default: read from stdin)"))
        .arg(Arg::with_name("min_k")
                 .short("m")
                 .long("min_k")
                 .default_value("3")
                 .takes_value(true)
                 .help("Smallest k to search"))
        .arg(Arg::with_name("max_k")
                 .short("M")
                 .long("max_k")
                 .default_value("7")
                 .takes_value(true)
                 .help("Largest k to search"))
        .get_matches();

    let seq_fp = cli.value_of("input").unwrap();
    let min_k: u32 = cli.value_of("min_k")
        .unwrap()
        .trim()
        .parse()
        .expect("Need an integer!");
    let max_k: u32 = cli.value_of("max_k")
        .unwrap()
        .trim()
        .parse()
        .expect("Need an integer!");

    let fq_records = read_fq(seq_fp);
    let alphabet = alphabets::dna::iupac_alphabet();
    let rank = RankTransform::new(&alphabet);

    fq_records
        .map(|r| r.expect("Error reading FASTQ record"))
        .map(|r| {
                 println!("{}", UniqueKs::new(min_k, max_k, &rank)
                     .count(r.seq())
                     .tabular("\t"))
             }).collect::<Vec<()>>();
}

fn read_fq(path: &str) -> Box<Iterator<Item = Result<fastq::Record, std::io::Error>>> {
    match path {
        "-" => Box::new(fastq::Reader::new(io::stdin()).records()),
        _ => Box::new(fastq::Reader::from_file(path).expect("Could not open specified file").records())
    }
}

struct UniqueKs<'a> {
    min_k: u32,
    max_k: u32,
    rank: &'a RankTransform,
    freqs: Vec<u32>,
}

impl<'a> UniqueKs<'a> {
    fn new(min_k: u32, max_k: u32, rank: &'a RankTransform) -> Self {
        UniqueKs {
            min_k,
            max_k,
            rank,
            freqs: Vec::new(),
        }
    }

    fn count(&self, text: &[u8]) -> Self {
        let mut freqs: Vec<u32> = Vec::new();
        for k in self.min_k..self.max_k {
            freqs.push(unique_qgrams(text, k, self.rank) as u32);
        }
        UniqueKs { freqs, ..*self }
    }

    fn tabular(&self, sep: &str) -> String {
        format!("{}", self.freqs.iter().format(sep))
    }
}

impl<'a> Display for UniqueKs<'a> {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        write!(f, "{}:{}, {:?}", self.min_k, self.max_k, self.freqs)
    }
}

fn unique_qgrams(text: &[u8], q: u32, rank: &RankTransform) -> usize {
    rank.qgrams(q, text)
        .collect::<FnvHashSet<usize>>()
        .len()
}


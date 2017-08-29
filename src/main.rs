extern crate bio;
extern crate clap;
extern crate fnv;

use std::vec::Vec;
use std::io;
use std::str::FromStr;

use bio::alphabets;
use bio::alphabets::RankTransform;
use bio::io::{fastq, fasta};

use fnv::FnvHashSet;

use clap::{App, Arg, ArgMatches};

fn main() {

    let args = App::new("Measure/filter sequence complexity based on unique k-mers")
            .arg(Arg::with_name("k")
                .short("k")
                .takes_value(true)
                .help("length of k-mer to use")
                .default_value("4"))
            .arg(Arg::with_name("fasta")
                .long("fasta")
                .short("f")
                .takes_value(false)
                .help("input is in fasta format"))
        .get_matches();

    calculate(&args);
}

fn calculate(args: &ArgMatches) {

    let k: u32 = args
        .value_of("k")
        .expect("missing k!")
        .trim()
        .parse()
        .expect("Need an integer!");

    let window_size = k*4;

    let freqs: Vec<KFrequencies>;
    let alphabet = alphabets::dna::iupac_alphabet();
    let rank = RankTransform::new(&alphabet);

    if args.is_present("fasta") {
        let fa_records = fasta::Reader::new(io::stdin()).records();
        freqs = fa_records
            .map(|r| r.expect("Error reading FASTQ record"))
            .map(|r| KFrequencies::new(k, &rank, r.seq(), r.id().unwrap()))
            .collect();
    } else {
        let fq_records = fastq::Reader::new(io::stdin()).records();
        freqs = fq_records
            .map(|r| r.expect("Error reading FASTQ record"))
            .map(|r| KFrequencies::new(k, &rank, r.seq(), r.id().unwrap()))
            .collect();
    }
    
    freqs
        .iter()
        .map(|kf| println!("{}\t{:?}\t{:?}", kf.id, kf.freq, kf.len))
        .collect::<Vec<()>>();

}

#[derive(Debug)]
struct KFrequencies {
    id: String,
    len: usize,
    freq: usize
}

#[derive(Debug)]
struct Interval {
    start: usize,
    end: usize
}

impl KFrequencies {
    fn new(k: u32, rank: &RankTransform, seq: &[u8], id: &str, ) -> Self {
        let len = seq.len();
        let freq = unique_qgrams(seq, k, rank);
        let intervals = find_lc_regions(seq, k, rank, 32);
        println!("{:?}", intervals);
        // let lc_regions = collapse_regions(starts, ends);
        // println!("{:?}", lc_regions);
        KFrequencies {
            id: String::from_str(id).unwrap(),
            len,
            freq
        }
    }
}

fn unique_qgrams(text: &[u8], q: u32, rank: &RankTransform) -> usize {
    rank.qgrams(q, text)
        .collect::<FnvHashSet<usize>>()
        .len()
}

fn find_lc_regions(text: &[u8], q: u32, rank: &RankTransform, window_size: usize) -> Vec<Interval> {
    // There are l-k+1 k-grams in a sequence of length l
    let qgrams = rank.qgrams(q, text).collect::<Vec<usize>>();
    let mut intervals: Vec<Interval> = Vec::new();
    for (idx, window) in qgrams.as_slice().windows(window_size).enumerate() {
        let n_unique = window.into_iter().collect::<FnvHashSet<&usize>>().len();
        let window_complexity = n_unique as f64 / window_size as f64;
        if window_complexity < 0.55 {
            let start = idx;
            let end = idx + (window_size - 1 + q as usize);
            intervals.push(Interval{ start, end });
        }
    }
    return intervals;
}

fn collapse_intervals(intervals: Vec<Interval>) -> Vec<Interval> {
    let mut collapsed: Vec<Interval> = Vec::new();
    let mut last = 0;
    let mut len = 0;
    for idx in starts {
        if idx == last+1 {
            len += 1;
        } else {
            if len > 0 {
                collapsed.push(len);    
            }        
            len = 0;
        }
        last = idx;
    }
    if len > 0 {
        collapsed.push(len);
    }
    return collapsed;
}
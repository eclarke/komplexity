extern crate bio;
extern crate clap;
extern crate fnv;

use std::vec::Vec;
use std::io;
use std::iter;
use std::collections::VecDeque;

use bio::alphabets;
use bio::alphabets::RankTransform;
use bio::io::{fastq, fasta};

use fnv::{FnvHashSet, FnvHashMap};

use clap::{App, Arg};

fn main() {

    let args = App::new("kz")
        .about("Calculate the overall complexity of a sequence. Reads from stdin and writes to stdout.")
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
        .arg(Arg::with_name("window_size")
            .long("window_size")
            .short("w")
            .takes_value(true)
            .default_value("32")
            .help("window size for masking"))
        .arg(Arg::with_name("threshold")
            .long("threshold")
            .short("t")
            .takes_value(true)
            .default_value("0.55")
            .help("complexity threshold ([0:1], 0 = most stringent, 1 = least)"))
        .arg(Arg::with_name("mask")
            .long("mask")
            .short("m")
            .takes_value(false)
            .help("use sliding window to mask low-complexity regions"))
        .arg(Arg::with_name("lower_case")
            .long("lower_case")
            .short("l")
            .takes_value(false)
            .help("mask using lower-case symbols rather than Ns"))
    .get_matches();

    let record_type = match args.is_present("fasta") {
        true => RecordType::Fasta,
        false => RecordType::Fastq,
    };

    let task = match args.is_present("mask") {
        true => Task::Mask,
        false => Task::Measure
    };

    let mask_type = match args.is_present("lower_case") {
        true => MaskType::LowerCase,
        false => MaskType::N
    };

    let k: u32 = args
        .value_of("k")
        .unwrap()
        .trim()
        .parse()
        .expect("k must be an integer");

    if k > 12 {
        // Because we use the extended IUPAC alphabet, we're restricted to
        // smaller ks (though this doesn't matter for our purposes)
        error_exit("-k must be less than or equal to 12")
    }

    let threshold: f64 = args
        .value_of("threshold")
        .unwrap()
        .trim()
        .parse()
        .expect("'--threshold' must be a number between 0-1");
    
    if threshold < 0.0 || threshold > 1.0 {
        error_exit("'--threshold' must be a number between 0-1");
    }
    
    let window_size: usize = args
        .value_of("window_size")
        .unwrap()
        .trim()
        .parse()
        .expect("'--window_size' must be an integer greater than 0");

    complexity(record_type, task, k, threshold, window_size, mask_type);
}

#[derive(Debug)]
struct Interval {
    start: usize,
    end: usize
}

enum RecordType {
    Fasta,
    Fastq
}

enum Task {
    Mask,
    Measure
}

enum MaskType {
    N,
    LowerCase
}

fn complexity(record_type: RecordType, task: Task, k: u32, threshold: f64, window_size: usize, mask_type: MaskType) {
    let alphabet = alphabets::dna::iupac_alphabet();
    let rank = RankTransform::new(&alphabet);
    
    match record_type {
        RecordType::Fasta => {
            let mut writer = fasta::Writer::new(io::stdout());
            let records = fasta::Reader::new(io::stdin()).records();
            records
                .map(|r| r.expect("Error reading FASTA record"))
                .map(|r| {
                    let id = r.id();
                    let seq = r.seq();
                    match task {
                        Task::Mask => {
                            let seq = mask_sequence(seq, &rank, k, threshold, window_size, &mask_type);
                            writer.write(id, r.desc(), &seq).unwrap();
                        },
                        Task::Measure => {
                            let length = seq.len();
                            let kmers = unique_kmers(seq, k, &rank);
                            println!("{}\t{}\t{}\t{}", id, length, kmers, kmers as f64 / length as f64);
                        }
                    } 
                    
                })
                .collect::<Vec<()>>();
        }, 
        RecordType::Fastq => {
            let mut writer = fastq::Writer::new(io::stdout());
            let records = fastq::Reader::new(io::stdin()).records();
            records
                .map(|r| r.expect("Error reading FASTQ record"))
                .map(|r| {
                    let id = r.id();
                    let seq = r.seq();
                    match task {
                        Task::Mask => {
                            let seq = mask_sequence(seq, &rank, k, threshold, window_size, &mask_type);
                            writer.write(id, r.desc(), &seq, r.qual()).unwrap();
                        },
                        Task::Measure => {
                            let length = seq.len();
                            let kmers = unique_kmers(seq, k, &rank);
                            println!("{}\t{}\t{}\t{}", id, length, kmers, kmers as f64 / length as f64);
                        }
                    } 
                })
                .collect::<Vec<()>>();
        }
    }
}

fn mask_sequence(seq: &[u8], rank: &RankTransform, k: u32, threshold: f64, window_size: usize, mask_type: &MaskType) -> Vec<u8> {
    // let intervals = lc_intervals(seq, k, rank, threshold, window_size);
    let intervals = lc2(seq, k, rank, threshold, window_size);
    mask_intervals(seq, intervals, mask_type)
}

fn unique_kmers(text: &[u8], k: u32, rank: &RankTransform) -> usize {
    rank.qgrams(k, text)
        .collect::<FnvHashSet<usize>>()
        .len()
}

fn lc2(text: &[u8], q: u32, rank: &RankTransform, threshold: f64, window_size: usize) -> Vec<Interval> {
    // Bounds checking
    let q = q as usize;

    let mut intervals: Vec<Interval> = Vec::new();
    let mut window: VecDeque<usize> = VecDeque::with_capacity(window_size);
    let mut kmer_iterator = rank.qgrams(q as u32, text).into_iter();
    let mut kmers: FnvHashMap<usize, usize> = FnvHashMap::default();
 
    // Init: fill window buffer
    for _ in 0..window_size {
        match kmer_iterator.next() {
            Some(kmer) => window.push_back(kmer),
            None => break
        }
    }
    // Count kmers in window
    for kmer in window.iter() {
        let n = kmers.entry(*kmer).or_insert(0);
        *n += 1;
    }

    let mut idx = 0;
    loop {
        let window_complexity = kmers.len() as f64 / window.len() as f64;
        if window_complexity < threshold {
            let start = idx;
            let end = idx + (window.len() - 1 + q as usize);
            intervals.push(Interval{ start, end });
        }
        match kmer_iterator.next() {
            Some(kmer) => {
                let prev = window.pop_front().unwrap();
                window.push_back(kmer);
                // Update kmer counts: remove 1 from leaving kmer
                let prev_n = *kmers.get(&prev).unwrap();
                if prev_n == 1 {
                    kmers.remove(&prev);
                } else {
                    kmers.insert(prev, prev_n-1);
                }
                // Update kmer counts: add 1 for entering kmer
                let next_n = kmers.entry(kmer).or_insert(0);
                *next_n += 1;
                // Update index
                idx += 1;
            },
            None => break,
        }
    }
    return collapse_intervals(intervals);
}

fn collapse_intervals(intervals: Vec<Interval>) -> Vec<Interval> {
    let mut collapsed: Vec<Interval> = Vec::new();
    let mut intervals = intervals.into_iter();
    if let Some(mut current) = intervals.next() {
        for interval in intervals {
            if interval.start < current.end {
                current.end = interval.end;
            } else {
                collapsed.push(current);
                current = interval;
            }
        }
        collapsed.push(current);
    }
    return collapsed;
}

fn mask_intervals(seq: &[u8], intervals: Vec<Interval>, mask_type: &MaskType) -> Vec<u8> {
    let mut new_seq: Vec<u8> = Vec::with_capacity(seq.len());
    let mut last = Interval{start: 0, end: 0};
    for interval in intervals {
        let intervening = &seq[last.end .. interval.start];
        new_seq.extend_from_slice(intervening);
        match *mask_type {
            MaskType::LowerCase => new_seq.append(&mut lowercase(&seq[interval.start..interval.end])),
            MaskType::N => new_seq.extend(iter::repeat(b'N').take(interval.end-interval.start)),
        };
        // for _ in interval.start..interval.end {
        //     new_seq.push(b'N');
        // }
        last = interval;
    }
    let end = &seq[last.end..];
    new_seq.extend_from_slice(end);
    return new_seq;
}

fn error_exit(msg: &str) {
    eprintln!("{}", msg);
    std::process::exit(1);
}

fn lowercase(seq: &[u8]) -> Vec<u8> {
    // ACGTRYSWKMBDHVNZ
    let mut new = Vec::with_capacity(seq.len());
    for b in seq {
        let b = match *b {
            b'A' => b'a',
            b'C' => b'c',
            b'G' => b'g',
            b'T' => b't',
            b'R' => b'r',
            b'Y' => b'y',
            b'S' => b's',
            b'W' => b'w',
            b'K' => b'k',
            b'M' => b'm',
            b'B' => b'b',
            b'D' => b'd',
            b'H' => b'h',
            b'V' => b'v',
            b'N' => b'n',
            b'Z' => b'z',
            _ => *b,
        };
        new.push(b);
    }
    return new;
}
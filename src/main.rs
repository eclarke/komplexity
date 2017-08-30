extern crate bio;
extern crate clap;
extern crate fnv;

use std::vec::Vec;
use std::io;

use bio::alphabets;
use bio::alphabets::RankTransform;
use bio::io::{fastq, fasta};

use fnv::{FnvHashSet, FnvHashMap};

use clap::{App, Arg};

fn main() {

    let args = App::new("kz")
            .about("Calculate the overall complexity of a sequence")
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
                .help("use sliding window to mask low-complexity regions")
            )

        .get_matches();

    let record_type = match args.is_present("fasta") {
        true => RecordType::Fasta,
        false => RecordType::Fastq,
    };

    let task = match args.is_present("mask") {
        true => Task::Mask,
        false => Task::Measure
    };

    let k: u32 = args
        .value_of("k")
        .unwrap()
        .trim()
        .parse()
        .expect("Need an integer!");

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

    complexity(record_type, task, k, threshold, window_size);
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

fn complexity(record_type: RecordType, task: Task, k: u32, threshold: f64, window_size: usize) {
    let alphabet = alphabets::dna::iupac_alphabet();
    let rank = RankTransform::new(&alphabet);
    match record_type {
        RecordType::Fasta => {
            let mut writer = fasta::Writer::new(io::stdout());
            let records = fasta::Reader::new(io::stdin()).records();
            records
                .map(|r| r.expect("Error reading FASTA record"))
                .map(|r| {
                    let id = r.id().unwrap();
                    let seq = r.seq();
                    match task {
                        Task::Mask => {
                            let seq = mask_sequence(seq, &rank, k, threshold, window_size);
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
                    let id = r.id().unwrap();
                    let seq = r.seq();
                    match task {
                        Task::Mask => {
                            let seq = mask_sequence(seq, &rank, k, threshold, window_size);
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

fn mask_sequence(seq: &[u8], rank: &RankTransform, k: u32, threshold: f64, window_size: usize) -> Vec<u8> {
    let intervals = lc_intervals(seq, k, rank, threshold, window_size);
    mask_intervals(seq, intervals)
}

fn unique_kmers(text: &[u8], k: u32, rank: &RankTransform) -> usize {
    rank.qgrams(k, text)
        .collect::<FnvHashSet<usize>>()
        .len()
}

fn lc_intervals(text: &[u8], q: u32, rank: &RankTransform, threshold: f64, window_size: usize) -> Vec<Interval> {
    let qgrams: Vec<usize> = rank.qgrams(q, text).into_iter().collect();
    let mut intervals: Vec<Interval> = Vec::new();
    let mut kmers: FnvHashMap<usize, usize> = FnvHashMap::default();
    // We keep track of the kmer appearing at the beginning of the window (`prev`)
    // and subtract 1 from its count (or remove it if < 1) when we move on 
    // to the next window. We add 1 to the count of the newest kmer we see 
    // (the last in this window).
    let mut prev: usize = 0;
    for (idx, window) in qgrams.as_slice().windows(window_size).enumerate() {
        {
            if idx == 0 {
                prev = window[0];
                window.iter().map(|k| {
                    let n = kmers.entry(*k).or_insert(0);
                    *n += 1;
                }).collect::<Vec<()>>();
            } else {
                let n = *kmers.get(&prev).unwrap();
                if n == 1 {
                    kmers.remove(&prev);
                } else {
                    kmers.insert(prev, n-1);
                }
                let next = window.last().unwrap();
                let n = kmers.entry(*next).or_insert(0);
                *n += 1;
                prev = window[0];  
            }
        }
        let n_unique = kmers.len();
        let window_complexity = n_unique as f64 / window_size as f64;
        if window_complexity < threshold {
            let start = idx;
            let end = idx + (window_size - 1 + q as usize);
            intervals.push(Interval{ start, end });
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

fn mask_intervals(seq: &[u8], intervals: Vec<Interval>) -> Vec<u8> {
    let mut new_seq: Vec<u8> = Vec::with_capacity(seq.len());
    let mut last = Interval{start: 0, end: 0};
    for interval in intervals {
        let intervening = &seq[last.end .. interval.start];
        new_seq.extend_from_slice(intervening);
        for _ in interval.start..interval.end {
            new_seq.push(b'N');
        }
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
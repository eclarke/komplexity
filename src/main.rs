extern crate bio;
extern crate clap;
extern crate fnv;

use std::vec::Vec;
use std::io::{self, BufRead};
use std::ops::Div;
use std::cmp::Ordering;

use bio::alphabets;
use bio::alphabets::RankTransform;
use bio::io::fastq;

use clap::{App, Arg, ArgMatches};

use fnv::FnvHashSet;

fn main() {

    let args = App::new("Measure/filter sequence complexity based on unique k-mers")
        .subcommand(App::new("calculate")
            .arg(Arg::with_name("k")
                .short("k")
                .takes_value(true)
                .help("length of k-mer to use")
                .default_value("4")))
        .subcommand(App::new("filter")
            .arg(Arg::with_name("threshold")
                .short("t")
                .long("threshold")
                .takes_value(true)
                .help("Minimum z-score to keep")
                .default_value("-1.5"))
            .arg(Arg::with_name("invert")
                .long("invert")
                .takes_value(false)
                .help("Keep seqs lower than threshold (for debugging)")))
        .get_matches();

    match args.subcommand() {
        ("calculate", Some(subargs)) => calculate(subargs),
        ("filter", Some(subargs)) => filter(subargs),
        (_, _) => {println!("{}",args.usage()); std::process::exit(1);}
    }
}

fn calculate(args: &ArgMatches) {

    let k: u32 = args
        .value_of("k")
        .expect("missing k!")
        .trim()
        .parse()
        .expect("Need an integer!");

    let fq_records = fastq::Reader::new(io::stdin()).records();
    let alphabet = alphabets::dna::iupac_alphabet();
    let rank = RankTransform::new(&alphabet);

    let freqs: Vec<KFrequencies> = fq_records
        .map(|r| r.expect("Error reading FASTQ record"))
        .map(|r| KFrequencies::new(k, &rank, r))
        .collect();

    freqs
        .iter()
        .map(|kf| println!("{}\t{:.4}", kf.id, kf.score))
        .collect::<Vec<()>>();

}

fn filter(args: &ArgMatches) {
    let threshold: f64 = args
        .value_of("threshold")
        .expect("missing threshold!")
        .trim()
        .parse()
        .expect("Need a numeric value!");

    let invert: bool = args.is_present("invert");

    let mut input: Vec<Score> = Vec::new();
    let stdin = io::stdin();

    for line in stdin.lock().lines() {
        match line {
            Ok(s) => input.push(Score::new(s)),
            Err(_) => break,
        }
    };

    let comp = if invert { Ordering::Less } else { Ordering::Greater};

    let mean: f64 = input.iter()
        .map(|s| s.score)
        .sum::<f64>()
        .div((input.len() as f64));
    let sd: f64 = input.iter()
        .map(|s| (s.score - mean).powi(2))
        .sum::<f64>()
        .div((input.len() - 1) as f64)
        .sqrt();
    input
        .iter()
        .filter(|s| ((s.score - mean)/sd).partial_cmp(&threshold) == Some(comp))
        .map(|s| println!("{}", s.id))
        .collect::<Vec<()>>();
}

struct Score {
    id: String,
    score: f64,
}

impl Score {
    fn new(s: String) -> Self {
        let parts: Vec<&str> = s.splitn(2, "\t").collect();
        let id: String = parts[0].trim().to_owned();
        let score: f64 = parts[1].trim().parse().expect(&format!("Error reading score from {}", s));
        return Score{ id, score }
    }
}

#[derive(Debug)]
struct KFrequencies {
    id: String,
    k: u32,
    freq: usize,
    len: usize,
    score: f64,
}

impl KFrequencies {
    fn new(k: u32, rank: &RankTransform, r: fastq::Record) -> Self {
        let freq = unique_qgrams(r.seq(), k, rank);
        let len = r.seq().len();
        let score = (freq as f64) / (len as f64);
        KFrequencies {
            id: r.id().unwrap_or_default().to_owned(),
            k,
            freq,
            len,
            score,
        }
    }
}

fn unique_qgrams(text: &[u8], q: u32, rank: &RankTransform) -> usize {
    rank.qgrams(q, text)
        .collect::<FnvHashSet<usize>>()
        .len()
}


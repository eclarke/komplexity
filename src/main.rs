extern crate bio;
extern crate clap;
extern crate fnv;
extern crate indicatif;

use std::thread;
use std::io::{Read, BufRead, BufReader};
use std::fs::File;
use std::sync::mpsc;
use std::sync::Arc;
use std::sync::Mutex;
use std::vec::Vec;

use bio::alphabets;
use bio::alphabets::RankTransform;
use bio::io::fastq::{self, Record};

use clap::{Arg, App};

use indicatif::ProgressBar;

use fnv::FnvHashSet;

const LF: u8 = '\n' as u8;

fn main() {

    let cli = App::new("Unique Kmers")
        .arg(Arg::with_name("seq")
                 .short("s")
                 .long("seq")
                 .required(true)
                 .takes_value(true)
                 .help("Path to FASTQ file"))
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
        .arg(Arg::with_name("threads")
                 .short("t")
                 .long("threads")
                 .default_value("1")
                 .takes_value(true)
                 .help("Threads to use"))
        .get_matches();

    let seq_fp = cli.value_of("seq").unwrap();
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
    let threads: u32 = cli.value_of("threads")
        .unwrap()
        .trim()
        .parse()
        .expect("Need an integer!");

    println!("Reading FASTQ file...");
    let records = fastq::Reader::from_file(seq_fp)
        .expect("Invalid FASTQ file!")
        .records();
    println!("Collecting lines...");

    let n_records = lines(seq_fp).expect("Error reading file.") / 4;

    #[derive(Debug)]
    struct UniqueKs {
        ks: Vec<u32>,
    };

    impl UniqueKs {
        fn new(min_k: u32, max_k: u32, text: &[u8], rank: &RankTransform) -> Self {
            let mut ks = Vec::new();
            for k in min_k..max_k {
                ks.push((unique_qgrams(text, k, rank) as u32));
            }
            return UniqueKs { ks };

        }
    }

    let alphabet = alphabets::dna::iupac_alphabet();
    let rank = Arc::new(RankTransform::new(&alphabet));

    let (tx, rx) = mpsc::channel();
    let mut freqs: Vec<UniqueKs> = Vec::new();
    let bar = ProgressBar::new((n_records as u64));
    let chunksize = ((n_records as f64) / (threads as f64)).ceil() as usize;

    println!("Reading all records into memory...");
    let full_records: Vec<Record> = records
        .map(|r| {
                 bar.inc(1);
                 r.expect("Error reading record")
             })
        .collect();
    bar.finish();
    let record_slice = &full_records[..];
    println!("Processing records...");

    let multibar = indicatif::MultiProgress::new();

    for chunk in record_slice.chunks(chunksize) {
        let thread_bar = multibar.add(ProgressBar::new(chunksize as u64));
        let thread_tx = tx.clone();
        let thread_rank = rank.clone();
        let thread_chunk = chunk.to_vec();
        // let record_chunk = &full_records[i as usize..(i as usize)+chunksize];
        thread::spawn(move || { 
            for record in thread_chunk {
                thread_bar.inc(1);
                thread_tx
                    .send(UniqueKs::new(min_k, max_k, record.seq(), &thread_rank))
                    .unwrap();     
            }
            thread_bar.finish();
        });
    }
    multibar.join_and_clear().unwrap();
    drop(tx);

    loop {
        match rx.recv() {
            Ok(r) => freqs.push(r),
            Err(_) => break,
        };
    }

    println!("{:?}", freqs.len());
    println!("{:?}", freqs[1]);
    println!("{:?}", freqs[254]);
}



fn unique_qgrams(text: &[u8], q: u32, rank: &RankTransform) -> usize {
    rank.qgrams(q, text)
        .collect::<FnvHashSet<usize>>()
        .len()
}

fn lines<P: AsRef<std::path::Path>>(f: P) -> Result<usize, std::io::Error> {
    let mut reader = BufReader::new(try!(File::open(f)));
    let mut line_count = 0;
    let mut raw_line: Vec<u8> = Vec::new();
    loop {
        match reader.read_until(LF, &mut raw_line) {
            Ok(n) if n == 0 => break,
            Err(e) => return Err(e),
            Ok(_) => {}
        };
        if *raw_line.last().unwrap() == LF {
            line_count += 1;
        };
    }
    return Ok(line_count);
}


// fn naive_qgrams(text: &[u8], q: u32) -> usize {
//     let mut set = FnvHashSet::default();
//     println!("{}", text.len());
//     for idx in 0..text.len() {
//         let end = if idx + (q as usize) > text.len() {
//             break;
//         } else {
//             idx + (q as usize)
//         };
//         set.insert(&text[idx..end]);
//         if idx % 1000 == 0 {
//             println!("{}", (idx as f64) / (text.len() as f64));
//         }
//     }
//     set.len()

// }


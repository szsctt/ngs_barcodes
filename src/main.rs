use std::collections::HashMap;
use std::fs::File;

use flate2::read::GzDecoder;
use bio::io::fastq;
use serde::Deserialize;
use std::io::BufReader;

#[derive(Debug, Deserialize, PartialEq)]
struct Barcode {
    #[serde(rename = "type")]
    barcode_type: String,
    start: usize,
    mismatches: u8,
    barcodes: HashMap<String, String>,
}

#[derive(Debug, Deserialize, PartialEq)]
struct BarcodeSet {
    #[serde(flatten)]
    set: HashMap<String, Barcode>,
}

fn main() {
    // Open and parse the YAML file containing the barcodes
    let barcode_file = File::open("barcodes.yaml").unwrap();
    let barcode_sets: Vec<BarcodeSet> = serde_yaml::from_reader(barcode_file).unwrap();

    // Initialize a HashMap to keep track of the counts of each barcode
    let mut counts: HashMap<String, u32> = HashMap::new();

    // Open the gzipped FASTQ file and create a new FASTQ reader
    let f = File::open("reads.fastq.gz").unwrap();
    let d = GzDecoder::new(f);
    let reader = BufReader::new(d);
    let fastq_reader = fastq::Reader::new(reader);

    // Iterate over the FASTQ records
    for result in fastq_reader.records() {
        let record = result.unwrap();
        let sequence = record.seq();

        // Initialize a vector to keep track of the barcodes found in each record
        let mut barcodes_found = Vec::new();

        // For each record, iterate over all barcode sets and barcodes
        // in each set, and check if the barcode is found at the expected
        // position in the sequence
        for barcode_set in &barcode_sets {
            for (set_name, set) in &barcode_set.set {
                for (bc_name, bc_seq) in &set.barcodes {
                    if sequence[set.start..set.start + bc_seq.len()] == *bc_seq.as_bytes() {
                        barcodes_found.push(format!("{}__{}", set_name, bc_name));
                    }
                }
            }
        }

        // If any barcodes were found in the record, increment the count for the
        // combination of barcodes found
        if !barcodes_found.is_empty() {
            barcodes_found.sort();
            *counts.entry(barcodes_found.join("__")).or_insert(0) += 1;
        }
    }

    // Print the counts of each barcode
    for (k, v) in &counts {
        println!("{} {}", k, v);
    }
}

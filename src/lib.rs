//! # Utilrs
//!
//! `utilrs` is a collection of utility functions used in GRAMEP.
//!
#![allow(unused_imports)]
#[macro_use(concat_string)]
extern crate concat_string;
use anyhow::{Error, Result};
use bio::alignment::distance::simd::*;
use bio::io::fasta;
use itertools::Itertools;
use pyo3::prelude::*;
use rayon::prelude::*;
use regex::Regex;
use rustc_hash::FxHashMap;
use spinners::{Spinner, Spinners};
use std::collections::HashSet;
use std::sync::Arc;

enum Value {
    Usize(u32),
    Str(String),
}

impl ToPyObject for Value {
    fn to_object(&self, py: Python) -> Py<PyAny> {
        match self {
            Value::Usize(u) => u.to_object(py),
            Value::Str(s) => s.to_object(py),
        }
    }
}

/// Checks if a sequence contains any characters that are not part of the specified dictionary.
///
/// # Arguments
///
/// * `seq` - A string slice representing the sequence to be checked.
/// * `dict` - A string slice representing the dictionary containing the allowed characters.
///
/// # Returns
///
/// * `true` if the sequence contains forbidden characters (i.e., characters not in the dictionary).
/// * `false` if all characters in the sequence are allowed.
///
/// # Examples
///
/// ```
/// use utilrs::contains_forbidden_chars;
/// // DNA sequence check (no forbidden characters, should return false)
/// assert_eq!(contains_forbidden_chars("ACGT", "DNA"), false);
///
/// // DNA sequence with a forbidden character (should return true)
/// assert_eq!(contains_forbidden_chars("ACGTX", "DNA"), true);
///
/// // RNA sequence check (no forbidden characters, should return false)
/// assert_eq!(contains_forbidden_chars("ACGU", "RNA"), false);
///
/// // RNA sequence with a forbidden character (should return true)
/// assert_eq!(contains_forbidden_chars("ACGUX", "RNA"), true);
///
/// // ALL dictionary check (only allowing NNNNN)
/// assert_eq!(contains_forbidden_chars("ANNNCTGTG", "ALL"), true);
/// ```
///
fn contains_forbidden_chars(seq: &str, dict: &str) -> bool {
    let alphabet = &get_alphabet(dict);
    let re = regex::Regex::new(alphabet).unwrap();
    re.is_match(seq)
}

/// Returns a regular expression pattern representing forbidden characters for a given dictionary.
///
/// # Arguments
///
/// * `dict` - A string slice that specifies the type of dictionary:
///   - `"DNA"`: Returns a pattern for forbidden characters in DNA sequences.
///   - `"RNA"`: Returns a pattern for forbidden characters in RNA sequences.
///   - Other values: Returns a pattern for forbidden characters in general sequences.
///
/// # Returns
///
/// * A `String` representing the regular expression pattern of forbidden characters.
///
/// # Examples
///
/// ```
/// use utilrs::get_alphabet;
/// // Testing the function with different dictionaries
/// assert_eq!(get_alphabet("DNA"), "B|D|E|F|H|I|J|K|L|M|N|O|P|Q|R|S|U|V|W|X|Y|Z");
/// assert_eq!(get_alphabet("RNA"), "B|D|E|F|H|I|J|K|L|M|N|O|P|Q|R|S|T|V|W|X|Y|Z");
/// assert_eq!(get_alphabet("General"), "B|D|E|F|H|I|J|K|L|M|O|P|Q|R|S|V|W|X|Y|Z");
/// assert_eq!(get_alphabet("Unknown"), "B|D|E|F|H|I|J|K|L|M|O|P|Q|R|S|V|W|X|Y|Z");
/// ```
///
fn get_alphabet(dict: &str) -> String {
    let alphabet = String::from(if dict == "DNA" {
        "B|D|E|F|H|I|J|K|L|M|N|O|P|Q|R|S|U|V|W|X|Y|Z"
    } else if dict == "RNA" {
        "B|D|E|F|H|I|J|K|L|M|N|O|P|Q|R|S|T|V|W|X|Y|Z"
    } else {
        "B|D|E|F|H|I|J|K|L|M|O|P|Q|R|S|V|W|X|Y|Z"
    });
    alphabet
}

/// Splits a sequence into overlapping substrings (k-mers) with a given step size.
///
/// # Arguments
///
/// * `seq` - A string slice representing the sequence to be split.
/// * `k` - The length of each substring (k-mer).
/// * `step` - The step size to move the starting index for each new substring.
///
/// # Returns
///
/// * A `Vec<String>` containing the k-mers extracted from the sequence.
///
/// # Examples
///
/// ```
/// use utilrs::split_seq;
/// // Example usage
/// let seq = "ABCDEFGH";
/// let k = 3;
/// let step = 1;
/// let kmers = split_seq(seq, k, step);
/// assert_eq!(kmers, vec!["ABC", "BCD", "CDE", "DEF", "EFG", "FGH"]);
///
/// // Example with a larger step
/// let step = 3;
/// let kmers = split_seq(seq, k, step);
/// assert_eq!(kmers, vec!["ABC", "DEF", "GHI"]);
///
/// // Example where step is larger than sequence length
/// let step = 10;
/// let kmers = split_seq(seq, k, step);
/// assert_eq!(kmers, vec!["ABC"]);
///
/// // Example where k is larger than sequence length
/// let k = 10;
/// let step = 1;
/// let kmers = split_seq(seq, k, step);
/// assert_eq!(kmers, vec![]);
/// ```
///
fn split_seq(seq: &str, k: usize, step: usize) -> Vec<String> {
    let end = seq.chars().count() - k + step;
    let mut index = 0;
    let mut kmers = Vec::new();
    while index + step < end {
        kmers.push(seq[index..index + k].to_owned());
        index += step;
    }
    kmers
}

/// Finds the positions and types of mutations (single nucleotide polymorphisms) between a pattern and a reference sequence.
///
/// # Arguments
///
/// * `pattern` - A slice of bytes representing the k-mer pattern to compare against the reference sequence.
/// * `ref_bytes` - A slice of bytes representing the reference sequence.
/// * `max_dist` - The maximum number of allowed mutations to find.
///
/// # Returns
///
/// * A `Vec<(String, usize)>` containing tuples where:
///   - The `String` represents the type of mutation (SNP) as a concatenation of the reference and pattern characters.
///   - The `usize` is the position of the mutation in the reference sequence.
///
/// # Examples
///
/// ```
/// use utilrs::get_kmer_mutation_index;
/// // Example usage
/// let pattern = b"ACGT";
/// let ref_bytes = b"AGCT";
/// let max_dist = 2;
/// let mutations = get_kmer_mutation_index(pattern, ref_bytes, max_dist);
/// assert_eq!(mutations, vec![("A".to_string(), 0), ("C".to_string(), 2)]);
///
/// // Example with more allowed mutations
/// let max_dist = 10;
/// let mutations = get_kmer_mutation_index(pattern, ref_bytes, max_dist);
/// assert_eq!(mutations, vec![("A".to_string(), 0), ("C".to_string(), 2)]);
///
/// // Example with no mutations
/// let pattern = b"ACGT";
/// let ref_bytes = b"ACGT";
/// let max_dist = 1;
/// let mutations = get_kmer_mutation_index(pattern, ref_bytes, max_dist);
/// assert_eq!(mutations, vec![]);
///
/// // Example with empty pattern
/// let pattern = b"";
/// let ref_bytes = b"ACGT";
/// let max_dist = 1;
/// let mutations = get_kmer_mutation_index(pattern, ref_bytes, max_dist);
/// assert_eq!(mutations, vec![]);
///
/// // Example with empty reference
/// let pattern = b"ACGT";
/// let ref_bytes = b"";
/// let max_dist = 1;
/// let mutations = get_kmer_mutation_index(pattern, ref_bytes, max_dist);
/// assert_eq!(mutations, vec![]);
/// ```
///
fn get_kmer_mutation_index(
    pattern: &[u8],
    ref_bytes: &[u8],
    max_dist: usize,
) -> Vec<(String, usize)> {
    let mut i = 0;
    let mut counter = 1;
    let mut mutations = Vec::new();

    while i < ref_bytes.len() && counter <= max_dist {
        if ref_bytes[i] != pattern[i] {
            let snp = concat_string!(
                String::from(ref_bytes[i] as char),
                String::from(pattern[i] as char)
            );
            mutations.push((snp, i));
            counter += 1;
        }
        i += 1;
    }
    mutations
}

/// Processes a FASTA file to count k-mer frequencies and returns a vector of frequency maps.
///
/// This function reads a FASTA file, processes sequences in batches, and counts the occurrences of k-mers of a specified length (`k`).
/// It uses multi-threading to handle large files efficiently. The result is a vector of `FxHashMap`s, where each map contains k-mer counts
/// and additional information based on the provided arguments.
///
/// # Arguments
///
/// * `seq_path` - The path to the FASTA file containing the sequences to process.
/// * `k` - The length of k-mers to count.
/// * `step` - The step size for extracting overlapping k-mers.
/// * `dict` - A string specifying the dictionary type ('DNA', 'RNA', or other) used to determine forbidden characters.
/// * `variants_kmers` - A vector of k-mers to track specifically in the output.
/// * `predict_data` - A boolean indicating whether to include sequence IDs (`true`) or a single class name (`false`) in the output.
/// * `batch_size` - The number of sequences to process per batch.
///
/// # Returns
///
/// A `Vec<FxHashMap<String, Value>>` where each `FxHashMap` contains:
/// - k-mer frequencies.
/// - An "ID" or "CLASS" entry depending on the `predict_data` flag.
///
fn count_freq(
    seq_path: String,
    k: usize,
    step: usize,
    dict: String,
    variants_kmers: Vec<String>,
    predict_data: bool,
    batch_size: usize,
) -> Vec<FxHashMap<String, Value>> {
    let class_name = seq_path.split(".fasta").collect::<Vec<&str>>()[0];
    let class_name = class_name.split("/").collect::<Vec<&str>>();
    let class_name = class_name.last().unwrap().to_string();

    let reader = fasta::Reader::from_file(seq_path).unwrap();
    let pool: rayon::ThreadPool = rayon::ThreadPoolBuilder::new().build().unwrap();
    let (tx, rx) = std::sync::mpsc::channel();

    // Process sequences in batches
    let mut batch = Vec::new();
    for result in reader.records() {
        let result_data = result.unwrap();
        let seq: &[u8] = result_data.seq();
        let seq_name = result_data.id().to_string();
        let seq_str = String::from_utf8(seq.to_owned()).unwrap();
        let seq_str = seq_str.to_uppercase();
        if !contains_forbidden_chars(&seq_str, &dict) {
            if predict_data {
                batch.push((seq_str, seq_name.clone()));
            } else {
                batch.push((seq_str, class_name.clone()));
            }
            if batch.len() >= batch_size {
                let variants_kmers = variants_kmers.clone();
                let tx = tx.clone();
                let batch_clone = batch.clone();
                batch.clear();
                pool.spawn(move || {
                    // let mut batch_hm: FxHashMap<String, u32> = FxHashMap::default();
                    for (seq, seq_class) in batch_clone {
                        let mut hm: FxHashMap<String, u32> = FxHashMap::default();
                        let mut kmer_hm: FxHashMap<String, Value> = FxHashMap::default();
                        let end = seq.chars().count() - k + step;
                        let mut index = 0;
                        while index + step < end {
                            *hm.entry(seq[index..index + k].to_owned()).or_insert(0) += 1;
                            index += step;
                        }
                        for kmer in &variants_kmers {
                            *kmer_hm.entry(kmer.to_string()).or_insert(Value::Usize(0)) =
                                Value::Usize(*hm.get(kmer).unwrap_or(&0));
                        }
                        if predict_data {
                            kmer_hm.insert("ID".to_string(), Value::Str(seq_class.clone()));
                        } else {
                            kmer_hm.insert("CLASS".to_string(), Value::Str(seq_class.clone()));
                        }
                        tx.send(kmer_hm).unwrap();
                    }
                });
            }
        }
    }
    // Process any remaining sequences in the last batch
    if !batch.is_empty() {
        let tx = tx.clone();
        let variants_kmers = variants_kmers.clone();
        pool.spawn(move || {
            // let mut batch_hm: FxHashMap<String, u32> = FxHashMap::default();
            // let variants_kmers = variants_kmers.clone();
            for (seq, seq_class) in batch {
                let mut hm: FxHashMap<String, u32> = FxHashMap::default();
                let mut kmer_hm: FxHashMap<String, Value> = FxHashMap::default();
                let end = seq.chars().count() - k + step;
                let mut index = 0;
                while index + step < end {
                    *hm.entry(seq[index..index + k].to_owned()).or_insert(0) += 1;
                    index += step;
                }
                for kmer in &variants_kmers {
                    *kmer_hm.entry(kmer.to_string()).or_insert(Value::Usize(0)) =
                        Value::Usize(*hm.get(kmer).unwrap_or(&0));
                }
                if predict_data {
                    kmer_hm.insert("ID".to_string(), Value::Str(seq_class.clone()));
                } else {
                    kmer_hm.insert("CLASS".to_string(), Value::Str(seq_class.clone()));
                }
                tx.send(kmer_hm).unwrap();
            }
        });
    }

    drop(tx); // Close all senders

    let samples = rx.iter().collect::<Vec<FxHashMap<String, Value>>>();

    samples
}

/// Computes the intersection of k-mer sets.
///
/// This function takes a subset of k-mer variants and a mapping from k-mer
/// identifiers to sets of k-mers. It finds the common k-mers across all
/// sets specified by the subset.
///
/// # Arguments
///
/// * `subset` - A vector of k-mer variant identifiers. Each identifier is used
///   to look up a set of k-mers in `variants_exclusive_kmers`.
/// * `variants_exclusive_kmers` - A map where the key is a k-mer variant
///   identifier and the value is a vector of k-mers associated with that variant.
///
/// # Returns
///
/// A vector of k-mers that are present in all sets corresponding to the
/// identifiers in `subset`.
///
/// # Examples
///
/// ```
/// use rustc_hash::FxHashMap;
/// use utilrs::get_set_intersection;
///
/// let mut variants_exclusive_kmers: FxHashMap<String, Vec<String>> = FxHashMap::default();
/// variants_exclusive_kmers.insert("A".to_string(), vec!["ACGT".to_string(), "ACGG".to_string()]);
/// variants_exclusive_kmers.insert("B".to_string(), vec!["ACGT".to_string(), "ACGA".to_string()]);
/// variants_exclusive_kmers.insert("C".to_string(), vec!["ACGT".to_string(), "ACGC".to_string()]);
///
/// let subset = vec!["A".to_string(), "B".to_string()];
/// let intersection = get_set_intersection(subset, variants_exclusive_kmers);
/// assert_eq!(intersection, vec!["ACGT".to_string()]);
/// ```
///
fn get_set_intersection(
    subset: Vec<String>,
    variants_exclusive_kmers: FxHashMap<String, Vec<String>>,
) -> Vec<String> {
    let mut select_kmers: Vec<&Vec<String>> = Vec::new();

    for variant in subset {
        if let Some(kmers) = variants_exclusive_kmers.get(&variant) {
            select_kmers.push(kmers);
        }
    }

    if select_kmers.is_empty() {
        return Vec::new();
    }

    let mut result: HashSet<&String> = select_kmers[0].iter().collect();

    for temp_var in select_kmers.iter().skip(1) {
        let unique: HashSet<&String> = temp_var.into_iter().collect();

        result = result.intersection(&unique).cloned().collect();
    }
    result.into_iter().cloned().collect()
}

/// Computes the Shannon entropy
///
/// This function calculates the Shannon entropy of a probability distribution.
///
/// # Arguments
///
/// * `probs` - A slice of probabilities.
///
/// # Returns
///
/// The Shannon entropy of the probability distribution.
///
/// # Examples
///
/// ```
/// use utilrs::entropy;
///
/// let probs = vec![0.25, 0.25, 0.25, 0.25];
/// let entropy = entropy(&probs);
/// assert_eq!(entropy, 2.0);
/// ```
///
fn entropy(probs: &[f64]) -> f64 {
    probs
        .iter()
        .filter(|&&p| p > 0.0) // Ignore zero probabilities
        .map(|&p| -p * p.log2())
        .sum()
}

/// Computes the maximum entropy threshold for a set of k-mers.
///
/// This function calculates the maximum entropy threshold for a set of k-mers.
/// The threshold is the frequency at which the entropy of the k-mer distribution
/// is maximized.
///
/// # Arguments
///
/// * `kmers` - A map where the key is a k-mer and the value is the frequency of that k-mer.
///
/// # Returns
///
/// The maximum entropy threshold.
///
/// # Examples
///
/// ```
/// use rustc_hash::FxHashMap;
/// use utilrs::max_entropy;
///
/// let mut kmers: FxHashMap<String, u32> = FxHashMap::default();
/// kmers.insert("ACGT".to_string(), 10);
/// kmers.insert("ACGA".to_string(), 5);
/// kmers.insert("ACGC".to_string(), 3);
/// kmers.insert("ACGG".to_string(), 2);
///
/// let threshold = max_entropy(&kmers);
/// assert_eq!(threshold, 5);
/// ```
///
fn max_entropy(kmers: &FxHashMap<String, u32>) -> u32 {
    let mut data: Vec<u32> = kmers.values().cloned().collect();
    let total: u32 = data.iter().sum();

    data.sort_by(|a, b| b.cmp(a));

    let normalized_data: Vec<f64> = data.iter().map(|&x| x as f64 / total as f64).collect();

    let entropy_curve: Vec<f64> = (1..normalized_data.len())
        .into_par_iter()
        .map(|s| {
            // Region A: probs[:s], total probability P_A
            let p_a: f64 = normalized_data[..s].iter().sum();
            let h_a = if p_a > 0.0 {
                let p_a_data: Vec<f64> = normalized_data[..s]
                    .iter()
                    .map(|&x| x as f64 / p_a as f64)
                    .collect();
                entropy(&p_a_data)
            } else {
                0.0
            };
            // Region B: probs[s:], total probability P_B
            let p_b: f64 = normalized_data[s..].iter().sum();
            let h_b = if p_b > 0.0 {
                let p_b_data: Vec<f64> = normalized_data[s..]
                    .iter()
                    .map(|&x| x as f64 / p_b as f64)
                    .collect();
                entropy(&p_b_data)
            } else {
                0.0
            };
            h_a + h_b
        })
        .collect();

    let (_max_entropy, max_entropy_idx) =
        entropy_curve
            .iter()
            .enumerate()
            .fold((f64::MIN, 0), |(max_val, max_idx), (idx, &val)| {
                if val > max_val {
                    (val, idx)
                } else {
                    (max_val, max_idx)
                }
            });

    let threshold = max_entropy_idx;
    let frequency = data[threshold];

    frequency
}

/// Selects k-mers with frequencies above a given threshold.
///
/// This function filters a map of k-mers by selecting only those with frequencies
/// above a specified threshold.
///
/// # Arguments
///
/// * `kmers` - A map where the key is a k-mer and the value is the frequency of that k-mer.
/// * `threshold` - The minimum frequency threshold for selecting k-mers.
///
/// # Returns
///
/// A map containing only the k-mers with frequencies above the threshold.
///
/// # Examples
///
/// ```
/// use rustc_hash::FxHashMap;
/// use utilrs::select_kmers;
///
/// let mut kmers: FxHashMap<String, u32> = FxHashMap::default();
/// kmers.insert("ACGT".to_string(), 10);
/// kmers.insert("ACGA".to_string(), 5);
/// kmers.insert("ACGC".to_string(), 3);
/// kmers.insert("ACGG".to_string(), 2);
///
/// let threshold = 5;
/// let selected_kmers = select_kmers(kmers, threshold);
/// assert_eq!(selected_kmers.len(), 2);
/// assert_eq!(selected_kmers.get("ACGT"), Some(&1));
/// assert_eq!(selected_kmers.get("ACGA"), Some(&1));
/// ```
///
fn select_kmers(kmers: FxHashMap<String, u32>, threshold: u32) -> FxHashMap<String, u32> {
    let mut selected_kmers: FxHashMap<String, u32> = FxHashMap::default();
    for (key, value) in kmers.iter() {
        if value >= &threshold {
            selected_kmers.insert(key.clone(), 1);
        }
    }
    selected_kmers
}

///
/// Processes a sequence to find mutations in k-mers.
///
/// This function processes a sequence to find mutations in k-mers that are exclusive to the sequence.
/// It uses parallel processing to find mutations in k-mers that are not present in a reference sequence.
///
/// # Arguments
///
/// * `seq` - A string slice representing the sequence to process.
/// * `ref_seq_str` - A reference sequence as a string slice.
/// * `exclusive_kmers_set` - A set of k-mers that are exclusive to the sequence.
/// * `k` - The length of k-mers to process.
/// * `step` - The step size for extracting overlapping k-mers.
/// * `max_dist` - The maximum number of allowed mutations in a k-mer.
/// * `ref_end` - The end index of the reference sequence.
///
/// # Returns
///
/// A vector of strings representing the mutations found in the sequence.
///
/// # Examples
///
/// ```
/// use std::sync::Arc;
/// use std::collections::HashSet;
/// use utilrs::process_sequence;
///
/// let seq = "ACGTACGTACGT";
/// let ref_seq_str = Arc::new("ACGTACGTACGT".to_string());
/// let exclusive_kmers_set = Arc::new(HashSet::new());
/// let k = 4;
/// let step = 1;
/// let max_dist = 1;
/// let ref_end = 12;
///
/// let mutations = process_sequence(seq, &ref_seq_str, &exclusive_kmers_set, k, step, max_dist, ref_end);
/// assert_eq!(mutations, vec![]);
/// ```
///
fn process_sequence(
    seq: &str,
    ref_seq_str: &Arc<String>,
    exclusive_kmers_set: &Arc<HashSet<String>>,
    k: usize,
    step: usize,
    max_dist: usize,
    ref_end: usize,
) -> Vec<String> {
    let splitted_seq: HashSet<String> = split_seq(seq, k, step).into_iter().collect();
    let var_kmers: Vec<&String> = exclusive_kmers_set.intersection(&splitted_seq).collect();

    let _sp = Spinner::new(Spinners::GrowVertical, String::new());

    var_kmers
        .into_par_iter()
        .flat_map(|var_kmer| {
            let pattern = var_kmer.as_bytes();
            (0..ref_end)
                .into_par_iter()
                .step_by(step)
                .filter_map(move |mutation_index| {
                    let ref_kmer = &ref_seq_str[mutation_index..mutation_index + k].as_bytes();
                    if levenshtein(pattern, ref_kmer) <= max_dist.try_into().unwrap() {
                        let ref_seq_kmer = &ref_seq_str[mutation_index..mutation_index + k];
                        Some(
                            get_kmer_mutation_index(pattern, ref_kmer, max_dist)
                                .into_par_iter()
                                .map(move |(snp, ind)| {
                                    concat_string!(
                                        (mutation_index + 1 + ind).to_string(),
                                        ":",
                                        snp,
                                        ":",
                                        ref_seq_kmer,
                                        ":",
                                        var_kmer
                                    )
                                }),
                        )
                    } else {
                        None
                    }
                })
                .flat_map(|x| x)
        })
        .collect()
}

/// Get the kmers from a sequence file
///
/// This function reads a FASTA file, processes sequences in batches, and counts the occurrences of k-mers of a specified length (`k`).
/// It uses multi-threading to handle large files efficiently. The result is a vector of `FxHashMap`s, where each map contains k-mer counts
/// and additional information based on the provided arguments.
///
/// # Arguments
///
/// * `seq_path` - The path to the FASTA file containing the sequences to process.
/// * `k` - The length of k-mers to count.
/// * `step` - The step size for extracting overlapping k-mers.
/// * `dict` - A string specifying the dictionary type ('DNA', 'RNA', or 'ALL') used to determine forbidden characters.
/// * `reference` - A boolean indicating if the reference is the input
/// * `batch_size` - The number of sequences to process per batch.
///
/// # Returns
///
/// A `FxHashMap<String, u32>` where each `FxHashMap` contains:
/// - k-mer frequencies.
///
///
/// Example
/// ```
/// use utilrs::get_kmers;
/// let seq_path = "path/to/sequences.fasta";
/// let k = 4;
/// let step = 1;
/// let dict = "DNA";
/// let reference = false;
/// let batch_size = 100;
/// let kmers = get_kmers(seq_path, k, step, dict, reference, batch_size);
/// ```
///
fn get_kmers_impl(
    seq_path: String,
    k: usize,
    step: usize,
    dict: String,
    reference: bool,
    batch_size: usize,
) -> Result<FxHashMap<String, u32>> {
    let _sp = Spinner::new(Spinners::GrowVertical, String::new());
    let reader = fasta::Reader::from_file(seq_path)?;
    let dict = Arc::new(dict);

    let sequences: Vec<_> = reader
        .records()
        .filter_map(|result| {
            let record = result.ok()?;
            let seq = String::from_utf8(record.seq().to_owned()).ok()?;
            if !contains_forbidden_chars(&seq, &dict) {
                Some(seq)
            } else {
                None
            }
        })
        .collect();

    let hm: FxHashMap<String, u32> = sequences
        .par_chunks(batch_size)
        .map(|batch| {
            let mut batch_hm: FxHashMap<String, u32> = FxHashMap::default();
            for seq in batch {
                let end = seq.len() - k + 1;
                for index in (0..end).step_by(step) {
                    *batch_hm
                        .entry(seq[index..index + k].to_string())
                        .or_insert(0) += 1;
                }
            }
            batch_hm
        })
        .reduce(
            || FxHashMap::default(),
            |mut acc, other| {
                for (key, value) in other {
                    *acc.entry(key).or_insert(0) += value;
                }
                acc
            },
        );

    let selected_kmers = if reference {
        select_kmers(hm, 0)
    } else {
        let threshold = max_entropy(&hm);
        select_kmers(hm, threshold)
    };

    Ok(selected_kmers)
}

#[pyfunction]
fn variants_intersection(
    variants_exclusive_kmers: FxHashMap<String, Vec<String>>,
    variants_names: Vec<String>,
    intersection_seletion: String,
) -> Py<PyAny> {
    let mut intersection_kmers: FxHashMap<String, Vec<String>> = FxHashMap::default();
    let mut intersection_kmers_sets: FxHashMap<String, Vec<String>> = FxHashMap::default();

    if intersection_seletion == "ALL" {
        for r in 1..variants_names.len() + 1 {
            for subset in variants_names.clone().into_iter().combinations(r) {
                if subset.len() > 1 {
                    let intersection_set =
                        get_set_intersection(subset.clone(), variants_exclusive_kmers.clone());
                    if intersection_set.len() > 0 {
                        for variant in &subset {
                            intersection_kmers_sets
                                .insert(variant.clone(), intersection_set.clone());
                        }
                        intersection_kmers.insert(subset.join("-"), intersection_set);
                    }
                }
            }
        }
    } else {
        let intersection_selects = intersection_seletion
            .split("-")
            .map(|s| s.to_string())
            .collect::<Vec<String>>();
        let intersection_set = get_set_intersection(
            intersection_selects.clone(),
            variants_exclusive_kmers.clone(),
        );
        if intersection_set.len() > 0 {
            for variant in &intersection_selects {
                intersection_kmers_sets.insert(variant.clone(), intersection_set.clone());
            }
            intersection_kmers.insert(intersection_seletion, intersection_set);
        }
    }

    return Python::with_gil(|py| (intersection_kmers, intersection_kmers_sets).to_object(py));
}

#[pyfunction]
fn write_ref(ref_path: String, variations: Vec<String>, save_path: String) -> () {
    if variations.is_empty() {
        return;
    }
    let reader = fasta::Reader::from_file(ref_path).unwrap();

    for result in reader.records() {
        let result_data = result.unwrap();
        let seq: &[u8] = result_data.seq();
        let seq_name = concat_string!(result_data.id().to_string(), String::from("_reference"));
        let seq_str = String::from_utf8(seq.to_owned()).unwrap();
        let mut seq_str = seq_str.to_uppercase();

        for variation in &variations {
            let variation = variation.split(":").collect::<Vec<&str>>();
            let position = variation[0].parse::<usize>().unwrap();
            let snp = variation[1].to_string();
            seq_str.replace_range(
                position..position + 1,
                &snp.chars().nth(1).unwrap().to_string(),
            );
        }

        let mut writer = fasta::Writer::to_file(save_path.clone()).unwrap();
        let record = fasta::Record::with_attrs(&seq_name, None, seq_str.as_bytes());
        let _ = writer.write_record(&record);
    }

    return;
}

#[pyfunction]
fn ref_length(ref_path: String) -> Option<usize> {
    let reader = fasta::Reader::from_file(ref_path).unwrap();
    for result in reader.records() {
        let seq_data = result.unwrap();
        return Some(seq_data.seq().len());
    }
    return None;
}

#[pyfunction]
fn get_freq_kmers(diffs: FxHashMap<String, Vec<String>>) -> Py<PyAny> {
    let mut freq_dict: FxHashMap<String, HashSet<(&str, String)>> = FxHashMap::default();
    let mut freqs: FxHashMap<&str, usize> = FxHashMap::default();
    let mut var_list: HashSet<&str> = HashSet::new();
    for (key, value) in diffs.iter() {
        for snp in value {
            let snp = snp.split(":").collect::<Vec<&str>>();
            let position = snp[0].to_string();
            let snp = snp[1].to_string();
            let position = concat_string!(position, String::from(":"), snp);
            freq_dict
                .entry(position)
                .and_modify(|freq| {
                    freq.insert((&key, snp.clone()));
                })
                .or_insert({
                    let mut freq_seq: HashSet<(&str, String)> = HashSet::new();
                    freq_seq.insert((&key, snp.clone()));
                    freq_seq
                });
        }
    }

    for (key, value) in freq_dict.iter() {
        freqs.insert(key, value.len());
        var_list.insert(key);
    }

    let var_list = Vec::from_iter(var_list);

    return Python::with_gil(|py| (freqs, var_list).to_object(py));
}

#[pyfunction]
fn kmers_analysis(
    seq_path: String,
    ref_path: String,
    exclusive_kmers: Vec<String>,
    k: usize,
    step: usize,
    max_dist: usize,
    batch_size: usize,
) -> PyResult<Py<PyAny>> {
    Python::with_gil(|py| {
        // Load reference sequence
        let ref_reader = fasta::Reader::from_file(&ref_path)
            .map_err(|e| PyErr::new::<pyo3::exceptions::PyIOError, _>(e.to_string()))?
            .records()
            .next()
            .ok_or_else(|| {
                PyErr::new::<pyo3::exceptions::PyValueError, _>("Reference sequence not found")
            })?
            .map_err(|e| PyErr::new::<pyo3::exceptions::PyValueError, _>(e.to_string()))?;
        let ref_seq_str = Arc::new(
            String::from_utf8(ref_reader.seq().to_owned())
                .map_err(|e| PyErr::new::<pyo3::exceptions::PyValueError, _>(e.to_string()))?
                .to_uppercase(),
        );
        let ref_end = ref_seq_str.chars().count() - k + step;

        let exclusive_kmers_set: Arc<HashSet<String>> =
            Arc::new(exclusive_kmers.into_iter().collect());

        // Create a vector for reading sequences
        let sequences = fasta::Reader::from_file(&seq_path)
            .map_err(|e| PyErr::new::<pyo3::exceptions::PyIOError, _>(e.to_string()))?
            .records()
            .map(|r| {
                r.map(|record| {
                    (
                        String::from_utf8(record.seq().to_owned())
                            .unwrap()
                            .to_uppercase(),
                        record.id().to_string(),
                    )
                })
            })
            .collect::<Result<Vec<_>, _>>()
            .map_err(|e| PyErr::new::<pyo3::exceptions::PyValueError, _>(e.to_string()))?;

        // Process sequences in batches, using parallel processing
        let variations_positions: FxHashMap<String, Vec<String>> = sequences
            .par_chunks(batch_size)
            .flat_map(|batch| {
                batch
                    .par_iter()
                    .map(|(seq, seq_name)| {
                        let variations = process_sequence(
                            seq,
                            &ref_seq_str,
                            &exclusive_kmers_set,
                            k,
                            step,
                            max_dist,
                            ref_end,
                        );
                        (seq_name.clone(), variations)
                    })
                    .collect::<Vec<_>>()
            })
            .collect();

        Ok(variations_positions.to_object(py))
    })
}

#[pyfunction]
fn get_kmers(
    seq_path: String,
    k: usize,
    step: usize,
    dict: String,
    reference: bool,
    batch_size: usize,
) -> PyResult<Py<PyAny>> {
    let result = get_kmers_impl(seq_path, k, step, dict, reference, batch_size);
    match result {
        Ok(selected_kmers) => Python::with_gil(|py| Ok(selected_kmers.to_object(py))),
        Err(e) => Err(PyErr::new::<pyo3::exceptions::PyRuntimeError, _>(
            e.to_string(),
        )),
    }
}

#[pyfunction]
fn load_sequences_classify(
    seqs_path: Vec<String>,
    k: usize,
    step: usize,
    dict: String,
    variants_kmers: Vec<String>,
    predict_data: bool,
    batch_size: usize,
) -> Py<PyAny> {
    let batch_size = batch_size.div_ceil(seqs_path.len());

    let samples: Vec<FxHashMap<String, Value>> = seqs_path
        .into_par_iter()
        .map(|seq_path| {
            count_freq(
                seq_path,
                k,
                step,
                dict.clone(),
                variants_kmers.clone(),
                predict_data,
                batch_size,
            )
        })
        .flatten()
        .collect::<Vec<FxHashMap<String, Value>>>();

    return Python::with_gil(|py| samples.to_object(py));
}

#[pymodule]
fn utilrs(_py: Python, m: &Bound<'_, PyModule>) -> PyResult<()> {
    m.add_function(wrap_pyfunction!(get_kmers, m)?)?;
    m.add_function(wrap_pyfunction!(kmers_analysis, m)?)?;
    m.add_function(wrap_pyfunction!(get_freq_kmers, m)?)?;
    m.add_function(wrap_pyfunction!(load_sequences_classify, m)?)?;
    m.add_function(wrap_pyfunction!(write_ref, m)?)?;
    m.add_function(wrap_pyfunction!(ref_length, m)?)?;
    m.add_function(wrap_pyfunction!(variants_intersection, m)?)?;
    Ok(())
}

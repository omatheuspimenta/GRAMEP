//! # Utilrs
//!
//! `utilrs` is a collection of utility functions used in GRAMEP.
//!
#![allow(unused_imports)]
#[macro_use(concat_string)]
extern crate concat_string;
use bio::alignment::distance::simd::*;
use bio::io::fasta;
use pyo3::prelude::*;
use rayon::prelude::*;
use regex::Regex;
use rustc_hash::FxHashMap;
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
    return alphabet;
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
    return kmers;
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
    return mutations;
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
    let num_threads: usize = num_cpus::get() - 1;

    let class_name = seq_path.split(".fasta").collect::<Vec<&str>>()[0];
    let class_name = class_name.split("/").collect::<Vec<&str>>();
    let class_name = class_name.last().unwrap().to_string();

    let reader = fasta::Reader::from_file(seq_path).unwrap();
    let pool: rayon::ThreadPool = rayon::ThreadPoolBuilder::new()
        .num_threads(num_threads)
        .build()
        .unwrap();
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

    return samples;
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
        freqs.insert(&key, value.len());
        var_list.insert(&key);
    }

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
) -> Py<PyAny> {
    let num_threads: usize = num_cpus::get() - 1;

    // Load reference sequence
    let reader = fasta::Reader::from_file(seq_path).unwrap();
    let ref_reader = fasta::Reader::from_file(ref_path)
        .unwrap()
        .records()
        .next()
        .unwrap()
        .unwrap();
    let ref_seq: &[u8] = ref_reader.seq();
    let ref_seq_str = String::from_utf8(ref_seq.to_owned()).unwrap();
    let ref_seq_str = Arc::new(ref_seq_str.clone().to_uppercase());
    let ref_end = ref_seq_str.chars().count() - k + step;

    let pool: rayon::ThreadPool = rayon::ThreadPoolBuilder::new()
        .num_threads(num_threads)
        .build()
        .unwrap();
    let (tx, rx) = std::sync::mpsc::channel();
    let exclusive_kmers_set: HashSet<String> = exclusive_kmers.into_iter().collect();

    let mut batch = Vec::new();
    for result in reader.records() {
        let result_data = result.unwrap();
        let seq: &[u8] = result_data.seq();
        let seq_name = result_data.id().to_string();
        let seq_str = String::from_utf8(seq.to_owned()).unwrap();
        let seq_str = seq_str.to_uppercase();

        batch.push((seq_str, seq_name));

        if batch.len() >= batch_size {
            let tx = tx.clone();
            let batch_clone = batch.clone();
            let exclusive_kmers_set = exclusive_kmers_set.clone();

            let ref_seq_str = Arc::clone(&ref_seq_str);

            batch.clear();

            pool.spawn(move || {
                let mut batch_hm: FxHashMap<String, Vec<String>> = FxHashMap::default();
                for (seq, seq_name) in batch_clone {
                    let splitted_seq: HashSet<String> =
                        split_seq(&seq, k, step).into_iter().collect();

                    let var_kmers = exclusive_kmers_set
                        .intersection(&splitted_seq)
                        .collect::<Vec<&String>>();

                    let mut variations = HashSet::new();

                    for var_kmer in var_kmers {
                        let pattern = var_kmer.as_bytes();
                        let mut mutation_index = 0;
                        let ref_step = 1;

                        while mutation_index + ref_step < ref_end {
                            let ref_kmer =
                                &ref_seq_str[mutation_index..mutation_index + k].as_bytes();
                            let ldist = levenshtein(pattern, &ref_kmer);

                            if ldist <= max_dist.try_into().unwrap() {
                                let ref_seq_kmer = &ref_seq_str[mutation_index..mutation_index + k];
                                let find_mut =
                                    get_kmer_mutation_index(pattern, &ref_kmer, max_dist);
                                for (snp, ind) in find_mut {
                                    let diff = concat_string!(
                                        (mutation_index + 1 + ind).to_string(),
                                        String::from(":"),
                                        snp,
                                        String::from(":"),
                                        ref_seq_kmer.to_string(),
                                        String::from(":"),
                                        var_kmer.to_string()
                                    );
                                    variations.insert(diff);
                                }
                            }
                            mutation_index += ref_step;
                        }
                    }
                    let variations_vec = Vec::from_iter(variations);
                    batch_hm.insert(seq_name, variations_vec);
                }
                tx.send(batch_hm).unwrap();
            });
        }
    }

    if !batch.is_empty() {
        let tx = tx.clone();

        let exclusive_kmers_set = exclusive_kmers_set.clone();

        pool.spawn(move || {
            let mut batch_hm: FxHashMap<String, Vec<String>> = FxHashMap::default();
            for (seq, seq_name) in batch {
                let splitted_seq: HashSet<String> = split_seq(&seq, k, step).into_iter().collect();

                let var_kmers = exclusive_kmers_set
                    .intersection(&splitted_seq)
                    .collect::<Vec<&String>>();

                let mut variations = HashSet::new();

                for var_kmer in var_kmers {
                    let pattern = var_kmer.as_bytes();
                    let mut mutation_index = 0;
                    let ref_step = 1;

                    while mutation_index + ref_step < ref_end {
                        let ref_kmer = &ref_seq_str[mutation_index..mutation_index + k].as_bytes();
                        let ldist = levenshtein(pattern, &ref_kmer);

                        if ldist <= max_dist.try_into().unwrap() {
                            let ref_seq_kmer = &ref_seq_str[mutation_index..mutation_index + k];
                            let find_mut = get_kmer_mutation_index(pattern, &ref_kmer, max_dist);
                            for (snp, ind) in find_mut {
                                let diff = concat_string!(
                                    (mutation_index + 1 + ind).to_string(),
                                    String::from(":"),
                                    snp,
                                    String::from(":"),
                                    ref_seq_kmer.to_string(),
                                    String::from(":"),
                                    var_kmer.to_string()
                                );
                                variations.insert(diff);
                            }
                        }
                        mutation_index += ref_step;
                    }
                }
                let variations_vec = Vec::from_iter(variations);
                batch_hm.insert(seq_name, variations_vec);
            }
            tx.send(batch_hm).unwrap();
        });
    }

    drop(tx);

    let mut variations_positions: FxHashMap<String, Vec<String>> = FxHashMap::default();
    for received_hm in rx {
        for (key, value) in received_hm.into_iter() {
            variations_positions.insert(key, value);
        }
    }

    return Python::with_gil(|py| variations_positions.to_object(py));
}

#[pyfunction]
fn get_kmers(
    seq_path: String,
    k: usize,
    step: usize,
    dict: String,
    batch_size: usize,
) -> Py<PyAny> {
    let num_threads: usize = num_cpus::get() - 1;
    let reader = fasta::Reader::from_file(seq_path).unwrap();
    let pool: rayon::ThreadPool = rayon::ThreadPoolBuilder::new()
        .num_threads(num_threads)
        .build()
        .unwrap();
    let (tx, rx) = std::sync::mpsc::channel();

    // Process sequences in batches
    let mut batch = Vec::new();
    for result in reader.records() {
        let result_data = result.unwrap();
        let seq: &[u8] = result_data.seq();
        let seq_str = String::from_utf8(seq.to_owned()).unwrap();
        if !contains_forbidden_chars(&seq_str, &dict) {
            batch.push(seq_str);
            if batch.len() >= batch_size {
                let tx = tx.clone();
                let batch_clone = batch.clone();
                batch.clear();
                pool.spawn(move || {
                    let mut batch_hm: FxHashMap<String, u32> = FxHashMap::default();
                    for seq in batch_clone {
                        let mut hm: FxHashMap<String, u32> = FxHashMap::default();
                        let end = seq.chars().count() - k + step;
                        let mut index = 0;
                        while index + step < end {
                            *hm.entry(seq[index..index + k].to_owned()).or_insert(0) += 1;
                            index += step;
                        }
                        for (key, value) in hm.into_iter() {
                            *batch_hm.entry(key).or_insert(0) += value;
                        }
                    }
                    tx.send(batch_hm).unwrap();
                });
            }
        }
    }
    // Process any remaining sequences in the last batch
    if !batch.is_empty() {
        let tx = tx.clone();
        pool.spawn(move || {
            let mut batch_hm: FxHashMap<String, u32> = FxHashMap::default();
            for seq in batch {
                let mut hm: FxHashMap<String, u32> = FxHashMap::default();
                let end = seq.chars().count() - k + step;
                let mut index = 0;
                while index + step < end {
                    *hm.entry(seq[index..index + k].to_owned()).or_insert(0) += 1;
                    index += step;
                }
                for (key, value) in hm.into_iter() {
                    *batch_hm.entry(key).or_insert(0) += value;
                }
            }
            tx.send(batch_hm).unwrap();
        });
    }

    drop(tx); // Close all senders

    // Merge HashMaps incrementally
    let mut hm: FxHashMap<String, u32> = FxHashMap::default();
    for received_hm in rx {
        for (key, value) in received_hm.into_iter() {
            *hm.entry(key).or_insert(0) += value;
        }
    }

    return Python::with_gil(|py| hm.to_object(py));
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
    let batch_size = batch_size.div_ceil(seqs_path.len()) as usize;

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
    Ok(())
}

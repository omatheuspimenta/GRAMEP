"""
kmers_utils
===========

This module provides utility functions for working with k-mers.

Contents:
    * kmers_freq: Get the frequency of k-mers in a sequence.
    * select_kmers: Select k-mers from a dictionary based on a frequency threshold.
    * kmers_difference: Calculate the difference in k-mers between two defaultdicts.
    * kmers_intersections: Calculate the intersection in k-mers between two defaultdicts.

Note:
This module assumes that the input sequence contains characters from the \
DNA alphabet ('ACTG') by default. Custom alphabets can be specified using \
the 'dictionary' parameter.

Todo:
    * Implement tests.
"""
from collections import Counter, defaultdict


def kmers_freq(sequence: str, word: int, step: int) -> defaultdict[str, int]:
    """
    Calculate the frequency of overlapping k-mers in a sequence.

    This function takes a sequence string and calculates the frequency of overlapping
    k-mers of the specified length and step size. It returns a defaultdict that maps
    k-mers to their frequency counts.

    Args:
        sequence (str): The input sequence as a string.
        word (int): The length of each k-mer.
        step (int): The step size for moving the sliding window.

    Returns:
        defaultdict[str, int]: A defaultdict mapping k-mers to their frequency counts.
    """

    index = 0
    kmers: defaultdict[str, int] = defaultdict(int)

    while (index + word) < len(sequence):
        kmers[str(''.join(sequence[index : index + word]))] += 1
        index += step
    return kmers


def select_kmers(threshold: int, kmers: defaultdict) -> defaultdict[str, int]:
    """
    Select k-mers from a dictionary based on a frequency threshold.

    This function takes a dictionary of k-mers and their frequency counts, and returns
    a new dictionary containing only those k-mers whose frequency exceeds the specified
    threshold.

    Args:
        threshold (int): The frequency threshold for selecting k-mers.
        kmers (defaultdict[str, int]): A dictionary mapping k-mers to their\
        frequency counts.

    Returns:
        defaultdict[str, int]: A dictionary containing selected k-mers \
        and their frequencies.
    """

    return defaultdict(int, {k: 1 for k, v in kmers.items() if v > threshold})


def kmers_difference(
    seq_kmers: defaultdict, ref_kmers: defaultdict
) -> list[str]:
    """
    Calculate the difference in k-mers between two defaultdicts.

    This function takes two defaultdicts containing k-mers and their frequency counts,
    calculates the k-mer differences between them, and returns a list of k-mers that
    are present in the 'seq_kmers' defaultdict but not in the 'ref_kmers' defaultdict.

    Args:
        seq_kmers (defaultdict): A defaultdict mapping k-mers to their frequency counts
                                 in the sequence data.
        ref_kmers (defaultdict): A defaultdict mapping k-mers to their frequency counts
                                 in the reference data.

    Returns:
        list[str]: A list of k-mers present in 'seq_kmers' but not in 'ref_kmers'.
    """
    seq = Counter(seq_kmers)
    ref = Counter(ref_kmers)

    return list(seq - ref)


def kmers_intersections(
    seq_kmers: defaultdict, ref_kmers: defaultdict
) -> list[str]:
    """
    Calculate the intersection in k-mers between two defaultdicts.

    This function takes two defaultdicts containing k-mers and their frequency counts,
    calculates the k-mer intersection between them, and returns a list of k-mers that
    are present in both defaultdicts.

    Args:
        seq_kmers (defaultdict): A defaultdict mapping k-mers to their frequency counts
                                 in the sequence data.
        ref_kmers (defaultdict): A defaultdict mapping k-mers to their frequency counts
                                 in the reference data.

    Returns:
        list[str]: A list of k-mers present in both defaultdicts.
    """
    seq = Counter(seq_kmers)
    ref = Counter(ref_kmers)

    return list(seq & ref)

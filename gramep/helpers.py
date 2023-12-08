"""helpers module.

This module houses a collection of auxiliary functions
designed to provide general-purpose utility and support
across different parts of the GRAMEP package. These functions
are not domain-specific and are meant to assist with
common tasks that might be needed in various contexts.

Contents:
    * split_sequence: Split a sequence into overlapping k-mers.
    * get_sequence_interval: Create a closed interval using the 'start' and 'end' \
    values from a row of data.
    * dicts_sum: Combine multiple dictionaries by summing their corresponding values.
    * dicts_append: Combine and append values from multiple dictionaries \
    into a defaultdict.
    * dicts_drop_duplicates: Remove duplicate values from a dictionary's \
    value lists and return a new defaultdict.
    * check_diffs: Check differences between reference and variant k-mers.
    * batched: Generate batches of items from an iterable.
    * make_report: Generate a report for variations and annotations in a sequence.
    * get_kmers_frequencies: Get k-mer frequencies and positions from differences.
    * get_colours: Get a dictionary of colours based on the specified palette.
    * next_colour: Get the next colour from the predefined colour cycle.
    * count_freq: Count the frequency of items in a sequence.
    * get_newick_list: Construct a Newick tree string from a SciPy hierarchical \
    clustering ClusterNode.
    * get_newick_str: Generate a Newick tree string from a SciPy hierarchical \
    clustering tree.



Todo:
    * Implement tests.
"""
from collections import defaultdict
from itertools import cycle, islice
from typing import Generator

import pandas as pd
from scipy.cluster.hierarchy import ClusterNode

from gramep.messages import Messages

message = Messages()
"""
Set the Message class for logging.
"""

COLOR_LIST = ['lightgrey', 'white']
"""
list[str]: A list of color values for cycling through in the color cycle.
"""
COLOR_CYCLE = cycle(COLOR_LIST)
"""
itertools.cycle: A cycle iterator that iterates through the COLOR_LIST repeatedly.
"""


def split_sequence(sequence: str, word: int, step: int) -> list[str]:
    """
    Split a sequence into overlapping k-mers of a specified length and step size.

    Args:
        sequence (str): The input sequence to be split.
        word (int): The length of each k-mer.
        step (int): The step size to move the sliding window.

    Returns:
        list[str]: A list of overlapping k-mers extracted from the sequence.
    """
    index = 0
    kmers = []
    while (index + word) < len(sequence):
        kmers.append(str(''.join(sequence[index : index + word])))
        index += step
    return kmers


def get_sequence_interval(row: pd.Series) -> pd.Interval:
    """
    Create a closed interval using the 'start' and 'end' values from a row of data.

    This function constructs a closed interval using the 'start' and 'end' values from
    a pandas Series row. The interval is closed on both ends.

    Args:
        row (pd.Series): A pandas Series containing 'start' and 'end' values.

    Returns:
        pd.Interval: A closed interval defined by the 'start' and 'end' values.
    """
    return pd.Interval(row.start, row.end, closed='both')


def dicts_sum(*dicts: defaultdict[str, int]) -> defaultdict[str, int]:
    """
    Combine multiple dictionaries by summing their corresponding values.

    This function takes a variable number of dictionaries and returns a new defaultdict
    where the values corresponding to the same keys are summed across all \
    input dictionaries.

    Args:
        *dicts (defaultdict[str, int]): Variable number of dictionaries to be combined.

    Returns:
        defaultdict[str, int]: A defaultdict containing the summed values for each key.
    """
    dict_sum = defaultdict(int)
    for dictionary in dicts:
        for key, value in dictionary.items():
            dict_sum[key] += value
    return dict_sum


def dicts_append(
    *dicts: defaultdict[str, list[str]]
) -> defaultdict[str, list[str]]:
    """
    Combine and append values from multiple dictionaries into a defaultdict.

    This function takes a variable number of dictionaries and creates a new defaultdict
    that combines and appends values from the input dictionaries based on their keys.

    Args:
        *dicts: Variable number of dictionaries to be combined and appended.

    Returns:
        defaultdict[str, list[str]]: A defaultdict containing combined and \
        appended values as lists for each key.
    """
    dict_append = defaultdict(list)
    for dictionary in dicts:
        for key, value in dictionary.items():
            for i in value:
                dict_append[key].append(i)
    return dict_append


def dicts_drop_duplicates(
    dictionary: defaultdict[str, list[str]]
) -> defaultdict[str, list[str]]:
    """
    Remove duplicate values from a dictionary's value lists and return \
    a new defaultdict.

    This function takes a dictionary and creates a new defaultdict \
    where duplicate values in the original dictionary's value lists are removed.

    Args:
        dictionary(defaultdict[str, list[str]]): A dictionary containing keys and\
            associated value lists.

    Returns:
        defaultdict[str, list[str]]: A defaultdict with duplicate-free value \
        lists for each key.
    """
    dict_no_dup = defaultdict(list)
    for key, value in dictionary.items():
        dict_no_dup[key] = list(set(value))
    return dict_no_dup


def check_diffs(ref_kmers: str, var_kmers: str, ref_index: int) -> list[str]:
    """
    Check differences between reference and variant k-mers.

    This function takes a reference k-mer string, a variant k-mer string,
    and the index of the reference k-mer in the context of the data. It checks
    for differences between the two k-mer strings and returns a list of strings
    representing the differences found.

    Args:
        ref_kmers (str): The reference k-mer string.
        var_kmers (str): The variant k-mer string.
        ref_index (int): The index of the reference k-mer in the data.

    Returns:
        list[str]: A list of strings representing differences found.
    """

    return [
        str(
            str(int(ref_index) + index + 1)
            + ':'
            + str(ref_kmers[index] + char)
            + ':'
            + str(ref_kmers)
            + ':'
            + str(var_kmers)
        )
        for index, char in enumerate(var_kmers)
        if not ref_kmers[index] == char
    ]


def batched(
    iterable: Generator[list[tuple[str, str]], None, None], n: int
) -> Generator[list[tuple[str, str]], None, None]:
    """
    Generate batches of items from an iterable.

    This function takes an iterable of items, such as a generator, and yields batches
    of items as lists of tuples. Each batch contains 'n' items from the input iterable.

    Args:
        iterable (Generator[list[tuple[str, str]], None, None]): An iterable of items
            represented as lists of tuples.
        n (int): The size of each batch.

    Yields:
        Generator[list[tuple[str, str]], None, None]: A generator that yields batches \
        of items as lists of tuples.

    Raises:
        ValueError: If 'n' is less than one.
    """
    if n < 1:
        raise ValueError('n must be at least one')
    it = iter(iterable)
    while batch := list(islice(it, n)):
        yield batch


def make_report(
    variation_position: list,
    seq_name: str,
    sequence_interval: pd.IntervalIndex,
    annotation_dataframe: pd.DataFrame,
) -> list[str]:
    """
    Generate a report for variations and annotations in a sequence.

    This function takes variation positions, sequence name, sequence intervals,
    and an annotation DataFrame, and generates a report detailing variations and
    annotations present in the sequence.

    Args:
        variation_position (list): A list of variation positions in the sequence.
        seq_name (str): The name of the sequence being analyzed.
        sequence_interval (pd.IntervalIndex): Interval index representing \
        sequence intervals.
        annotation_dataframe (pd.DataFrame): DataFrame containing sequence annotations.

    Returns:
        list[str]: A list of strings representing the generated report.
    """

    report_list = []
    for var in variation_position:
        variation_index, snp_value, ref_value, var_value = var.split(':')
        if sequence_interval is not None:
            annotation_index = []
            for interval in sequence_interval:
                if int(variation_index) in interval:
                    annotation_index.append(
                        list(
                            sequence_interval[
                                sequence_interval == interval
                            ].index
                        )
                    )
            annotation_index = [
                auxIndex1[auxIndex2]
                for auxIndex1 in list(set(map(tuple, annotation_index)))
                for auxIndex2 in range(len(auxIndex1))
            ]

            if len(annotation_index) > 0:
                for index_element in annotation_index:
                    report_line = str()
                    if (
                        annotation_dataframe.loc[index_element]['product']
                        is None
                    ):
                        report_line = str(
                            str(seq_name)
                            + ';'
                            + str(
                                annotation_dataframe.loc[index_element]['Name']
                            )
                            + ';'
                            + str(
                                annotation_dataframe.loc[index_element][
                                    'start'
                                ]
                            )
                            + ';'
                            + str(
                                annotation_dataframe.loc[index_element]['end']
                            )
                            + ';'
                            + str(
                                annotation_dataframe.loc[index_element]['type']
                            )
                            + ';'
                            + str(variation_index)
                            + ';'
                            + str(ref_value)
                            + ';'
                            + str(var_value)
                            + ';'
                            + str(snp_value[0])
                            + ';'
                            + str(snp_value[1])
                        )
                    report_list.append(report_line)
            else:
                report_line = str(
                    str(seq_name)
                    + ';'
                    + str('NO_ANNOTATIONS_AVAILABLE')
                    + ';'
                    + str('NO_ANNOTATIONS_AVAILABLE')
                    + ';'
                    + str('NO_ANNOTATIONS_AVAILABLE')
                    + ';'
                    + str('NO_ANNOTATIONS_AVAILABLE')
                    + ';'
                    + str(variation_index)
                    + ';'
                    + str(ref_value)
                    + ';'
                    + str(var_value)
                    + ';'
                    + str(snp_value[0])
                    + ';'
                    + str(snp_value[1])
                )
                report_list.append(report_line)
        else:
            report_line = str(
                str(seq_name)
                + ';'
                + str('NO_ANNOTATIONS_AVAILABLE')
                + ';'
                + str('NO_ANNOTATIONS_AVAILABLE')
                + ';'
                + str('NO_ANNOTATIONS_AVAILABLE')
                + ';'
                + str('NO_ANNOTATIONS_AVAILABLE')
                + ';'
                + str(variation_index)
                + ';'
                + str(ref_value)
                + ';'
                + str(var_value)
                + ';'
                + str(snp_value[0])
                + ';'
                + str(snp_value[1])
            )
            report_list.append(report_line)

    return list(set(report_list))


def get_kmers_frequencies(
    diffs_positions: defaultdict,
) -> tuple[defaultdict[str, list[str]], list[str]]:
    """
    Get k-mer frequencies and positions from differences.

    This function takes a defaultdict containing differences positions and \
    returns a tuple containing a defaultdict of k-mer frequencies and a list of k-mers.

    Args:
        diffs_positions (defaultdict): A defaultdict mapping \
        differences to their positions.

    Returns:
        tuple[defaultdict[str, list[str]], list[str]]: A tuple containing a \
        defaultdict of k-mer frequencies and a list of k-mers.
    """
    freq_dict_temp = defaultdict(list)
    freq_dict = defaultdict(int)
    var_list = []

    for key, value in diffs_positions.items():
        snps = value
        for snp in snps:
            pos, var, _, _ = snp.split(':')
            ref = var[0]
            base = var[1]
            freq_dict_temp[str(str(pos) + ':' + str(ref) + str(base))].append(
                (key, str(str(ref) + str(base)))
            )

    freq_dict_temp = dicts_drop_duplicates(freq_dict_temp)

    for key, value in freq_dict_temp.items():
        freq_dict[key] = len(value)
        var_list.append(key)
    var_list = list(set(var_list))
    return freq_dict, var_list


def get_colours(colour_palette: str) -> dict[str, dict[str, str]]:
    """
    Get a dictionary of colours based on the specified palette.

    This function takes a colour palette name and returns a dictionary mapping
    predefined keys to nested dictionaries containing colour values based on \
    the chosen palette.

    Args:
        colour_palette (str): The name of the colour palette to use.

    Returns:
        dict[str, dict[str, str]]: A dictionary mapping predefined keys to \
        nested dictionaries containing colour values.
    """

    palettes = {
        'classic': {
            'A': 'steelblue',
            'C': 'indianred',
            'T': 'darkseagreen',
            'G': 'skyblue',
        },
        'wes': {
            'A': '#CC8B3C',
            'C': '#456355',
            'T': '#541F12',
            'G': '#B62A3D',
        },
        'primary': {
            'A': 'green',
            'C': 'goldenrod',
            'T': 'steelblue',
            'G': 'indianred',
        },
        'purine-pyrimidine': {
            'A': 'indianred',
            'C': 'teal',
            'T': 'teal',
            'G': 'indianred',
        },
        'greyscale': {
            'A': '#CCCCCC',
            'C': '#999999',
            'T': '#666666',
            'G': '#333333',
        },
        'blues': {
            'A': '#3DB19D',
            'C': '#76C5BF',
            'T': '#423761',
            'G': 'steelblue',
        },
        'verity': {
            'A': '#EC799A',
            'C': '#df6eb7',
            'T': '#FF0080',
            'G': '#9F0251',
        },
    }
    if colour_palette not in palettes:
        exit(message.error_color_pallete())
    else:
        color_dict = palettes[colour_palette]

    return color_dict


def next_colour() -> str:
    """
    Get the next colour from the predefined colour cycle.

    This function returns the next colour from the predefined colour cycle,
    which is used to cycle through a list of colour values.

    Returns:
        str: The next colour from the colour cycle.
    """
    return next(COLOR_CYCLE)


def count_freq(seq: list) -> defaultdict[str, int]:
    """
    Count the frequency of items in a sequence.

    This function takes a sequence and returns a defaultdict containing \
    the frequency count of each unique item in the sequence.

    Args:
        seq (list): The sequence of items to count frequencies for.

    Returns:
        defaultdict[str, int]: A defaultdict containing the frequency count of \
        each unique item.
    """
    freq = defaultdict(int)
    for item in seq:
        freq[item] += 1
    return freq


def get_newick_list(
    node: ClusterNode,
    newick: list[str],
    parentdist: float,
    leaf_names: list[str],
) -> list[str]:
    """
    Construct a Newick tree string from a SciPy hierarchical clustering ClusterNode.

    This recursive function aids in building a Newick output string from a scipy.cluster.hierarchy.to_tree input
    with user-specified leaf node names.

    Notes:
        This function is intended for use with the `to_newick` function.

    Args:
        node (scipy.cluster.hierarchy.ClusterNode): The root node obtained from scipy.cluster.hierarchy.to_tree,
            representing the hierarchical clustering linkage matrix.
        parentdist (float): The distance of the parent node of the current `node`.
        newick (list of string): An accumulator list for the Newick string output.
            It needs to be reversed and concatenated (i.e., `''.join(newick)`) for the final output.
        leaf_names (list of string): The names of leaf nodes.

    Returns:
        list of string: The `newick` list containing Newick output strings.
    """
    if node.is_leaf():
        return newick + [f'{leaf_names[node.id]}:{parentdist - node.dist}']

    if len(newick) > 0:
        newick.append(f'):{parentdist - node.dist}')
    else:
        newick.append(');')
    newick = get_newick_list(node.get_left(), newick, node.dist, leaf_names)
    newick.append(',')
    newick = get_newick_list(node.get_right(), newick, node.dist, leaf_names)
    newick.append('(')
    return newick


def get_newick_str(tree: ClusterNode, leaf_names: list[str]) -> str:
    """
    Generate a Newick tree string from a SciPy hierarchical clustering tree.

    This function converts a SciPy ClusterNode tree to a Newick format string. To use this function, first obtain
    the root ClusterNode by applying scipy.cluster.hierarchy.to_tree on a hierarchical clustering linkage matrix.

    Args:
        tree (scipy.cluster.hierarchy.ClusterNode): The root node obtained from scipy.cluster.hierarchy.to_tree,
            representing the hierarchical clustering linkage matrix.
        leaf_names (list of string): The names of leaf nodes.

    Returns:
        str: The Newick output string.
    """
    newick_list = get_newick_list(tree, [], tree.dist, leaf_names)
    return ''.join(newick_list[::-1])

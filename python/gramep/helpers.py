"""helpers module.

This module houses a collection of auxiliary functions
designed to provide general-purpose utility and support
across different parts of the GRAMEP package. These functions
are not domain-specific and are meant to assist with
common tasks that might be needed in various contexts.

Contents:
    * get_sequence_interval: Create a closed interval using the 'start' and 'end' \
    values from a row of data.
    * make_report: Generate a report for variations and annotations in a sequence.
    * get_colours: Get a dictionary of colours based on the specified palette.
    * next_colour: Get the next colour from the predefined colour cycle.
    * get_newick_list: Construct a Newick tree string from a SciPy hierarchical \
    clustering ClusterNode.
    * get_newick_str: Generate a Newick tree string from a SciPy hierarchical \
    clustering tree.

Todo:
    * Implement tests.
"""
from itertools import cycle

import pandas as pd
from gramep.messages import Messages
from scipy.cluster.hierarchy import ClusterNode

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

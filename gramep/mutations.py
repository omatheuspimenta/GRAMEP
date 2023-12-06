"""mutations module.

This module contains the functions used to get the mutates from sequences using the
maximum entropy principle.

Contents:
    * get_mutations: Perform mutation analysis on sequence data.
    * get_variants_intersection: Get the intersection of variants.

Todo:
    * Implement tests.
"""
from collections import defaultdict
from sys import exit

from gramep.analysis import (  # , intersection_sequences
    kmers_analysis,
    variants_analysis,
)
from gramep.data_io import (
    annotation_dataframe,
    buffer_sequences,
    load_exclusive_kmers_file,
    load_sequences,
    reference_sequence,
    save_diffs_positions,
    save_exclusive_kmers,
    save_intersection_kmers,
    write_frequencies,
    write_report,
)
from gramep.graphics import plot_graphic
from gramep.helpers import get_kmers_frequencies, split_sequence
from gramep.kmers_utils import kmers_difference, kmers_intersections
from gramep.messages import Messages

message = Messages()
"""
Set the Message class for logging.
"""


def get_mutations(
    reference_path: str,
    sequence_path: str,
    save_path: str,
    word: int,
    step: int,
    annotation_path: str | None = None,
    snps_max: int = 1,
    dictonary: str = 'ACTG',
    create_report: bool = False,
    save_kmers: bool = False,
    load_exclusive_kmers: bool = False,
    path_exclusive_kmers: str | None = None,
    chunk_size: int = 100,
):
    """
    Perform mutation analysis on sequence data.

    This function performs mutation analysis on the provided sequence data.\
    It calculates variations, k-mer frequencies, and other relevant\
    information based on the input parameters.

    Args:
        reference_path (str): The path to the reference sequence data file.
        sequence_path (str): The path to the sequence data file.
        annotation_path (str): The path to the annotation data file.
        save_path (str): The path to save the generated results.
        word (int): The length of each k-mer.
        step (int): The step size for moving the sliding window.
        snps_max (int, optional): The maximum number of allowed SNPs within \
        an exclusive k-mer. Default is 1.
        dictonary (str, optional): The DNA dictionary for k-mer analysis. \
        Default is 'ACTG'.
        create_report (bool, optional): Whether to create a report. Default is False.
        save_kmers (bool, optional): Whether to save k-mers to a file. Default is False.
        load_exclusive_kmers (bool, optional): Whether to load exclusive k-mers \
        from a file. Default is False.
        path_exclusive_kmers (str | None, optional): The path to the exclusive \
        k-mers file. Default is None.
        chunk_size (int, optional): The chunk size for loading sequences. \
        Default is 100.

    Returns:
        Message class: A message confirming the analysis was completed.
    """

    message.info_start_objetive('get-mutations method')
    # Load reference sequence
    ref_sequence = reference_sequence(seq_path=reference_path)
    ref_sequence_split = split_sequence(
        sequence=ref_sequence, word=word, step=step
    )
    # Buffer sequences
    message.info_loading_kmers()
    sequences = buffer_sequences(sequence_path=sequence_path)
    reference = buffer_sequences(sequence_path=reference_path, reference=True)

    # Check if report will be generated
    annotation_df, sequence_interval = None, None
    if create_report:
        if annotation_path is not None:
            annotation_df, sequence_interval = annotation_dataframe(
                annotation_path=annotation_path
            )
        else:
            message.warning_annotation_file()
            annotation_df, sequence_interval = None, None

    # Load exclusive kmers
    if load_exclusive_kmers:
        seq_kmers_exclusive = load_exclusive_kmers_file(
            path_exclusive_kmers=path_exclusive_kmers
        )
        message.info_kmers_load()
    else:
        seq_kmers = load_sequences(
            seq_records=sequences,
            word=word,
            step=step,
            dictonary=dictonary,
            reference=False,
            chunk_size=chunk_size,
        )
        ref_kmers = load_sequences(
            seq_records=reference,
            word=word,
            step=step,
            dictonary=dictonary,
            reference=True,
            chunk_size=chunk_size,
        )

        seq_kmers_exclusive = kmers_difference(seq_kmers, ref_kmers)
        seq_kmers_intersections = kmers_intersections(seq_kmers, ref_kmers)

        if save_kmers:
            save_exclusive_kmers(
                sequence_path=sequence_path,
                seq_kmers_exclusive=seq_kmers_exclusive,
                save_path=save_path,
            )
            save_intersection_kmers(
                sequence_path=sequence_path,
                seq_kmers_intersections=seq_kmers_intersections,
                save_path=save_path,
            )

    # Analize kmers
    message.info_get_kmers()
    diffs_positions, report = kmers_analysis(
        seq_list=sequences,
        word=word,
        step=step,
        snps_max=snps_max,
        ref_sequence_split=ref_sequence_split,
        annotation_dataframe=annotation_df,
        seq_kmers_exclusive=seq_kmers_exclusive,
        sequence_interval=sequence_interval,
        create_report=create_report,
        chunk_size=chunk_size,
    )
    if diffs_positions is None:
        variations = []
        save_diffs_positions(
            sequence_path=sequence_path,
            ref_sequence=ref_sequence,
            variations=variations,
            save_path=save_path,
        )
        write_report(
            report=[], sequence_path=sequence_path, save_path=save_path
        )
        message.error_no_exclusive_kmers()
        exit(1)

    freq_kmers, variations = get_kmers_frequencies(diffs_positions)

    if save_kmers:
        save_diffs_positions(
            sequence_path=sequence_path,
            ref_sequence=ref_sequence,
            variations=variations,
            save_path=save_path,
        )

    if create_report:
        write_report(
            report=report, sequence_path=sequence_path, save_path=save_path
        )

    write_frequencies(
        freq_kmers=freq_kmers, sequence_path=sequence_path, save_path=save_path
    )

    plot_graphic(
        variations=variations,
        ref_sequence=ref_sequence,
        freq_kmers=freq_kmers,
        sequence_name=sequence_path,
        save_path=save_path,
    )
    return message.info_done()


def get_variants_intersection(
    save_path: str, intersection_seletion: str = 'ALL'
) -> defaultdict[str, list[str]]:
    """
    Get variants intersection based on selection criteria.

    This function retrieves variants intersection data based on the specified \
    selection criteria.
    The function reads data from the provided save path and performs intersection \
    calculations
    according to the chosen selection option.

    Args:
        save_path (str): The path to the directory containing data to process.
        intersection_seletion (str, optional): The selection criteria for variants \
        intersection. Options: 'ALL' (default).

    Returns:
        defaultdict[str, list[str]]: A dictionary mapping sequence IDs to lists of \
        variants based on the specified selection criteria.
    """

    message.info_start_objetive('get_variants_intersection method')
    variants_intersections = variants_analysis(
        save_path, intersection_seletion
    )
    message.info_done()
    return variants_intersections


# def get_sequence_intersection(reference_path: str,
#                               save_path: str,
#                               word: int,
#                               step: int,
#                               intersection_seletion: str = 'ALL'):
#     """

#     """
#     message.info_start_objetive('get_sequence_intersection method')

#     intersection_kmers = intersection_sequences(save_path, word, step, intersection_seletion)

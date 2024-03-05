"""
analysis
========================

This module provides functions for analyzing exclusive k-mers obtained using the \
    maximum entropy principle.

Contents:
    * run_exclusive_kmers: Run analysis of exclusive k-mers using sequence data.
    * kmers_analysis: Perform k-mers analysis and optionally generate a report.
    * variants_analysis: Perform variants analysis based on intersection selection.


Todo:
    * Implement tests.
"""
from collections import defaultdict
from itertools import combinations, count, groupby

import numpy as np
import pandas as pd
from Bio import File, SeqRecord
from joblib import Parallel, delayed
from joblib_progress import joblib_progress
from Levenshtein import distance
from matplotlib import pyplot as plt
from rich.progress import (
    BarColumn,
    MofNCompleteColumn,
    Progress,
    SpinnerColumn,
    TaskProgressColumn,
    TextColumn,
    TimeElapsedColumn,
    TimeRemainingColumn,
)
from upsetplot import from_contents, plot

from gramep.data_io import load_variants_exclusive  # , load_intersection_kmers
from gramep.helpers import (
    batched,
    check_diffs,
    dicts_append,
    dicts_drop_duplicates,
    make_report,
    split_sequence,
)
from gramep.messages import Messages

message = Messages()
"""
Set the Message class for logging.
"""


def run_exclusive_kmers(
    sequence: SeqRecord.SeqRecord,
    word: int,
    step: int,
    snps_max: int,
    ref_sequence_split: list[str],
    seq_kmers_exclusive: list[str],
) -> defaultdict[str, list[str]]:
    """
    Run analysis of exclusive k-mers using sequence data.

    This function takes sequence data, k-mer parameters, SNP maximum count, reference
    sequence split, and exclusive k-mers. It analyzes the exclusive k-mer
    in the context of the provided data and parameters, and returns a
    defaultdict mapping exclusive k-mers to lists of k-mers.

    Args:
        sequence (SeqRecord.SeqRecord): A SeqRecord.SeqRecord object \
        containing the sequence.
        word (int): The length of each k-mer.
        step (int): The step size for moving the sliding window.
        snps_max (int): The maximum number of SNPs allowed.
        ref_sequence_split (list[str]): A list of splitted reference sequence.
        seq_kmers_exclusive (list[str]): A list of exclusive k-mers.

    Returns:
        defaultdict[str, list[str]]: A defaultdict mapping exclusive k-mers to lists
                                     of k-mers that contains information about the \
                                     localization and mutations found.
    """

    variation_position: defaultdict[str, list[str]] = defaultdict(list[str])
    ref_seq = np.array(ref_sequence_split)

    sequence_str = str(sequence.seq).upper()
    sequence_name = str(sequence.name)
    variation_position[sequence_name] = []
    seq = np.array(split_sequence(sequence_str, word, step))

    var_kmers = [
        str(value) for value in np.intersect1d(seq_kmers_exclusive, seq)
    ]

    for var_adj in var_kmers:
        levenshtein_distance = np.array(
            [distance(var_adj, ref_adj) for ref_adj in ref_seq]
        )
        lev_dist_snps = np.array(
            np.where(levenshtein_distance <= snps_max)[0].tolist()
        )
        lev_dist_ref_index = lev_dist_snps * step
        lev_dist_ind_group = [
            min(list(g))
            for _, g in groupby(
                lev_dist_snps, lambda n, c=count(): n - next(c)
            )
        ]

        var_ref_kmers = [
            (str(ref_seq[ind_ref_split]), var_adj, ind_ref)
            for (ind_ref_split, ind_ref) in zip(
                lev_dist_ind_group, lev_dist_ref_index
            )
        ]
        if len(var_ref_kmers) != 0:
            variation_position[sequence_name].extend(
                np.concatenate(
                    [
                        check_diffs(*var_ref_adj)
                        for var_ref_adj in var_ref_kmers
                    ]
                )
            )

    variation_position[sequence_name] = np.unique(
        variation_position[sequence_name]
    ).tolist()

    return variation_position


def kmers_analysis(
    seq_list: File._IndexedSeqFileDict,
    word: int,
    step: int,
    snps_max: int,
    ref_sequence_split: list[str],
    annotation_dataframe: pd.DataFrame,
    seq_kmers_exclusive: list[str],
    sequence_interval: pd.Series,
    create_report: bool = False,
    chunk_size: int = 100,
) -> tuple[defaultdict[str, list[str]], np.ndarray] | tuple[
    None, None
] | tuple[defaultdict[str, list[str]], None]:
    """
    Perform k-mers analysis and optionally generate a report.

    This function performs k-mers analysis on the provided sequence data, exclusive \
        k-mers, and annotations. It calculates exclusive adjacencies, checks \
        differences, and returns results in a tuple. If 'create_report' is \
        set to True, a report is generated.

    Args:
        seq_list (File._IndexedSeqFileDict): A File._IndexedSeqFileDict object \
        containing the sequences in FASTA format.
        word (int): The length of each k-mer.
        step (int): The step size for moving the sliding window.
        snps_max (int): The maximum number of SNPs allowed within an exclusive \
            adjacency.
        ref_sequence_split (list[str]): A list of reference sequence k-mers splitted.
        annotation_dataframe (pd.DataFrame): DataFrame containing sequence annotations.
        seq_kmers_exclusive (list[str]): A list of exclusive k-mers.
        sequence_interval (pd.Series): Series containing sequence intervals.
        create_report (bool, optional): Whether to generate a report. Default is False.
        chunk_size (int, optional): The chunk size for loading sequences. \
        Default is 100.

    Returns:
        tuple[defaultdict[str, list[str]], np.ndarray]: A tuple containing results \
            of k-mers analysis and optionally a generated report.
    """

    progress = Progress(
        SpinnerColumn(spinner_name='dots'),
        TaskProgressColumn(),
        TextColumn('[progress.description]{task.description}'),
        BarColumn(),
        MofNCompleteColumn(),
        TimeElapsedColumn(),
        '<',
        TimeRemainingColumn(),
    )

    if len(seq_kmers_exclusive) == 0:
        message.error_no_exclusive_kmers()
        return None, None

    diffs_positions_temp = defaultdict(list)
    seq_list_generator = (name for name in seq_list.keys())
    with progress:
        task = progress.add_task(
            '[cyan]Getting SNPs positions ...', total=len(seq_list)
        )
        for batch_seqs in batched(seq_list_generator, chunk_size):
            diffs_positions_list = Parallel(n_jobs=-2)(
                delayed(run_exclusive_kmers)(
                    seq_list[seq_name],
                    word,
                    step,
                    snps_max,
                    ref_sequence_split,
                    seq_kmers_exclusive,
                )
                for seq_name in batch_seqs
            )

            diffs_positions_append = dicts_append(*diffs_positions_list)
            del diffs_positions_list
            diffs_positions_rmdup = dicts_drop_duplicates(
                diffs_positions_append
            )
            for key, value in diffs_positions_rmdup.items():
                diffs_positions_temp[key].extend(value)
            del diffs_positions_append, diffs_positions_rmdup
            progress.update(task, advance=len(batch_seqs))

    diffs_positions = dicts_drop_duplicates(diffs_positions_temp)
    if create_report:
        with joblib_progress(
            'Creating report ...', total=len(diffs_positions)
        ):
            report_list = Parallel(n_jobs=-2)(
                delayed(make_report)(
                    diffs_positions[key],
                    key,
                    sequence_interval,
                    annotation_dataframe,
                )
                for key in diffs_positions.keys()
            )
        report = np.hstack(report_list)

        return diffs_positions, report
    else:
        return diffs_positions, None


def variants_analysis(
    save_path: str, intersection_seletion: str = 'ALL'
) -> defaultdict[str, list[str]]:
    """
    Perform variants analysis based on intersection selection.

    This function performs variants analysis based on the specified intersection \
        selection criteria.
    It reads variant data from the provided file and returns a defaultdict \
        containing analysis results.

    Args:
        save_path (str): The path to the file containing variant data.
        intersection_seletion (str, optional): Criteria for selecting which variants \
        to intersect. To specify the variants for intersection, provide them \
        separated by '-'. For example: 'variant1-variant2-variant3'. Default is 'ALL'.

    Returns:
        defaultdict[str, list[str]]: A defaultdict containing analysis results.
    """
    variants_exclusive_kmers, variants_names = load_variants_exclusive(
        save_path
    )
    intersection_kmers = defaultdict(str)
    intersection_kmers_sets = defaultdict(str)

    if intersection_seletion == 'ALL':
        for r in range(len(variants_names) + 1):
            for subset in combinations(variants_names, r):
                if len(subset) > 1:
                    intersection_selects = list(subset)
                    intersection_set = list(
                        set.intersection(
                            *(
                                set(variants_exclusive_kmers[k])
                                for k in intersection_selects
                            )
                        )
                    )
                    if len(intersection_set) > 0:
                        intersection_kmers[
                            '-'.join(intersection_selects)
                        ] = intersection_set
                        for variant in intersection_selects:
                            intersection_kmers_sets[variant] = intersection_set
    else:
        intersection_selects = intersection_seletion.split(sep='-')
        intersection_kmers[intersection_seletion] = list(
            set.intersection(
                *(
                    set(variants_exclusive_kmers[k])
                    for k in intersection_selects
                )
            )
        )
        for variant in intersection_selects:
            intersection_kmers_sets[variant] = intersection_kmers[
                intersection_seletion
            ]

    with open(save_path + '/intersections.txt', 'a') as export_file:
        export_file.write('INTERSECTIONS\n')
        for k, v in intersection_kmers.items():
            export_file.write(str(k) + ': ' + str(', '.join(v)) + '\n')

    plot(from_contents(intersection_kmers_sets))
    plt.savefig(save_path + '/intersections.png')

    message.info_intersections_saved(save_path + '/intersections.txt')
    return intersection_kmers


# def intersection_sequences(save_path: str,
#                               word: int,
#                               step: int,
#                               intersection_seletion):

#     exclusive_kmers, variants_names = load_intersection_kmers(save_path=save_path)

#     intersection_kmers = defaultdict(str)

#     if intersection_seletion == 'ALL':
#         for r in range(len(variants_names) + 1):
#             for subset in combinations(variants_names, r):
#                 if len(subset) > 1:
#                     intersection_selects = list(subset)
#                     intersection_set = list(
#                         set.intersection(
#                             *(
#                                 set(exclusive_kmers[k])
#                                 for k in intersection_selects
#                             )
#                         )
#                     )
#                     if len(intersection_set) > 0:
#                         intersection_kmers[
#                             '-'.join(intersection_selects)
#                         ] = intersection_set

#     else:
#         intersection_selects = intersection_seletion.split(sep='-')
#         intersection_kmers[intersection_seletion] = list(
#             set.intersection(
#                 *(
#                     set(exclusive_kmers[k])
#                     for k in intersection_selects
#                 )
#             )
#         )

#     return intersection_kmers

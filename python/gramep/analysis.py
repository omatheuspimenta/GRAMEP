"""
analysis
========================

This module provides functions for analyzing exclusive k-mers obtained using the \
    maximum entropy principle.

Contents:
    * mutations_analysis: Perform k-mers analysis and optionally generate a report.
    * variants_analysis: Perform variants analysis based on intersection selection.

Todo:
    * Implement tests.
"""
from collections import defaultdict

import numpy as np
import pandas as pd
from gramep.data_io import load_variants_exclusive
from gramep.helpers import make_report
from gramep.messages import Messages
from gramep.utilrs import kmers_analysis, variants_intersection
from joblib import Parallel, delayed
from joblib_progress import joblib_progress
from matplotlib import pyplot as plt
from rich.progress import (
    BarColumn,
    Progress,
    SpinnerColumn,
    TaskProgressColumn,
    TextColumn,
    TimeElapsedColumn,
)
from upsetplot import from_contents, plot

message = Messages()
"""
Set the Message class for logging.
"""


def mutations_analysis(
    seq_path: str,
    ref_path: str,
    seq_kmers_exclusive: list[str],
    word: int,
    step: int,
    snps_max: int,
    annotation_dataframe: pd.DataFrame,
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
        seq_path (str): The path to the file containing the sequences in FASTA format.
        ref_path (str): The path to the reference sequence data file.
        seq_kmers_exclusive (list[str]): A list of exclusive k-mers.
        word (int): The length of each k-mer.
        step (int): The step size for moving the sliding window.
        snps_max (int): The maximum number of SNPs allowed.
        annotation_dataframe (pd.DataFrame): DataFrame containing sequence annotations.
        sequence_interval (pd.Series): Series containing sequence intervals.
        create_report (bool, optional): Whether to generate a report. Default is False.
        chunk_size (int, optional): The chunk size for loading sequences. \
        Default is 100.

    Returns:
        tuple[defaultdict[str, list[str]], np.ndarray]: A tuple containing results \
            of k-mers analysis and optionally a generated report.
    """

    progress = Progress(
        SpinnerColumn(),
        TaskProgressColumn(),
        TextColumn('[progress.description]{task.description}'),
        BarColumn(),
        TimeElapsedColumn(),
    )

    if len(seq_kmers_exclusive) == 0:
        message.error_no_exclusive_kmers()
        return None, None

    with progress:
        progress.add_task('[cyan]Getting SNPs positions ...', total=None)

        diffs_positions = kmers_analysis(
            seq_path=seq_path,
            ref_path=ref_path,
            exclusive_kmers=seq_kmers_exclusive,
            k=word,
            step=step,
            max_dist=snps_max,
            batch_size=chunk_size,
        )

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

    intersection_kmers, intersection_kmers_sets = variants_intersection(
        variants_exclusive_kmers, variants_names, intersection_seletion
    )

    with open(save_path + '/intersections.txt', 'a') as export_file:
        export_file.write('INTERSECTIONS\n')
        for k, v in intersection_kmers.items():
            export_file.write(str(k) + ': ' + str(', '.join(v)) + '\n')

    plot(from_contents(intersection_kmers_sets))
    plt.savefig(save_path + '/intersections.png')

    message.info_intersections_saved(save_path + '/intersections.txt')
    return intersection_kmers

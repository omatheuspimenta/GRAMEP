"""data_io module

This module provides functions for loading and saving data.

Contents:
    * annotation_dataframe: Process annotation information from a GFF3 file into \
    a DataFrame and a sequence interval Series.
    * load_sequences: Load and select the most informative kmers from sequences \
    * save_exclusive_kmers: Save exclusive k-mers to a file.
    * save_intersection_kmers: Save intersection k-mers to a file.
    * save_diffs_positions: Save variations and their positions to a file.
    * write_report: Write a report to a file.
    * write_frequencies: Write k-mer frequencies to a file.
    * load_variants_exclusive: Load variants and exclusive k-mers from a folder.
    * load_variants_kmers: Load and return unique exclusive k-mers from saved files.
    * save_data: Save data to files in the specified directory.
    * save_ranges: Save MinMaxScaler object ranges to a file.
    * load_ranges: Load MinMaxScaler object ranges from a file.
    * load_model: Load a RandomForestClassifier model from a file.
    * save_model: Save a RandomForestClassifier model to a file.
    * save_metrics: Save accuracy and metrics data to a file.
    * save_confusion_matrix: Save a confusion matrix plot to a file.
    * save_predict_data: Save predict data to a file.
    * load_exclusive_kmers_file: Load exclusive k-mers from a file.
    * load_mutations: Load and return the unique mutations.
    * load_reports: Load and return the reports dirnames.
    * load_seqs_phylogenic: Loads the sequences and returns the mutations in binary format.
    * load_dataframe_phylo: Load the dataframe from sequences and mutations.
    * save_newick: Save the newick string to a file.
    
Todo:
    * Implement tests.
"""
import pickle
from collections import defaultdict
from os import scandir
from pathlib import Path

import numpy as np
import numpy.typing as npt
import pandas as pd
import polars as pl
from gffpandas import gffpandas as gffpd
from gramep.entropy_utils import maxentropy
from gramep.helpers import get_sequence_interval
from gramep.kmers_utils import select_kmers
from gramep.messages import Messages
from gramep.utilrs import get_kmers, write_ref
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
from seaborn import heatmap
from sklearn.ensemble import RandomForestClassifier
from sklearn.preprocessing import MinMaxScaler

message = Messages()
"""
Set the Message class for logging.
"""


def annotation_dataframe(
    annotation_path: str,
) -> tuple[pd.DataFrame, pd.Series]:
    """
    Process annotation information from a GFF3 file into a DataFrame and a sequence \
    interval Series.

    Args:
        annotation_path (str): The path to the GFF3 annotation file.

    Returns:
        tuple: A tuple containing two elements:
            - df_annotation (pd.DataFrame): A DataFrame with processed annotation data.
            - seq_interval (pd.Series): A Series containing sequence intervals \
            generated from the data.
    """
    df_annotation = gffpd.read_gff3(annotation_path).attributes_to_columns()
    df_annotation = df_annotation[['type', 'start', 'end', 'Name', 'product']][
        df_annotation['Name'].notna()
    ]
    seq_interval = df_annotation.apply(get_sequence_interval, axis=1)
    return df_annotation, seq_interval


def load_sequences(
    file_path: str,
    word: int,
    step: int,
    dictonary: str = 'DNA',
    reference: bool = False,
    chunk_size: int = 100,
) -> defaultdict[str, int]:
    """
    Load and select the most informative kmers from sequences from a file into \
    overlapping k-mers.

    This function reads sequences from a file and select the most informative
    k-mers of the specified length and step size. Optionally, the function can
    use a custom dictionary of characters and return a list of k-mers or the
    reference sequence as a single string.

    Args:
        file_path (str): The path to the file containing sequences.
        word (int): The length of each k-mer.
        step (int): The step size for moving the sliding window.
        dictonary (str, optional): A string containing characters to consider \
        in sequences. Default is 'DNA' (DNA alphabet).
        reference (bool, optional): If True, return the reference sequence.\
        Default is False.
        chunk_size (int, optional): The chunk size for loading sequences. \
        Default is 100.

    Returns:
        defaultdict[str, int]: A defaultdict mapping most informative k-mers.
    """
    message.info_start()

    progress = Progress(
        SpinnerColumn(),
        TaskProgressColumn(),
        TextColumn('[progress.description]{task.description}'),
        BarColumn(),
        TimeElapsedColumn(),
    )

    if reference:
        loading_text = 'Loading reference sequence ...'
    else:
        loading_text = 'Loading sequences ...'

    with progress:
        progress.add_task(f'[cyan]{loading_text}', total=None)
        kmers = get_kmers(file_path, word, step, dictonary, chunk_size)

    # Check if sequence is reference
    if reference:
        return select_kmers(0, kmers)

    # Entropy
    message.info_entropy()
    _, _, sequenceFrequency, _ = maxentropy(kmers)

    return select_kmers(sequenceFrequency, kmers)


def save_exclusive_kmers(
    sequence_path: str,
    seq_kmers_exclusive: list[str],
    save_path: str,
):
    """
    Save exclusive k-mers to a file.

    This function takes a sequence data file path, a list of exclusive k-mers,
    and a save path. It writes the exclusive k-mers to the specified file.

    Args:
        sequence_path (str): The path to the sequence data file.
        seq_kmers_exclusive (list[str]): A list of exclusive k-mers.
        save_path (str): The path to save the exclusive k-mers.

    Returns:
        Message class
    """
    seq_name = sequence_path.split('/')[-1].split('.')[0]
    # Check if path exists
    dirPath = save_path + '/' + seq_name
    path = Path(dirPath)
    path.mkdir(mode=0o777, parents=True, exist_ok=True)
    save_file_name = dirPath + '/' + seq_name + '_ExclusiveKmers.sav'
    with open(save_file_name, 'wb') as exclusive_kmers_file:
        pickle.dump(seq_kmers_exclusive, exclusive_kmers_file)
    save_file_name = dirPath + '/' + seq_name + '_ExclusiveKmers.txt'
    with open(save_file_name, 'w') as exclusive_kmers_file:
        exclusive_kmers_file.write(str(seq_kmers_exclusive))
    return message.info_kmers_saved(dirPath)


def save_intersection_kmers(
    sequence_path: str,
    seq_kmers_intersections: list[str],
    save_path: str,
):
    """
    Save intersection k-mers to a file.

    This function takes a sequence data file path, a list of intersection k-mers,
    and a save path. It writes the intersection k-mers to the specified file.

    Args:
        sequence_path (str): The path to the sequence data file.
        seq_kmers_intersections (list[str]): A list of intersection k-mers.
        save_path (str): The path to save the intersection k-mers.

    Returns:
        Message class
    """
    seq_name = sequence_path.split('/')[-1].split('.')[0]
    # Check if path exists
    dirPath = save_path + '/' + seq_name
    path = Path(dirPath)
    path.mkdir(mode=0o777, parents=True, exist_ok=True)
    save_file_name = dirPath + '/' + seq_name + '_IntersectionKmers.sav'
    with open(save_file_name, 'wb') as exclusive_kmers_file:
        pickle.dump(seq_kmers_intersections, exclusive_kmers_file)
    save_file_name = dirPath + '/' + seq_name + '_IntersectionKmers.txt'
    with open(save_file_name, 'w') as exclusive_kmers_file:
        exclusive_kmers_file.write(str(seq_kmers_intersections))

    return message.info_kmers_saved(dirPath)


def save_diffs_positions(
    sequence_path: str,
    ref_path: str,
    variations: list[str],
    save_path: str,
):
    """
    Save variations and their positions to a file.

    This function takes a sequence data file path, a list of variations, and a \
    save path.
    It writes the variations and their positions to the specified file in \
        multiple formats.

    Args:
        sequence_path (str): The path to the sequence data file.
        ref_path (str): The reference sequence.
        variations (list[str]): A list of variations.
        save_path (str): The path to save the variations and their positions.

    Returns:
        Message class
    """
    seq_name = sequence_path.split('/')[-1].split('.')[0]

    # Check if path exists
    dirPath = save_path + '/' + seq_name
    path = Path(dirPath)
    path.mkdir(mode=0o777, parents=True, exist_ok=True)
    save_file_name = dirPath + '/' + seq_name + '_variations.sav'
    with open(save_file_name, 'wb') as exclusive_kmers_file:
        pickle.dump(variations, exclusive_kmers_file)
    save_file_name = dirPath + '/' + seq_name + '_variations.txt'
    with open(save_file_name, 'w') as exclusive_kmers_file:
        exclusive_kmers_file.write(str(variations))
    save_file_name = dirPath + '/' + seq_name + '_variations.bed3'
    with open(save_file_name, 'w') as exclusive_kmers_file:
        for variation in variations:
            position, mutation = variation.split(':')
            exclusive_kmers_file.write(
                str(
                    str(seq_name)
                    + '\t'
                    + str(position)
                    + '\t'
                    + str(int(position) + 1)
                    + '\n'
                )
            )

    save_file_name = dirPath + '/' + seq_name + '_reference.fasta'

    write_ref(
        ref_path=ref_path, variations=variations, save_path=save_file_name
    )

    return message.info_kmers_saved(dirPath)


def write_report(report: np.ndarray, sequence_path: str, save_path: str):
    """
    Write a report to a file.

    This function takes a report as a NumPy array, a sequence data file path,
    and a save path. It writes the report to the specified file.

    Args:
        report (np.ndarray): A NumPy array representing the report.
        sequence_path (str): The path to the sequence data file.
        save_path (str): The path to save the report.

    Returns:
        Message class
    """
    seq_name = sequence_path.split('/')[-1].split('.')[0]
    # Check if path exists
    dirPath = save_path + '/' + seq_name
    path = Path(dirPath)
    path.mkdir(mode=0o777, parents=True, exist_ok=True)
    save_file_name = dirPath + '/' + seq_name + '_report.csv'

    with open(save_file_name, 'a') as out_file:
        out_file.write(
            'sequence_id;annotation_name;start;end;type;modification_localization_in_reference;reference_kmer;exclusive_variant_kmer;reference_snp;variant_snp\n'
        )
    pd.Series(list(report)).to_csv(
        save_file_name, header=False, index=False, mode='a'
    )

    return message.info_report_saved(dirPath)


def write_frequencies(
    freq_kmers: defaultdict, sequence_path: str, save_path: str
):
    """
    Write k-mer frequencies to a file.

    This function takes a defaultdict containing k-mer frequencies, a sequence \
    data file path,
    and a save path. It writes the k-mer frequencies to the specified file.

    Args:
        freq_kmers (defaultdict): A defaultdict containing k-mer frequencies.
        sequence_path (str): The path to the sequence data file.
        save_path (str): The path to save the k-mer frequencies.

    Returns:
        Message class
    """

    seq_name = sequence_path.split('/')[-1].split('.')[0]
    # Check if path exists
    dirPath = save_path + '/' + seq_name
    path = Path(dirPath)
    path.mkdir(mode=0o777, parents=True, exist_ok=True)
    save_file_name = dirPath + '/' + seq_name + '_FreqExclusiveKmers.csv'

    with open(save_file_name, 'a') as out_file:
        out_file.write('position;reference_value;variant_value;frequency\n')
        for key, value in freq_kmers.items():
            position, var = key.split(':')
            out_file.write(
                str(
                    str(position)
                    + ';'
                    + str(var[0])
                    + ';'
                    + str(var[1])
                    + ';'
                    + str(value)
                    + '\n'
                )
            )

    return message.info_freq_saved(dirPath)


def load_variants_exclusive(
    save_path: str,
) -> tuple[defaultdict[str, list[str]], list[str]]:
    """
    Load variants and exclusive variants from a folder.

    This function reads and loads variants and exclusive variants from the specified file.
    It returns a tuple containing a defaultdict of variant information and a \
    list of exclusive variants.

    Args:
        save_path (str): The path to the file containing variants folders and \
        exclusive variants.

    Returns:
        tuple[defaultdict[str, list[str]], list[str]]: A tuple containing a \
        defaultdict of variant information and a list of exclusive variants.
    """
    subfolders = [f.path for f in scandir(save_path) if f.is_dir()]
    loads = defaultdict(list)
    variants_names = []
    for dirname in subfolders:
        with open(
            dirname + '/' + dirname.split(sep='/')[-1] + '_variations.sav',
            'rb',
        ) as f:
            loads[dirname.split(sep='/')[-1]] = pickle.load(f)
        variants_names.append(dirname.split(sep='/')[-1])
    return loads, variants_names


def load_variants_kmers(save_path: str) -> list[str]:
    """
    Load and return unique exclusive k-mers from saved files.

    This function reads exclusive k-mers data from saved files located within \
    subfolders of the specified
    'save_path'. It returns a NumPy array containing the unique set of exclusive \
    k-mers across all files.

    Args:
        save_path (str): The path to the directory containing subfolders with \
        saved exclusive k-mers files.

    Returns:
        list: A list containing the unique set of exclusive k-mers.
    """
    subfolders = [f.path for f in scandir(save_path) if f.is_dir()]
    loads = []
    for dirname in subfolders:
        with open(
            dirname + '/' + dirname.split(sep='/')[-1] + '_ExclusiveKmers.sav',
            'rb',
        ) as f:
            seq_kmers_exclusive = [''.join(item) for item in pickle.load(f)]
            loads.extend(seq_kmers_exclusive)

    return np.unique(loads).tolist()


def save_data(
    data_frame: pd.DataFrame, class_names_to_save: np.ndarray, dir_path: str
):
    """
    Save data to files in the specified directory.

    This function saves the provided data frame and class names array to \
    files in the specified
    directory. The data frame is saved in CSV format, and the class names array \
    is saved as a text file.

    Args:
        data_frame (pd.DataFrame): The pandas DataFrame to be saved.
        class_names_to_save (np.ndarray): The array of class names to be saved.
        dir_path (str): The path to the directory where files will be saved.

    Returns:
        Message class
    """

    path_dir = str(Path(dir_path).parent) + '/classify/results'
    path = Path(path_dir)
    path.mkdir(mode=0o777, parents=True, exist_ok=True)

    data_frame['CLASS'] = class_names_to_save
    data_frame.to_csv(path_dir + '/dataframe.csv', encoding='utf-8')
    return message.info_dataframe_saved(path_dir)


def save_ranges(ranges: MinMaxScaler, dir_path: str):
    """
    Save MinMaxScaler object ranges to a file.

    This function saves the MinMaxScaler object's data range information to a file
    in the specified directory. The saved data can be used for scaling operations later.

    Args:
        ranges (MinMaxScaler): The MinMaxScaler object containing data \
        range information.
        dir_path (str): The path to the directory where the data range \
        file will be saved.

    Returns:
        None
    """

    path_dir = str(Path(dir_path).parent) + '/classify/model'
    path = Path(path_dir)
    path.mkdir(mode=0o777, parents=True, exist_ok=True)

    file_name = path_dir + '/ranges.sav'
    with open(file_name, 'wb') as ranges_file:
        pickle.dump(ranges, ranges_file)
    return message.info_ranges_saved(path_dir)


def load_ranges(load_ranges_path: str) -> MinMaxScaler:
    """
    Load MinMaxScaler object ranges from a file.

    This function loads a MinMaxScaler object's data range information from a \
    file and returns the MinMaxScaler object.

    Args:
        load_ranges_path (str): The path to the file containing MinMaxScaler \
        object ranges.

    Returns:
        MinMaxScaler: A MinMaxScaler object with loaded data range information.
    """

    with open(load_ranges_path, 'rb') as handle:
        return pickle.load(handle)


def load_model(load_model_path: str) -> RandomForestClassifier:
    """
    Load a RandomForestClassifier model from a file.

    This function loads a trained RandomForestClassifier model from a file and \
    returns the model.

    Args:
        load_model_path (str): The path to the file containing the trained model.

    Returns:
        RandomForestClassifier: A trained RandomForestClassifier model.
    """
    with open(load_model_path, 'rb') as handle:
        return pickle.load(handle)


def save_model(model: RandomForestClassifier, dir_path: str):
    """
    Save a RandomForestClassifier model to a file.

    This function saves the provided RandomForestClassifier model to a \
    file in the specified
    directory. The saved model can be loaded and used for predictions later.

    Args:
        model (RandomForestClassifier): The RandomForestClassifier model to be saved.
        dir_path (str): The path to the directory where the model file will be saved.

    Returns:
        Message class
    """

    path_dir = str(Path(dir_path).parent) + '/classify/model'
    path = Path(path_dir)
    path.mkdir(mode=0o777, parents=True, exist_ok=True)

    file_name = path_dir + '/model.sav'
    with open(file_name, 'wb') as model_file:
        pickle.dump(model, model_file)
    return message.info_model_saved(path_dir)


def save_metrics(acc: str, metrics: str, dir_path: str):
    """
    Save accuracy and metrics data to a file.

    This function saves the provided accuracy and metrics data to a file in \
    the specified directory.
    The saved data can be used for analysis and reporting.

    Args:
        acc (str): The accuracy data to be saved.
        metrics (str): The metrics data to be saved.
        dir_path (str): The path to the directory where the data file will be saved.

    Returns:
        Message class
    """

    path_dir = str(Path(dir_path).parent) + '/classify/results'
    path = Path(path_dir)
    path.mkdir(mode=0o777, parents=True, exist_ok=True)

    file_name = path_dir + '/metrics.txt'
    with open(file_name, 'a') as metrics_file:
        metrics_file.write(
            'GRAMEP - Genome vaRiation Analysis from the Maximum Entropy Principle.\n'
        )
        metrics_file.write('Results from validation\n')
        metrics_file.write('Accuracy:')
        metrics_file.write(acc)
        metrics_file.write('\nMetrics:\n')
        metrics_file.write(metrics)
    return message.info_metrics_saved(path_dir)


def save_confusion_matrix(
    conf_mtx: np.ndarray,
    name_class: np.ndarray,
    vmax: int,
    dir_path: str,
    is_conf_mtx: bool = True,
):
    """
    Save a confusion matrix plot to a file.

    This function saves a confusion matrix plot, based on the provided \
    confusion matrix and class names, to a file in the specified directory. \
    The plot includes the color scale, with the maximum value set \
    by the 'vmax' parameter.

    Args:
        conf_mtx (np.ndarray): The confusion matrix data to be plotted.
        name_class (np.ndarray): The array of class names corresponding to the matrix.
        vmax (int): The maximum value for the color scale in the plot.
        dir_path (str): The path to the directory where the plot file will be saved.
        is_conf_mtx (bool): Boolean. If True, print the Confusion Matrix title. Else, distance text.

    Returns:
        Message class
    """

    # if len(conf_mtx) > 200:
    #     figsize = (int(len(conf_mtx) / 4), int(len(conf_mtx) / 4))
    if is_conf_mtx:
        path_dir = str(Path(dir_path).parent) + '/classify/results'
        path = Path(path_dir)
        path.mkdir(mode=0o777, parents=True, exist_ok=True)
        text = 'Heatmap for Random Forest\nClassification Algorithm'
        save_name = '/confusion_matrix.pdf'
        fmt = 'd'
        figsize = (20, 15)
    else:
        text = 'Heatmap for Distance Matrix'
        save_name = '/heatmap.pdf'
        fmt = '.2f'
        path_dir = str(Path(dir_path).parent) + '/phylogenics/'
        path = Path(path_dir)
        path.mkdir(mode=0o777, parents=True, exist_ok=True)
        figsize = (25, 25)

    cm_fig = plt.figure(figsize=figsize)
    axes = plt.axes()
    # x_axis_labels = name_class
    # y_axis_labels = name_class
    heatmap(
        conf_mtx,
        vmin=0,
        vmax=vmax,
        annot=True,
        fmt=fmt,
        ax=axes,
        xticklabels=name_class,
        yticklabels=name_class,
    )
    axes.set_title(text, fontsize=20, pad=15)
    cm_fig.savefig(path_dir + save_name, dpi=300)
    if is_conf_mtx:
        return message.result_confusion_matrix(path_dir)
    return message.result_heatmap(path_dir)


def save_predict_data(
    id_values: pd.Series, predicted_data: pd.Series, dir_path: str
) -> str:
    """
    Save predicted data to a CSV file.

    This function saves the provided 'ID' values and predicted data as a CSV \
    file in the specified
    directory. It returns a message confirming the successful saving of the predictions.

    Args:
        id_values (pd.Series): A pandas Series containing 'ID' values.
        predicted_data (pd.Series): A pandas Series containing predicted data.
        dir_path (str): The path to the directory for saving the CSV file.

    Returns:
        str: A message confirming the successful saving of the predictions.
    """

    pd.DataFrame({'ID': id_values, 'Predicted': predicted_data}).to_csv(
        dir_path + '/predict_data.csv', encoding='utf-8', index=False
    )
    message.info_done()
    return message.info_predictions_saved(dir_path)


def load_exclusive_kmers_file(path_exclusive_kmers: str) -> list[str]:
    """
    Load exclusive k-mers from a file.

    This function loads a list of exclusive k-mers from a file located at \
    the specified path.
    The function supports loading k-mers from both plain text files \
    and binary (pickle) files.

    Args:
        path_exclusive_kmers (str): The path to the file containing exclusive k-mers.

    Returns:
        list[str]: A list of exclusive k-mers loaded from the file.

    Raises:
        ValueError: If the k-mers have different lengths in the file.
    """
    if path_exclusive_kmers.split('.')[-1] != 'sav':
        with open(path_exclusive_kmers, 'r') as file_exclusive_kmers:
            seq_kmers_exclusive = file_exclusive_kmers.read().splitlines()
        it = iter(seq_kmers_exclusive)
        the_len = len(next(it))
        if not all(len(kmer_len) == the_len for kmer_len in it):
            return message.error_different_kmers_len
        else:
            return seq_kmers_exclusive
    else:
        with open(path_exclusive_kmers, 'rb') as file_exclusive_kmers:
            return pickle.load(file_exclusive_kmers)


def load_mutations(save_path: str) -> npt.NDArray[np.str_]:
    """
    Load and return the unique mutations.

    This function loads the mutations from the specified directory and
    returns them as a NumPy array.

    Args:
        save_path (str): The path to the directory containing subfolders with \
        saved mutations .sav files.

    Returns:
        np.ndarray: A NumPy array containing the unique set of mutations.
    """
    subfolders = [f.path for f in scandir(save_path) if f.is_dir()]
    loads = []
    for dirname in subfolders:
        with open(
            dirname + '/' + dirname.split(sep='/')[-1] + '_variations.sav',
            'rb',
        ) as f:
            mutations = [''.join(item) for item in pickle.load(f)]
            loads.extend(mutations)

    return np.unique(loads)


def load_reports(save_path: str) -> list[str]:
    """
    Load and return the reports dirnames.

    This function loads the reports directory and
    returns them as a list of strings.

    Args:
        save_path (str): The path to the directory containing subfolders with \
        reports in _report.csv format.

    Returns:
        list[str]: A list of strings containing the reports dirnames.
    """
    subfolders = [f.path for f in scandir(save_path) if f.is_dir()]
    loads_names = [
        dirname + '/' + dirname.split(sep='/')[-1] + '_report.csv'
        for dirname in subfolders
    ]
    return loads_names


def load_seqs_phylogenic(
    loadname: str, variants_mutations: np.ndarray
) -> defaultdict[str, int]:
    """
    Loads the sequences and returns the mutations in binary format.

    This function loads the sequences from the report in specified directory and
    returns them as a dictionary of mutations in binary format.

    Args:
        loadname (str): The path to the report file.
        variants_mutations (np.ndarray): The unique set of mutations.

    Returns:
        defaultdict[str, int]: A dictionary of mutations in binary format.
    """

    report = pl.scan_csv(
        loadname,
        separator=';',
        infer_schema_length=0,
    )
    report = report.drop_nulls()
    report_mutations = (
        report.lazy()
        .with_columns(
            (
                pl.col('modification_localization_in_reference')
                + ':'
                + pl.col('reference_snp')
                + pl.col('variant_snp')
            ).alias('mutations')
        )
        .collect()
        .group_by('sequence_id')
        .agg(pl.col('mutations'))
    )
    report_group = dict(report_mutations.iter_rows())
    ret_sample = []
    for k, v in report_group.items():
        sample = defaultdict(int, {k: 0 for k in variants_mutations})
        for item in np.unique(v):
            sample[item] += 1
        sample['ID'] = k + '_' + loadname.split('/')[-2]
        ret_sample.append(sample)

    return ret_sample


def load_dataframe_phylo(save_path: str) -> pl.DataFrame:
    """
    Load the dataframe from sequences and mutations.

    This function loads the sequences from the report in specified directory and
    returns them as a dataframe, where the columns are the mutations and the rows
    are the sequences.

    Args:
        save_path (str): The path to the directory containing subfolders with \
        reports in _report.csv format.

    Returns:
        pl.DataFrame: A dataframe containing the mutations in binary format.
    """
    variants_mutations = load_mutations(save_path=save_path)
    loads_names = load_reports(save_path=save_path)

    with joblib_progress(
        '[cyan]Loading reports files ...', total=len(list(loads_names))
    ):
        feat_list = Parallel(n_jobs=-2)(
            delayed(load_seqs_phylogenic)(filename, variants_mutations)
            for filename in loads_names
        )

    return pl.DataFrame(list(np.hstack(feat_list))).lazy()


def save_newick(newick_str: str, save_path: str):
    """
    Save a Newick tree string to a file.

    This function saves a Newick tree string to a file.

    Args:
        newick_str (str): The Newick tree string.
        save_path (str): The path to the file to save the Newick tree string to.
    """
    path_dir = str(Path(save_path).parent) + '/phylogenics/'
    path = Path(path_dir)
    path.mkdir(mode=0o777, parents=True, exist_ok=True)
    path_save = path_dir + '/newick.tree'
    with open(path_save, 'w') as f:
        f.write(newick_str)

    return message.info_newick_saved(path_dir)

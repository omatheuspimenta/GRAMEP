"""
classify_utils
===================

This module provides functions for classifying biological sequences into variants.

Contents:
    * extract_features: Extract features from sequences and save to a DataFrame.
    * process_dataframe: Process a DataFrame and optionally save data and model.
    * sequence_classification: Perform sequence classification based on provided \
    data and options.
    * classify: Perform sequence classification pipeline.
    * extract_features_to_predict: Extract features from sequences for prediction \
    and return as a DataFrame.
    * process_dataframe_predict: Process a DataFrame for prediction using MinMaxScaler.
    * predict_data: Predict classes using a trained RandomForestClassifier model.
    * predict: Predict sequence classes using a trained model.

Note:
    This module is designed to work with biological sequences and their classifications,
    allowing researchers to quickly classify and analyze sequence variants.

Todo:
    * Implement tests.
"""

from fnmatch import fnmatch
from os import listdir

import numpy as np
import pandas as pd
from gramep.data_io import (
    load_model,
    load_ranges,
    load_variants_kmers,
    save_confusion_matrix,
    save_data,
    save_metrics,
    save_model,
    save_predict_data,
    save_ranges,
)
from gramep.messages import Messages
from gramep.utilrs import load_sequences_classify
from rich.progress import (
    BarColumn,
    Progress,
    SpinnerColumn,
    TaskProgressColumn,
    TextColumn,
    TimeElapsedColumn,
)
from sklearn.ensemble import RandomForestClassifier
from sklearn.metrics import (
    accuracy_score,
    classification_report,
    confusion_matrix,
)
from sklearn.model_selection import (
    RepeatedStratifiedKFold,
    cross_validate,
    train_test_split,
)
from sklearn.preprocessing import LabelEncoder, MinMaxScaler

message = Messages()
"""
Set the Message class for logging.
"""


def extract_features(
    word: int,
    step: int,
    save_path: str,
    dir_path: str,
    dictonary: str,
    chunk_size: int = 100,
) -> pd.DataFrame:
    """
    Extract features from sequences and save to a DataFrame.

    This function extracts features from sequences located in the specified \
    directory, and
    then saves the extracted features to a pandas DataFrame. The extracted \
    features are based
    on the exclusive kmers.

    Args:
        word (int): The length of each k-mer.
        step (int): The step size for moving the sliding window.
        save_path (str): The path to save the extracted features DataFrame.
        dir_path (str): The path to the directory containing sequence data.
        dictonary (str): The DNA dictionary for k-mer analysis.
        chunk_size (int, optional): The chunk size for loading sequences. \
        Default is 100.

    Returns:
        pd.DataFrame: A pandas DataFrame containing the extracted features.
    """

    message.info_start_objetive('Extracting the features ...')

    progress = Progress(
        SpinnerColumn(),
        TaskProgressColumn(),
        TextColumn('[progress.description]{task.description}'),
        BarColumn(),
        TimeElapsedColumn(),
    )

    variants_kmers = load_variants_kmers(save_path=save_path)

    file_list = [
        dir_path + '/' + name
        for name in listdir(dir_path)
        if fnmatch(name, '*.fasta')
    ]

    with progress:
        progress.add_task('[cyan]Loading sequences ...', total=None)
        feat_list = load_sequences_classify(
            file_list, word, step, dictonary, variants_kmers, False, chunk_size
        )

    data_frame = pd.DataFrame(feat_list)

    return data_frame


def process_dataframe(
    data_frame: pd.DataFrame,
    dir_path: str = None,
    should_save_data: bool = False,
    should_save_model: bool = False,
) -> tuple[pd.DataFrame, np.ndarray]:
    """
    Process a DataFrame and optionally save data and model.

    This function takes a DataFrame and performs processing on it. It optionally saves
    the processed data and a model based on the specified flags. It returns a tuple
    containing the processed DataFrame and a numpy array.

    Args:
        data_frame (pd.DataFrame): The DataFrame to be processed.
        dir_path (str, optional): The directory path for saving data and model. \
        Default is None.
        should_save_data (bool, optional): Whether to save processed data. \
        Default is False.
        should_save_model (bool, optional): Whether to save a model. Default is False.

    Returns:
        tuple[pd.DataFrame, np.ndarray]: A tuple containing the processed \
        DataFrame and a numpy array.
    """

    message.info_processing_dataframe()

    class_values = data_frame['CLASS']
    class_names_to_save = data_frame['CLASS']
    name_class = np.unique(class_values).tolist()

    data_frame.drop(columns=['CLASS'], axis=1, inplace=True)
    data_frame.replace([np.inf, -np.inf], 0, inplace=True)
    data_frame.replace(np.nan, 0, inplace=True)
    df_col_names = data_frame.columns

    # MinMax Scaler
    minMax_scaler = MinMaxScaler()
    minMax_scaler.fit(data_frame)
    df_minmax = minMax_scaler.transform(data_frame)
    data_frame = pd.DataFrame(df_minmax)
    del df_minmax
    label_encdr = LabelEncoder()
    class_values = label_encdr.fit_transform(class_values)

    data_frame.columns = df_col_names
    data_frame['CLASS'] = class_values

    if should_save_data:
        save_data(
            data_frame=data_frame,
            class_names_to_save=class_names_to_save,
            dir_path=dir_path,
        )
    if should_save_model:
        save_ranges(ranges=minMax_scaler, dir_path=dir_path)
    return data_frame, name_class


def sequence_classification(
    data_frame: pd.DataFrame,
    name_class: np.ndarray,
    dir_path: str,
    should_save_model: bool = False,
    should_save_confusion_matrix: bool = False,
) -> None:
    """
    Perform sequence classification based on provided data and options.

    This function performs sequence classification using the provided data frame \
    and class names.
    It allows for optional saving of a trained model and confusion matrix plot \
    based on specified flags.

    Args:
        data_frame (pd.DataFrame): The data frame containing sequence data and features.
        name_class (np.ndarray): The array of class names corresponding to the data.
        dir_path (str): The path to the directory for saving model and plot files.
        should_save_model (bool, optional): Whether to save the trained model. \
        Default is False.
        should_save_confusion_matrix (bool, optional): Whether to save the \
        confusion matrix plot. Default is False.

    Returns:
        None
    """

    message.info_classifying()

    # Split the data into training and test sets
    x_axis = data_frame.drop(columns='CLASS', axis=1)
    y_axis = data_frame['CLASS']
    x_train, x_test, y_train, y_test = train_test_split(
        x_axis, y_axis, test_size=0.2, stratify=y_axis
    )

    # Create a random forest classifier with default parameters
    rf_classifier = RandomForestClassifier(n_estimators=100)
    rf_classifier.fit(x_train, y_train)

    if should_save_model:
        save_model(model=rf_classifier, dir_path=dir_path)

    # Make predictions on the test set
    y_pred = rf_classifier.predict(x_test)

    # Evaluate the model - 10-fold
    cross_val = RepeatedStratifiedKFold(
        n_splits=10, n_repeats=10, random_state=7
    )
    n_scores = cross_validate(
        rf_classifier,
        x_axis,
        y_axis,
        scoring='accuracy',
        cv=cross_val,
        n_jobs=-1,
        error_score='raise',
    )
    message.result_mean_kfold(str(np.mean(n_scores['test_score'])))

    # Print the classification report
    acc = str(accuracy_score(y_test, y_pred))
    metrics = str(
        classification_report(y_test, y_pred, target_names=name_class)
    )
    message.result_accuracy(acc)
    message.result_metrics(metrics)
    save_metrics(acc=acc, metrics=metrics, dir_path=dir_path)
    del acc, metrics

    if should_save_confusion_matrix:
        conf_mtx = confusion_matrix(y_true=y_test, y_pred=y_pred)
        vmax = max(np.unique(y_test, return_counts=True)[1])
        save_confusion_matrix(
            conf_mtx=conf_mtx,
            name_class=name_class,
            vmax=vmax,
            dir_path=dir_path,
        )
    return


def classify(
    word: int,
    step: int,
    save_path: str,
    dir_path: str,
    dictonary: str = 'DNA',
    should_save_data: bool = True,
    should_save_model: bool = True,
    should_save_confusion_matrix: bool = True,
    chunk_size: int = 100,
):
    """
    Perform sequence classification pipeline.

    This is the main function for the classification module. It performs \
    sequence classification
    using the specified parameters and options. The function includes feature \
    extraction, model training, and optional saving of data, model, and \
    confusion matrix plot.

    Args:
        word (int): The length of each k-mer.
        step (int): The step size for moving the sliding window.
        save_path (str): The path to save the processed data and model files.
        dir_path (str): The path to the directory containing sequence data.
        dictonary (str): The DNA dictionary for k-mer analysis. Default is 'DNA'.
        should_save_data (bool, optional): Whether to save processed data. \
        Default is True.
        should_save_model (bool, optional): Whether to save the trained model. \
        Default is True.
        should_save_confusion_matrix (bool, optional): Whether to save the \
        confusion matrix plot. Default is True.
        chunk_size (int, optional): The chunk size for loading sequences. \
        Default is 100.

    Returns:
        Message class: A message confirming the classification pipeline has completed.
    """

    data_frame = extract_features(
        word=word,
        step=step,
        save_path=save_path,
        dir_path=dir_path,
        dictonary=dictonary,
        chunk_size=chunk_size,
    )

    # Process the feature matrix, ie, do MinMax scaler
    df_process, name_class = process_dataframe(
        data_frame=data_frame,
        dir_path=save_path,
        should_save_data=should_save_data,
        should_save_model=should_save_model,
    )
    sequence_classification(
        data_frame=df_process,
        name_class=name_class,
        dir_path=save_path,
        should_save_model=should_save_model,
        should_save_confusion_matrix=should_save_confusion_matrix,
    )
    return message.info_done()


def extract_features_to_predict(
    word: int,
    step: int,
    save_path: str,
    predict_seq_path: str,
    dictonary: str = 'DNA',
    chunk_size: int = 100,
) -> pd.DataFrame:
    """
    Extract features from sequences for prediction and return as a DataFrame.

    This function extracts features from sequences located in the specified \
    file for prediction purposes.
    The extracted features are based on the specified word length, step size, \
    and DNA dictionary.
    The extracted features are returned as a pandas DataFrame.

    Args:
        word (int): The length of each k-mer.
        step (int): The step size for moving the sliding window.
        save_path (str): The path to save the extracted features for prediction.
        predict_seq_path (str): The path to the file containing sequences \
        for prediction.
        dictonary (str): The DNA dictionary for k-mer analysis. Default is 'DNA'.
        chunk_size (int, optional): The chunk size for loading sequences. \
        Default is 100.

    Returns:
        pd.DataFrame: A pandas DataFrame containing the extracted features \
        for prediction.
    """

    message.info_start_prediction()
    progress = Progress(
        SpinnerColumn(),
        TaskProgressColumn(),
        TextColumn('[progress.description]{task.description}'),
        BarColumn(),
        TimeElapsedColumn(),
    )

    variants_kmers = load_variants_kmers(save_path=save_path)

    with progress:
        progress.add_task('[cyan]Loading sequences ...', total=None)

        features = load_sequences_classify(
            [predict_seq_path],
            word,
            step,
            dictonary,
            variants_kmers,
            True,
            chunk_size,
        )

    data_frame = pd.DataFrame(features)
    return data_frame


def process_dataframe_predict(
    data_frame: pd.DataFrame, load_ranges_path: str
) -> pd.DataFrame:
    """
    Process a DataFrame for prediction using MinMaxScaler.

    This function processes the provided DataFrame for prediction using the \
    MinMaxScaler object
    loaded from the specified path. It scales the data and restores column names \
    and 'ID' values.

    Args:
        data_frame (pd.DataFrame): The DataFrame to be processed for prediction.
        load_ranges_path (str): The path to load the MinMaxScaler object ranges.

    Returns:
        pd.DataFrame: A processed DataFrame with scaled values and restored 'ID' column.
    """

    minMax_scaler = load_ranges(load_ranges_path=load_ranges_path)

    id_values = data_frame['ID']

    data_frame.drop(columns=['ID'], axis=1, inplace=True)
    data_frame.replace([np.inf, -np.inf], 0, inplace=True)
    data_frame.replace(np.nan, 0, inplace=True)
    df_col_names = data_frame.columns

    df_minmax = minMax_scaler.transform(data_frame)
    data_frame = pd.DataFrame(df_minmax)
    del df_minmax

    data_frame.columns = df_col_names
    data_frame['ID'] = id_values

    return data_frame


def predict_data(
    data_frame: pd.DataFrame, load_model_path: str
) -> tuple[pd.Series, pd.Series]:
    """
    Predict classes using a trained RandomForestClassifier model.

    This function predicts classes for the provided data frame using the trained
    RandomForestClassifier model loaded from the specified path. It returns two
    pandas Series: one containing the predicted classes and the other containing
    the 'ID' values.

    Args:
        data_frame (pd.DataFrame): The data frame containing features for prediction.
        load_model_path (str): The path to load the trained \
        RandomForestClassifier model.

    Returns:
        tuple[pd.Series, pd.Series]: A tuple containing a Series of predicted \
        classes and a Series of 'ID' values.
    """

    rf_classifier = load_model(load_model_path=load_model_path)

    id_values = data_frame['ID']
    x_axis = data_frame.drop(columns='ID', axis=1)

    message.info_prediction()
    return rf_classifier.predict(x_axis), id_values


def predict(
    word: int,
    step: int,
    save_path: str,
    predict_seq_path: str,
    dir_path: str,
    dictonary: str,
    load_ranges_path: str,
    load_model_path: str,
    chunk_size: int = 100,
):
    """
    Predict sequence classes using a trained model.

    This function performs sequence class prediction using the specified \
    parameters and a trained model.
    It extracts features from sequences in the specified file for prediction \
    and scales them using the
    MinMaxScaler object loaded from the given path. The prediction results are \
    saved as a CSV file.

    Args:
        word (int): The length of each k-mer.
        step (int): The step size for moving the sliding window.
        save_path (str): The path to save the processed data and predictions.
        predict_seq_path (str): The path to the file containing sequences \
        for prediction.
        dir_path (str): The path to the directory containing additional files.
        dictonary (str): The DNA dictionary for k-mer analysis.
        load_ranges_path (str): The path to load the MinMaxScaler object ranges.
        load_model_path (str): The path to load the trained model.
        chunk_size (int, optional): The chunk size for loading sequences. \
        Default is 100.

    Returns:
        str: A message confirming the successful prediction and saving of results.
    """

    data_frame = extract_features_to_predict(
        word,
        step,
        save_path,
        predict_seq_path,
        dictonary,
        chunk_size,
    )

    data_frame = process_dataframe_predict(
        data_frame, load_ranges_path=load_ranges_path
    )

    predicted_data, id_values = predict_data(
        data_frame, load_model_path=load_model_path
    )

    return save_predict_data(id_values, predicted_data, dir_path)

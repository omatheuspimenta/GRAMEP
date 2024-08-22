"""
phylogenics module.

This module contains functions for generating phylogenetic tree visualizations from variant data.

Contents:
    * get_distance_linkage: Calculate distance matrix and linkage tree from a DataFrame.
    * get_phylogenic_tree: Generate and save a phylogenetic tree visualization.
    * get_phylogenic: Generate phylogenetic tree visualizations and save results.

Todo:
    * Implement tests.
"""
from pathlib import Path

import numpy as np
import numpy.typing as npt
import polars as pl
from gramep.data_io import (
    load_dataframe_phylo,
    save_confusion_matrix,
    save_newick,
)
from gramep.helpers import get_newick_str
from gramep.messages import Messages
from scipy.cluster.hierarchy import ClusterNode, linkage, to_tree
from sklearn.metrics import pairwise_distances
from toytree import save, tree

message = Messages()
"""
Set the Message class for logging.
"""


def get_distance_linkage(
    data_frame: pl.DataFrame,
) -> tuple[npt.NDArray[np.float64], ClusterNode, list[str]]:
    """
    Calculate distance matrix and linkage tree from a DataFrame.

    This function takes a DataFrame containing variant data, calculates the Hamming distance matrix,
    and constructs a linkage tree using hierarchical clustering. It returns the distance matrix,
    the linkage tree, and a list of variant names.

    Args:
        data_frame (pl.DataFrame): A DataFrame containing variant data.

    Returns:
        Tuple[np.ndarray, ClusterNode, list[str]]: A tuple containing:
            - The distance matrix as a NumPy array.
            - The linkage tree as a ClusterNode.
            - A list of variant names.

    Note:
        The 'ID' column in the DataFrame is assumed to represent variant names.
    """
    message.info_creating_distance_linkage()
    data_frame = data_frame.collect()
    variants_names = data_frame.select('ID').to_series().to_list()
    data_frame.drop_in_place(name='ID')

    df = data_frame.to_numpy()
    dist = pairwise_distances(df, metric='hamming', n_jobs=-2)

    link = linkage(df, method='complete', metric='hamming')
    link_tree = to_tree(link)

    return dist, link_tree, variants_names


# def get_phylogenic_tree(newick_str: str, save_path: str):
#     """
#     Generate and save a phylogenetic tree visualization.

#     This function takes a Newick tree string, converts it to an ultrametric tree,
#     and generates a phylogenetic tree visualization. The resulting visualization is saved as a PDF file.

#     Args:
#         newick_str (str): The Newick tree string representing the phylogenetic tree.
#         save_path (str): The directory path where the phylogenetic tree visualization will be saved.

#     Returns:
#         Message class: A message confirming the analysis was completed.

#     """
#     t = tree(newick_str)
#     canvas, _, _ = t.draw(
#         tip_labels_align=False,
#         tip_labels_style={'font-size': '3px', '-toyplot-anchor-shift': '5px'},
#         scale_bar=True,
#         shrink=10,
#         node_sizes=1,
#     )
#     path_dir = str(Path(save_path).parent) + '/phylogenics/'
#     path = Path(path_dir)
#     path.mkdir(mode=0o777, parents=True, exist_ok=True)

#     save(canvas, path_dir + '/tree.html')

#     return message.info_phylogenic_tree_saved(path_dir)


def get_phylogenic(save_path: str, save_heatmap: bool = False):
    """
    Generate phylogenetic tree visualizations and save results.

    This function orchestrates the process of generating phylogenetic tree visualizations from variant data.
    It retrieves a DataFrame from the specified 'save_path', calculates the distance matrix and linkage tree
    using hierarchical clustering, and saves the heatmap plot and Newick tree string.

    Args:
        save_path (str): The directory path where the phylogenetic visualizations will be saved.

    Returns:
        Message class: A message confirming the analysis was completed.

    Notes:
        - The DataFrame is obtained using the 'load_dataframe_phylo' function.
        - The distance matrix and linkage tree are obtained using the 'get_distance_linkage' function.
        - The heatmap is saved using the 'save_confusion_matrix' function.
        - The Newick tree string is saved using the 'save_newick' function.
        - The phylogenetic tree visualizations are generated using the 'get_phylogenic_tree' function.
    """
    message.info_start_objetive('phylogenics tree construction method ...')

    data_frame = load_dataframe_phylo(save_path=save_path)
    dist, linkage_tree, variants_names = get_distance_linkage(
        data_frame=data_frame
    )

    if save_heatmap:
        message.info_creating_heatmap()
        save_confusion_matrix(
            conf_mtx=dist,
            name_class=variants_names,
            vmax=np.max(dist),
            dir_path=save_path,
            is_conf_mtx=False,
        )

    newick_str = get_newick_str(linkage_tree, variants_names)
    save_newick(newick_str, save_path)
    # get_phylogenic_tree(newick_str, save_path)
    return message.info_done()

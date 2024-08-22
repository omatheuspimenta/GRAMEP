"""
entropy_utils
=============

This module provides utility functions for working with entropy and the maximum entropy
principle.

Contents:
    * entropy: Calculate the entropy of a probability distribution.
    * maxentropy: Calculate the maximum entropy threshold for a histogram.

Todo:
    * Implement tests.
    * Implement in Rust.
"""
from collections import defaultdict

import numpy as np
import numpy.typing as npt


def entropy(x: npt.NDArray[np.float64]) -> npt.NDArray[np.float64]:
    """
    Calculates the entropy of a array.

    Parameters:
        x: The probabilities p0 and p1.

    Returns:
        entropy_values: The calculated entropy.

    """
    # non_zero_indices:Tuple[np.int64] = np.nonzero(x)
    # return -(x[non_zero_indices] * np.log2(x[non_zero_indices]))
    return -(x[x != 0] * np.log2(x[x != 0]))


def maxentropy(
    kmers: defaultdict[str, int]
) -> tuple[np.float64, np.int64, np.int64, list[np.float64]]:
    """
    Kapur's threshold method.

    Kapur, J.N., P.K. Sahoo, and A.K.C. Wong. “A New Method for Gray-Level
    Picture Thresholding Using the Entropy of the Histogram.”
    Computer Vision, Graphics, and Image Processing 29,
    no. 3 (March 1985): 273–85. [DOI](https://doi.org/10.1016/0734-189X(85)90125-2).


    Parameters:
        kmers: frequency of occurrence for each kmers.

    Returns:
        maxEntropy: the maximum entropy value.
        threshold: the threshold value.
        frequency: the frequency at which the threshold occurs in the histogram.
        entropyCurve: all calculated entropies.

    """

    data: npt.NDArray[np.int64] = np.array(
        list(kmers.values()), dtype=np.int64
    )
    total_pixel: np.int64 = np.sum(data)
    # descending_data = -np.sort(-data.flatten())
    descending_data: npt.NDArray[np.int64] = -np.sort(-data)
    probs: npt.NDArray[np.float64] = descending_data / total_pixel

    p0: npt.NDArray[np.float64] = np.cumsum(probs)
    p1: npt.NDArray[np.float64] = np.cumsum(probs)[::-1]

    h0: npt.NDArray[np.float64] = entropy(p0)
    h1: npt.NDArray[np.float64] = entropy(p1)
    entropy_curve: list[np.float64] = [
        h0_i + h1_i for h0_i, h1_i in zip(h0, h1)
    ]

    max_entropy: np.float64 = np.max(entropy_curve)
    threshold: np.int64 = np.argmax(entropy_curve) + 1
    frequency: np.int64 = descending_data[threshold - 1]

    return max_entropy, threshold, frequency, entropy_curve

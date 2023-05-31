from ChromProcess import chromate
import numpy as np


def create_bins(time_array: np.ndarray, width: float) -> list[tuple[float, float]]:
    """
    Create bins using a mean-shift agglomeration algorithm.

    Parameters
    ----------
    time_array: np.ndarray

    Returns
    -------
    (rt_bin_borders, new_axis): np.ndarray, np.ndarray
    """

    sorted_times = np.sort(time_array)
    clusters = chromate.cluster(sorted_times, width)
    bin_borders = []
    for c in clusters:
        if len(c) == 1:
            lower = c[0] - width / 10
            upper = c[0] + width / 10
        else:
            lower = np.min(c)
            upper = np.max(c)

        bin_borders.append((lower, upper))

    bin_borders.sort(key=lambda x: x[0])

    return bin_borders

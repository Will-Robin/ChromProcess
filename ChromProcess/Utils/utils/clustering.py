import numpy as np


def cluster(values, bound=0.1):
    """
    Use an agglomerative algorithm to find clusters in a 1D array using bound
    to determine when a point is far enough from an average value to warrant
    creating a new cluster.

    Parameters
    ----------
    values: 1D numpy array
    bound: float

    Yields
    -------
    cluster: list
    """

    values = np.unique(np.sort(values))

    cluster = []
    for m in range(0, len(values)):
        if len(cluster) > 0:
            clust_av = np.average(cluster)

            if abs(values[m] - clust_av) > bound:
                yield cluster
                cluster = []

        cluster.append(values[m])

    yield cluster


def cluster_indices(values, bound=0.1):
    """
    Use an agglomerative algorithm to find clusters in a 1D array using bound
    to determine when a point is far enough from an average value to warrant
    creating a new cluster.

    This function returns the indices of the clusters in values.

    Parameters
    ----------
    values: 1D numpy array
    bound: float

    Yields
    -------
    cluster: list
    """

    sortedvalues = np.unique(np.sort(values))

    cluster = []
    for m in range(0, len(sortedvalues)):
        if len(cluster) == 0:
            pass
        else:
            clust_av = np.average(sortedvalues[cluster])

            if abs(sortedvalues[m] - clust_av) > bound:
                yield cluster
                cluster = []

        cluster.append(m)

    yield cluster


def create_bins(time_array: np.ndarray, width: float) -> tuple[np.ndarray]:
    """
    Create bins using a mean-shift agglomeration algorithm.

    Parameters
    ----------
    time_array: np.ndarray

    Returns
    -------
    (rt_bin_borders, new_axis): np.ndarray, np.ndarray
    """

    clusters = [c for c in cluster(time_array, bound=width)]
    bin_borders = []
    for c in clusters:
        mean = np.mean(c)
        if len(c) == 1:
            bin_borders.append((mean - width / 2, mean + width / 2))
        else:
            std = np.std(c, ddof=1)
            if std < width or np.isnan(std):
                mean = np.mean(c)
                bin_borders.append((mean - width / 2, mean + width / 2))
            else:
                bin_borders.append((min(c), max(c)))

    bin_borders.sort(key=lambda x: x[0])

    return bin_borders

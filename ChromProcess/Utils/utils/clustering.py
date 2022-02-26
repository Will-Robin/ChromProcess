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

    values = np.sort(values)

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

    sortedvalues = np.sort(values)

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

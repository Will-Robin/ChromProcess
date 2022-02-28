import numpy as np
from ChromProcess.Utils.utils import clustering as clust


def is_float(thing):
    """
    Test if an thing (e.g. str) can be converted to a float.

    Parameters
    ----------
    thing: any type

    Returns
    -------
    bool
    """
    try:
        float(thing)
        return True
    except ValueError:
        return False


def is_int(x):
    """
    Test if variable can be converted to an integer.

    Parameters
    ----------
    x: any type

    Returns
    -------
    bool
    """
    try:
        int(x)
        return True
    except ValueError:
        return False


def indices_from_boundary(data, start, end):
    """
    Get the indices of elements of the data array between start and end.

    Parameters
    ----------
    data: ndarray
        Array from which indice will be found
    start: float
        Lower boundary.
    end: float
        Upper boundary.

    Returns
    -------
    indices: ndarray
        Indices where start < data < end.
    """
    pre_indices = np.where((data >= start) & (data <= end))

    indices = pre_indices[0]

    return indices


def bin_dictionary(value_dict, stdev=0.001):
    """
    Combine dictionary values using a clustering algorithm.

    Parameters
    ----------
    value_dict: dict
        Keys are floats.

    Returns
    -------
    out_log: dict
        Binned dictionary
    """

    # create clusters of values for sorting
    sorted_values = sorted([*value_dict])

    clusters = []
    for c in clust.cluster(sorted_values, bound=stdev):
        clusters.append(c)

    out_log = {}
    for c in clusters:

        position = np.round(np.average(c), 2)

        out_log[position] = []

        for o in range(0, len(sorted_values)):
            value = sorted_values[o]
            if value in c:
                if len(out_log[position]) == 0:
                    out_log[position] = value_dict[value]
                else:
                    for count, p in enumerate(out_log[position]):
                        if p == 0:
                            out_log[position][count] = value_dict[value][count]
                        else:
                            pass

    return out_log


def peak_dict_to_spreadsheet(peak_dict, series_values, series_unit):
    """
    Convert a dictionary of peak series to a spreadsheet-like grid and a
    header.

    Parameters
    ----------
    peak_dict: dict
        Dictionary of peak series.
    series_values: list
        list of series values

    Returns
    -------
    peak_header, peak_grid_transposed: (list, list of lists)
    """

    peak_grid = []
    peak_header = [series_unit]

    peak_grid.append(series_values)

    # writing data
    peak_names = [*peak_dict]

    for s in peak_names:
        peak_header.append(s)
        peak_grid.append(peak_dict[s])

    peak_grid_transposed = [list(i) for i in zip(*peak_grid)]

    return peak_header, peak_grid_transposed


def peak_indices_to_times(time, picked_peaks):
    """
    Converts peak indices to times.

    Parameters
    ----------
    time: ndarray
        array of time, should match the indices.
    picked_peaks: dict
        dictionary containing list of indices of peak start, center, and end.

    Returns
    -----------
    peak_features: list
        list of lists containing times of peak start, center, and end.
    """

    peak_features = []
    for x in range(0, len(picked_peaks["Peak_indices"])):

        rt_ind = picked_peaks["Peak_indices"][x]
        start_ind = picked_peaks["Peak_start_indices"][x]
        end_ind = picked_peaks["Peak_end_indices"][x]

        retention_time = time[rt_ind]
        start = time[start_ind]
        end = time[end_ind]

        peak_params = [start, retention_time, end]

        peak_features.append(peak_params)

    return peak_features

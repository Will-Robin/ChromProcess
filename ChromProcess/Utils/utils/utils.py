import numpy as np
from ChromProcess.Utils.utils import clustering as clust

def is_float(thing):
    '''
    Test if an thing (e.g. str) can be converted to a float.
    '''
    try:
        float(thing)
        return True
    except ValueError:
        return False

def is_int(x):
    '''
    Test if variable can be converted to an integer.

    Parameters
    ----------
    x: any type

    Returns
    -------
    bool
    '''
    try:
        int(x)
        return True
    except ValueError:
            return False

def indices_from_boundary(time, start, end):
    '''
    Get the indices of elements of the time array between start and end.

    Parameters
    ----------
    time: ndarray
        Array from which indice will be found
    start: float
        Lower boundary.
    end: float
        Upper boundary.

    Returns
    -------
    indices: ndarray
        Indices where start < time < end.
    '''
    pre_indices = np.where(
                            (time >= p[0]) &
                            (time <= p[-1])
                            )
    indices = pre_indices[0]

    return indices

def bin_dictionary(value_dict, stdev = 0.001):
    '''
    Parameters 
    ----------
    value_dict: dict
        Keys are floats.

    Returns
    -------
    out_log: dict
        Binned dictionary
    '''

    # create clusters of values for sorting
    sorted_values = sorted([*value_dict])

    clusters = []
    for c in clust.cluster(sorted_values, bound = stdev):
        clusters.append(c)

    out_log = {}
    for c in clusters:

        position = np.round(np.average(c),2)

        out_log[position] = []

        for o in range(0,len(sorted_values)):
            value = sorted_values[o]
            if value in c:
                if len(out_log[position]) == 0:
                    out_log[position] = value_dict[value]
                else:
                    for count,p in enumerate(out_log[position]):
                        if p == 0:
                            out_log[position][count] = value_dict[values][count]
                        else:
                            pass

    return out_log


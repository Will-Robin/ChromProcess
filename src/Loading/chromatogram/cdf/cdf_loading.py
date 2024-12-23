"""
Functions for loading data from .cdf files (ANDI format).
"""

import numpy as np
from netCDF4 import Dataset
from pathlib import Path


def load_from_cdf(filename: str | Path, keys: list[str]) -> dict[str, np.ndarray]:
    """
    Extracts data from a .cdf file using the Dataset function from the netCDF4
    library.

    Parameters
    ----------
    filename: str
        File name of a .cdf file.
    keys: list[str]
        Key to a set of data in the .cdf file.

    Returns
    -------
    data_container: dict of numpy arrays
        Container for the extracted data.
    """

    data_container = dict()

    f = Dataset(filename, "r")
    f.set_auto_mask(False)
    for key in keys:
        data_container[key] = f.variables[key][:]
    f.close()

    return data_container

"""
This module contains 'low level' functions used by the rest of ChromProcess in
data analysis and processing operations.
"""

from .peak_finding import pick_peaks
from .signal_processing import deconvolution

from .utils import sum_error_prop
from .utils import mult_div_error_prop
from .utils import cluster, create_bins
from .utils import linear, inverse_linear
from .utils import quadratic, inverse_quadratic
from .utils import residual_squared_error
from .utils import inverse_quadratic_standard_error
from .utils import is_float
from .utils import is_int
from .utils import peak_dict_to_spreadsheet
from .utils import peak_indices_to_times
from .utils import bin_dictionary
from .utils import indices_from_boundary

__all__ = [
    "pick_peaks",
    "deconvolution",
    "sum_error_prop",
    "mult_div_error_prop",
    "cluster",
    "create_bins",
    "linear",
    "inverse_linear",
    "quadratic",
    "inverse_quadratic",
    "residual_squared_error",
    "inverse_quadratic_standard_error",
    "is_float",
    "is_int",
    "peak_dict_to_spreadsheet",
    "peak_indices_to_times",
    "bin_dictionary",
    "indices_from_boundary",
]

"""
General purpose functionality.
"""

from .error_propagation import sum_error_prop
from .error_propagation import mult_div_error_prop

from .clustering import cluster, create_bins

from .functions import linear, inverse_linear
from .functions import quadratic, inverse_quadratic
from .functions import residual_squared_error
from .functions import inverse_quadratic_standard_error

from .utils import is_float
from .utils import is_int
from .utils import peak_dict_to_spreadsheet
from .utils import peak_indices_to_times
from .utils import bin_dictionary
from .utils import indices_from_boundary

__all__ = [
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

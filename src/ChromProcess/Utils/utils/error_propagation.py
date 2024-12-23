import numpy as np


def sum_error_prop(errors):
    """
    Calculate error propagation for a sum operation.

    Parameters
    ----------
    errors: list of float

    Returns
    -------
    error: float
    """
    sq_err = [x**2 for x in errors]

    error = np.sqrt(sum(sq_err))

    return error


def mult_div_error_prop(averages, errors):
    """
    Calculate error propagation for a multiplication operation.

    Parameters
    ----------
    averages: list of float
    errors: list of float

    Returns
    -------
    error: float
    """

    sq_err = [(b / a) ** 2 for a, b in zip(averages, errors)]
    error = np.sqrt(sum(sq_err))

    return error

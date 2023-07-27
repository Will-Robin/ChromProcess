import numpy as np
from scipy.stats import norm
from scipy.optimize import curve_fit


def fit_pdf(time, signal, peaks, expected_heights, expected_widths, baseline):
    """
    Fitting sums of gaussian peaks to data using supplied peak indices.

    Parameters
    ----------
    time: np.ndarray
        time values
    signal: np.ndarray
        signal to be fit to
    peaks: list of peak positions
        list peak positions
    expected_heights: list | None
        Initial guess for peak magnitude.
    expected_widths: list | None
        Initial guess for peak widths (standard deviation).
    baseline: float
        Initial guess for the baseline.

    Returns
    -------
    popt:
        list of fitted values [[magnitude, positions, widths, baseline],]
    pcov:
        correlation matrix
    """

    # Create a time axis for each peak
    rep_time = np.tile(time, (len(peaks), 1)).T

    guess = np.hstack([expected_heights, peaks, expected_widths, [baseline]])

    bounds = (
        np.hstack(
            [
                np.full(len(peaks), signal.min()),
                np.full(len(peaks), time.min()),
                np.full(len(peaks), expected_widths.min()),
                0.0,
            ]
        ),
        np.hstack(
            [
                np.full(len(peaks), signal.max()),
                np.full(len(peaks), time.max()),
                np.full(len(peaks), expected_widths.max()),
                np.std(time, ddof=1),
            ]
        ),
    )

    popt, pcov = curve_fit(
        pdf_wrapper, rep_time, signal, p0=guess, bounds=bounds, method="trf"
    )

    return popt, pcov


def pdf_wrapper(time, *params):
    """
    Wrapper for `pdf` to unpack the parameter array for compatibility.

    time: np.ndarray
    *params:
        This will take the form of a tuple.
        Should be divisible into 4xn rows.
    """

    baseline = params[-1]

    unpack_params = np.reshape(params[:-1], (3, -1))

    return pdf(time, unpack_params[0], unpack_params[1], unpack_params[2], baseline)


def pdf(time, magnitude, positions, stdev, baseline):
    """
    Creates a summed, weighted multi-modal gaussian.

    Make sure time, stdevs, scale, magnitude, baseline are all the same
    dimension.

    Suggestion courtesy of Mathieu G. Baltussen.

    Parameters
    ----------
    time: np.ndarray
    magnitude: np.ndarray
    positions: np.ndarray
    stdev: np.ndarray
    baseline: float

    Returns
    -------
    np.ndarray
    """

    res = np.sum(
        magnitude * norm.pdf(time, loc=positions, scale=stdev) + baseline,
        axis=1,
    )

    return res

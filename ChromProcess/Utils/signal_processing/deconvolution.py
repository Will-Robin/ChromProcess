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

    # magnitude bounds
    mag_lower_bound = signal.min()
    mag_upper_bound = signal.max()

    # position bounds
    position_lower_bound = time.min()
    position_upper_bound = time.max()

    # width bounds
    width_lower_bound = expected_widths.min()
    width_upper_bound = expected_widths.max()

    # baseline bounds
    baseline_lower_bound = 0.0
    baseline_upper_bound = np.std(signal, ddof=1)

    bounds = (
        np.hstack(
            [
                np.full(len(peaks), mag_lower_bound),
                np.full(len(peaks), position_lower_bound),
                np.full(len(peaks), width_lower_bound),
                baseline_lower_bound,
            ]
        ),
        np.hstack(
            [
                np.full(len(peaks), mag_upper_bound),
                np.full(len(peaks), position_upper_bound),
                np.full(len(peaks), width_upper_bound),
                baseline_upper_bound,
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

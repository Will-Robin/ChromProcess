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
    mag_lower_range = np.full(len(peaks), min(signal.min(), min(peaks)))
    mag_upper_range = np.full(len(peaks), max(signal.max(), max(peaks)))

    # position bounds
    position_lower_range = np.full(len(peaks), time.min())
    position_upper_range = np.full(len(peaks), time.max())

    # width bounds
    width_lower_range = expected_widths * 0.5
    width_upper_range = expected_widths * 2.0

    # baseline bounds
    baseline_lower_bound = 0.0
    baseline_upper_bound = signal.max()

    bound_labels = ["magnitude", "position", "width", "baseline"]
    bounds = (
        np.hstack(
            [
                mag_lower_range,
                position_lower_range,
                width_lower_range,
                baseline_lower_bound,
            ]
        ),
        np.hstack(
            [
                mag_upper_range,
                position_upper_range,
                width_upper_range,
                baseline_upper_bound,
            ]
        ),
    )

    # Check bounds
    lower_bound_check = guess - bounds[0]
    if lower_bound_check.min() < 0.0:
        for x in range(0, len(lower_bound_check)):
            if lower_bound_check[x] < 0:
                print(
                    f"""
                An initial guess is lower than the lower bounds:
                      name: {bound_labels[x]}
                      value: {bounds[0][x]}
                      guess: {guess[x]}
                """
                )

    upper_bound_check = bounds[1] - guess
    if upper_bound_check.min() < 0.0:
        for x in range(0, len(upper_bound_check)):
            if upper_bound_check[x] < 0:
                print(
                    f"""
                An initial guess is higher than the upper bounds:
                      name: {bound_labels[x]}
                      value: {bounds[1][x]}
                      guess: {guess[x]}
                """
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

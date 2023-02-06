import numpy as np
from scipy.stats import norm
from scipy.optimize import curve_fit

from ChromProcess.Classes import Chromatogram


def deconvolute_region(chromatogram, region, num_peaks=1):
    """
    Deconvolute a region of a chromatogram.

    The fitting procedure assumes that the baselinne in the selected region of
    the chromatogram is constant.

    Parameters
    ----------
    chromatogram: Chromatogram
        Chromatogram
    region: list
        region of chromatogram under operation [lower bound, upper bound]
    num_peaks: int
        Number of peaks expected in the region.

    Returns
    -------
    result: ndarray
        Array of fitted values for each peak:
        [[magnitude, positions, widths, baseline],]
    """

    lower, upper = (region[0], region[1])

    # Cut out portion of signal
    idx = np.where((chromatogram.time > lower) & (chromatogram.time < upper))[0]

    time = chromatogram.time[idx]
    signal = chromatogram.signal[idx]

    # Find any peaks that have already been picked in the region.
    relevant_peaks = [
        chromatogram.peaks[p] for p in chromatogram.peaks if lower <= p <= upper
    ]
    peaks = [p.retention_time for p in relevant_peaks]

    # Create default values if no peaks have been found in the region.
    if len(peaks) > num_peaks:
        print("Deconvolution warning: deconvoluting more peaks than specified.")

    height_init = np.array([p.height for p in relevant_peaks])

    width_init = np.array(
        [np.std(chromatogram.time[p.indices], ddof=1) for p in relevant_peaks]
    )

    # Pad out any extra peaks
    if len(peaks) < num_peaks:
        pad_len = num_peaks - len(peaks)

        peaks = np.hstack((peaks, np.linspace(time.min(), time.max(), pad_len)))

        width_init = np.pad(
            width_init,
            (0, pad_len),
            mode="constant",
            constant_values=(0.0, time.std(ddof=1)),
        )

        height_init = np.pad(
            height_init,
            (0, pad_len),
            mode="constant",
            constant_values=(0.0, signal.mean()),
        )

    # Scale heights
    height_init /= len(height_init)

    baseline = signal.min()

    popt, pcov = fit_pdf(
        time,
        signal,
        peaks,
        height_init,
        width_init,
        baseline,
    )

    # Format output
    # -> [[magnitude, position, width, baseline],]
    result = np.vstack((np.reshape(popt[:-1], (3, -1)), np.full(len(peaks), popt[-1])))

    return result.T


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

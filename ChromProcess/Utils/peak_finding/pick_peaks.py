import numpy as np
from scipy.signal import find_peaks
from ChromProcess.Utils.signal_processing import signal_processing as sig
from chromate import find_peaks_wrapper


def pick_peaks(signal, distance=1, threshold=0.1, prominence=0.1, wlen=1001):
    """
    Peak detection routine.

    Parameters
    ----------
    y: np.ndarray
        1D amplitude data to search for peaks.

    Returns
    -------
    peak_indices: dict[str, np.ndarray]
        dict containing the indexes of the peaks that were detected and their
        start and end indices.
    """

    height = signal.max() * threshold

    smoothed = sig.adjacent_average(signal, wlen)

    peaks, starts, ends = find_peaks_wrapper(smoothed, height)

    return {
        "Peak_indices": peaks,
        "Peak_start_indices": starts,
        "Peak_end_indices": ends,
    }

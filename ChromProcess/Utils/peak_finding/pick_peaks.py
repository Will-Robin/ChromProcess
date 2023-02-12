import numpy as np
from scipy.signal import find_peaks
from ChromProcess.Utils.signal_processing import signal_processing as sig


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

    height = [signal.max() * threshold, signal.max()]

    smoothed = sig.adjacent_average(signal, wlen)

    peak_indices, _ = find_peaks(
        smoothed,
        threshold=None,
        distance=distance,
        height=height,
        prominence=prominence,
        wlen=wlen,
        width=None,
        rel_height=0.5,
        plateau_size=None,
    )

    diff = np.hstack(([0.0], np.diff(smoothed)))

    transition_is = diff > 0.0
    transition_idx = np.where(transition_is[1:] != transition_is[:-1])[0]

    idx = np.searchsorted(transition_idx, peak_indices, side="left")

    peak_starts = []
    peak_pos = []
    peak_ends = []
    for x in range(0, len(peak_indices)):
        i = idx[x]
        if i + 1 <= len(transition_idx) - 1:
            peak_ends.append(transition_idx[i + 1])
        else:
            peak_ends.append(-1)

        peak_starts.append(transition_idx[i - 1])
        peak_pos.append(peak_indices[x])

    return {
        "Peak_indices": peak_pos,
        "Peak_start_indices": peak_starts,
        "Peak_end_indices": peak_ends,
    }

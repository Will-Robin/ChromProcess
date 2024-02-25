import numpy as np
from typing import Optional
from scipy.signal import find_peaks
from ChromProcess.Utils.signal_processing import signal_processing as sig


def find_peak_bounds(
    signal: np.ndarray, peak_indices: list[int]
) -> tuple[list[int], list[int]]:
    """
    Find the start and end of a set of peaks in signal given a set of peak
    indices.

    Parameters
    ----------
    signal: np.ndarray
        Signal containing peaks.
    peak_indices: list[int]
        Indices of peaks in signal

    Returns
    -------
    peak_starts, peak_ends: tuple(list[int], list[int])
        Peak start and end indices.
    """

    # Find peak starts and ends
    diff = np.hstack(([0.0], np.diff(signal)))

    transition_is = diff > 0.0
    transition_idx = np.where(transition_is[1:] != transition_is[:-1])[0]

    idx = np.searchsorted(transition_idx, peak_indices, side="left")

    peak_starts = []
    peak_ends = []
    for x in range(0, len(peak_indices)):
        i = idx[x]
        if i + 1 <= len(transition_idx) - 1:
            peak_ends.append(transition_idx[i + 1])
        else:
            peak_ends.append(-1)

        peak_starts.append(transition_idx[i - 1])

    return peak_starts, peak_ends


def pick_peaks(
    signal: np.ndarray,
    smooth_width: int = 1,
    threshold: float = 0.0,
    distance: Optional[int] = 1,
    prominence: Optional[float] = None,
    wlen: Optional[int] = None,
    width: Optional[float] = None,
    rel_height: Optional[float] = None,
    plateau_size: Optional[int] = None,
) -> dict[str, list[int]]:
    """
    Peak detection routine.

    Parameters
    ----------
    signal: np.array
        Signal containing peaks
    smooth_width: int
        Width for smoothing window in number of data points.
    threshold: float
        Peaks below this fraction of the highest intensity of the chromatogram
        relative to the minimum of the signal will not be picked.
    distance: int | None
        See documentation for scipy.signal.find_peaks
    prominence: float
        See documentation for scipy.signal.find_peaks
    width: float
        See documentation for scipy.signal.find_peaks
    wlen: int
        See documentation for scipy.signal.find_peaks
    rel_height: float | None
        See documentation for scipy.signal.find_peaks
    plateau_size: float | None
        See documentation for scipy.signal.find_peaks

    Returns
    -------
    peak_indices: dict[str, list[int]]
        dict containing the indexes of the peaks that were detected and their
        start and end indices.
    """

    height = signal.min() + (signal.max() - signal.min()) * threshold

    if smooth_width > 1:
        smoothed = sig.adjacent_average(signal, smooth_width)
    else:
        smoothed = signal

    peak_indices, _ = find_peaks(
        smoothed,
        threshold=None,
        distance=distance,
        height=height,
        prominence=prominence,
        wlen=wlen,
        width=width,
        rel_height=rel_height,
        plateau_size=plateau_size,
    )

    peak_starts, peak_ends = find_peak_bounds(smoothed, peak_indices)

    return {
        "Peak_indices": peak_indices,
        "Peak_start_indices": peak_starts,
        "Peak_end_indices": peak_ends,
    }

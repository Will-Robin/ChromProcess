import numpy as np

from ChromProcess import Classes
from ChromProcess.Utils.peak_finding.pick_peaks import find_peak_boundaries_look_ahead
from ChromProcess.Loading.peak_collection.peak_collection_from_csv import (
    peak_collection_from_csv,
)


def peak_indices_from_file(peak_file, chromatogram):
    """
    Get peak retention time indices from a PeakCollections file.

    Parameters
    ----------
    chromatogram: Classes.Chromatogram object
        Chromatogram object to which the peaks should be added.
    peak_file: string
        Location of peak collection csv

    Returns
    -------
    peak_indices: ndarray
        array of peak indices
    """

    peak_collection = peak_collection_from_csv(peak_file, round_digits=7)
    peak_retention_times = [p.retention_time for p in peak_collection.peaks]

    # Find the indices of the retention times in the
    # chromatogram time axis.
    peaks_indices = np.empty(0, dtype="int64")
    for p_rt in peak_retention_times:
        # Find the closest retention time value in to p_rt
        delta_rt = np.abs(chromatogram.time - p_rt)
        idx = delta_rt.argmin()
        peaks_indices = np.append(peaks_indices, idx)

    return peaks_indices


def peak_from_csv(peak_file, chromatogram, look_ahead=12):
    """
    Load sets of peak indices from a set of retention times stored in a
    PeakCollections file rather than directly from a chromatogram.

    Parameters
    ----------
    chromatogram: Classes.Chromatogram object
        Chromatogram object to which the peaks should be added
    peak_file: string
        Location of peak collection csv
    peak_window: int
        Number of spaces the function will search through to
        find the start or end of a peak

    Returns
    -------
    peaks: list of Peak objects
    """

    peaks_indices = peak_indices_from_file(peak_file, chromatogram)

    peak_starts, peak_ends = find_peak_boundaries_look_ahead(
        chromatogram.signal, peaks_indices, look_ahead=look_ahead
    )

    picked_peaks = {
        "Peak_indices": peaks_indices,
        "Peak_start_indices": peak_starts,
        "Peak_end_indices": peak_ends,
    }

    peaks = []
    for x in range(0, len(picked_peaks["Peak_indices"])):

        pk_idx = picked_peaks["Peak_indices"][x]
        start_idx = picked_peaks["Peak_start_indices"][x]
        end_idx = picked_peaks["Peak_end_indices"][x]

        retention_time = chromatogram.time[pk_idx]
        start = chromatogram.time[start_idx]
        end = chromatogram.time[end_idx]

        indices = np.where((chromatogram.time >= start) & (chromatogram.time <= end))[0]

        peaks.append(Classes.Peak(retention_time, start, end, indices=indices))

    return peaks

import numpy as np

from ChromProcess.Utils.utils.utils import peak_indices_to_times
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
        delta_rt = np.abs(chromatogram.time - float(p_rt))
        nearest_value_idx = delta_rt.argmin()
        #search the 5 nearest value for the highest peak, in case the time input was not exact
        idx = nearest_value_idx - 2 + chromatogram.signal[nearest_value_idx-2:nearest_value_idx+2].argmax() 
        peaks_indices = np.append(peaks_indices, idx)
    return peaks_indices

def peak_boundaries_from_file(chromatogram, peak_file):
    peak_collection = peak_collection_from_csv(peak_file,round_digits=7)
    peak_start_retention_times = [p.start for p in peak_collection.peaks]

    peak_starts = np.empty(0,dtype='int64') #the functions written to find the start and end of the peak relies on indexes, so we need to convert the retention times to their corresponding indexes.
    for ps_rt in peak_start_retention_times:
        delta_rt = np.abs(chromatogram.time-float(ps_rt))
        idx = delta_rt.argmin()
        peak_starts = np.append(peak_starts, idx)
    
    peak_end_retention_times = [p.end for p in peak_collection.peaks]

    peak_ends = np.empty(0,dtype='int64') #the functions written to find the start and end of the peak relies on indexes, so we need to convert the retention times to their corresponding indexes.
    for pe_rt in peak_end_retention_times:
        delta_rt = np.abs(chromatogram.time-float(pe_rt))
        idx = delta_rt.argmin()
        peak_ends = np.append(peak_ends, idx)
    
    return peak_starts, peak_ends




def peak_from_csv(peak_file, chromatogram, look_ahead = 12):
    '''
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
    peak_features: list
        list of list containing times of peak start, center, and end.
    '''

    peaks_indices = peak_indices_from_file(peak_file, chromatogram)

    peak_starts, peak_ends = find_peak_boundaries_look_ahead(
        chromatogram.signal, peaks_indices, look_ahead=look_ahead
    )

    picked_peaks = {
        "Peak_indices": peaks_indices,
        "Peak_start_indices": peak_starts,
        "Peak_end_indices": peak_ends,
    }

    peak_features = peak_indices_to_times(chromatogram.time, picked_peaks)

    return peak_features

import numpy as np
from ChromProcess.Utils.utils.utils import peak_indices_to_times
from ChromProcess.Utils.peak_finding.pick_peaks import find_peak_boundaries_look_ahead
from ChromProcess.Loading.peak_collection.peak_collection_from_csv import peak_collection_from_csv

def peak_rt_from_file(peak_file, chromatogram):
    '''
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
    '''

    peak_collection = peak_collection_from_csv(peak_file,round_digits=7)
    peak_retention_times = [p.retention_time for p in peak_collection.peaks]

    # Find the indices of the retention times in the 
    # chromatogram time axis.
    peaks_indices = np.empty(0,dtype='int64') 
    for p_rt in peak_retention_times:
        # This covers any exact matches
        idx = np.where(chromatogram.time == p_rt)[0]

        if len(idx) == 0: 
            # If the exact value of p_rt is not in the time axis, 
            # find the closest time value.
            delta_rt = np.abs(chromatogram.signal - p_rt)
            # Create in list; similar behaviour to return of np.where()[0]
            idx = [delta_rt.argmin()]

        peaks_indices = np.append(peaks_indices, idx[0])

    return peaks_indices

def peak_from_csv(peak_file, chromatogram, look_ahead = 12):
    '''
    Get the peak retention times from a PeakCollections file rather 
    than directly from a chromatogram.

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

    peaks_indices = peak_rt_from_file(
                                    peak_file,
                                    chromatogram
                                    )

    peak_starts, peak_ends = find_peak_boundaries_look_ahead(
                                                            chromatogram.signal,
                                                            peaks_indices,
                                                            look_ahead = 1
                                                            )
    picked_peaks = {
                    'Peak_indices': peaks_indices,
                    'Peak_start_indices': peak_starts,
                    'Peak_end_indices': peak_ends
                    }

    peak_features = peak_indices_to_times(
                                        chromatogram.time,
                                        picked_peaks
                                        )

    return peak_features

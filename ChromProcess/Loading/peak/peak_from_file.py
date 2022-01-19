import numpy as np


def peak_rt_from_file(chromatogram, Peakfile):
    '''
    Get peak retention time indices from a PeakCollections file.
    chromatogram: Chromatogram object
        Chromatogram object to which the peaks should be added
    Peakfile: string
        Location of peak collection csv
    
    Returns
    -------
    peak_indices: ndarray  
        array of peak indices
    '''
    import pandas as pd
    df = pd.read_csv(Peakfile,index_col=False)
    peak_retention_times = df.iloc[3:,0].values #This function is very specific to the file type
    peaks_indices = np.empty(0,dtype='int64') #the functions written to find the start and end of the peak relies on indexes, so we need to convert the retention times to their corresponding indexes.
    for p_rt in peak_retention_times:
        idx = np.where(chromatogram.time == float(p_rt)) #This covers any exact matches
        if len(idx[0]) == 0: #no exact match found, this should mean that a value was manually inserted, the function will now search the closest time value.
            search_idx = np.searchsorted(chromatogram.time, float(p_rt), "left")
            max_signal = np.max(chromatogram.signal[search_idx-4:search_idx+4]) #search the 9 nearest value for the highest peak, this gives a bit of play for time input
            idx = np.where(chromatogram.signal == max_signal)
        peaks_indices = np.append(peaks_indices, int(idx[0]))
    return peaks_indices


def peak_from_file(chromatogram, Peakfile, peak_window = 12):
    '''
    Get the peak retention times from a PeakCollections file rather than directly from a chromatogram.
    Peak boundaries will be generated. 

    chromatogram: Chromatogram object
        Chromatogram object to which the peaks should be added
    Peakfile: string
        Location of peak collection csv
    peak_window: int
        Number of spaces the function will search through to find the start or end of a peak

    Returns
    -------
    peak_features: list
        list of list containing times of peak start, center, and end.
    '''
    import ChromProcess.Utils.peak_finding.pick_peaks as pf
    import ChromProcess.Utils.utils.utils as utils

    peaks_indices =  peak_rt_from_file(chromatogram, Peakfile)
    peak_starts, peak_ends = pf.find_peak_boundaries(chromatogram.signal, peaks_indices, peak_window=1)
    picked_peaks = {'Peak_indices':peaks_indices, 'Peak_start_indices':peak_starts, 'Peak_end_indices':peak_ends}
    peak_features = utils.peak_indices_to_times(chromatogram.time, picked_peaks)
    
    return peak_features
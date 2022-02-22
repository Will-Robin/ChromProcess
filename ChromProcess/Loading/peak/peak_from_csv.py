import numpy as np

from ChromProcess.Loading.peak_collection.peak_collection_from_csv import peak_collection_from_csv


def peak_rt_from_file(chromatogram, peak_file):
    '''
    Get peak retention time indices from a PeakCollections file.
    chromatogram: Chromatogram object
        Chromatogram object to which the peaks should be added
    peak_file: string
        Location of peak collection csv
    
    Returns
    -------
    peak_indices: ndarray  
        array of peak indices
    '''
    peak_collection = peak_collection_from_csv(peak_file,round_digits=7)
    peak_retention_times = [p.retention_time for p in peak_collection.peaks]

    peaks_indices = np.empty(0,dtype='int64') #the functions written to find the start and end of the peak relies on indexes, so we need to convert the retention times to their corresponding indexes.
    for p_rt in peak_retention_times:
        idx = np.where(chromatogram.time == float(p_rt)) #This covers any exact matches
        if len(idx[0]) == 0: #no exact match found, this should mean that a value was manually inserted, the function will now search the closest time value.
            search_idx = np.searchsorted(chromatogram.time, float(p_rt), "left")
            max_signal = np.max(chromatogram.signal[search_idx-4:search_idx+4]) #search the 9 nearest value for the highest peak, this gives a bit of play for time input
            idx = np.where(chromatogram.signal == max_signal)
        peaks_indices = np.append(peaks_indices, int(idx[0]))
    return peaks_indices

def peak_boundaries_from_file(chromatogram, peak_file):
    peak_collection = peak_collection_from_csv(peak_file,round_digits=7)
    peak_start_retention_times = [p.start for p in peak_collection.peaks]

    peak_starts = np.empty(0,dtype='int64') #the functions written to find the start and end of the peak relies on indexes, so we need to convert the retention times to their corresponding indexes.
    for ps_rt in peak_start_retention_times:
        idx = np.where(chromatogram.time == float(ps_rt))[0] #This covers any exact matches
        if len(idx) == 0: #no exact match found, this should mean that a value was manually inserted, the function will now search the closest time value.
            idx = np.searchsorted(np.round(chromatogram.time,7), float(ps_rt), "left") #The rounding should not be necessary, but without it searchsorted can give wrong values
        peak_starts = np.append(peak_starts, int(idx))
    
    peak_end_retention_times = [p.end for p in peak_collection.peaks]

    peak_ends = np.empty(0,dtype='int64') #the functions written to find the start and end of the peak relies on indexes, so we need to convert the retention times to their corresponding indexes.
    for pe_rt in peak_end_retention_times:
        idx = np.where(chromatogram.time == float(pe_rt))[0] #This covers any exact matches
        if len(idx) == 0: #no exact match found, this should mean that a value was manually inserted, the function will now search the closest time value.
            idx =  np.searchsorted(np.round(chromatogram.time,7), float(pe_rt), "left")
        peak_ends = np.append(peak_ends, int(idx))
    
    return peak_starts, peak_ends




def peak_from_csv(chromatogram, peak_file, look_ahead = 12):
    '''
    Get the peak retention times from a PeakCollections file rather than directly from a chromatogram.
    Peak boundaries will be generated. 

    chromatogram: Chromatogram object
        Chromatogram object to which the peaks should be added
    peak_file: string
        Location of peak collection csv
    peak_window: int
        Number of spaces the function will search through to find the start or end of a peak

    Returns
    -------
    peak_features: list
        list of list containing times of peak start, center, and end.
    '''
    from ChromProcess.Utils.peak_finding.pick_peaks import find_peak_boundaries_look_ahead
    from ChromProcess.Utils.utils.utils import peak_indices_to_times
    
    peaks_indices =  peak_rt_from_file(chromatogram, peak_file)
    peak_starts, peak_ends = find_peak_boundaries_look_ahead(chromatogram.signal, peaks_indices, look_ahead=1)
    picked_peaks = {'Peak_indices':peaks_indices, 'Peak_start_indices':peak_starts, 'Peak_end_indices':peak_ends}
    peak_features = peak_indices_to_times(chromatogram.time, picked_peaks)
    
    return peak_features
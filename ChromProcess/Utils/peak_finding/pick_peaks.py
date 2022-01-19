import numpy as np
from ChromProcess.Utils.signal_processing import signal_processing as sig


def find_peak_boundaries(signal, peaks_indices, peak_window=1):
    '''
    Find the start and end of a peak using the peak rt and signal.
    The function does a stepwise search checking if the next value is lower than the current, ending if it is not.
    To search noisy data which might have a some small peaks on the slope of a peak the peak_window can be used, extending the search range.
    Note: The function is cannot search beyond the signal boundaries, when using region_peak_picks this could cause early peak boundaries to be cut off early.
    signal: 1D array
        The signal in which to search for the peak start and end. The smoothed signal can be used.
    peaks_indices: ndarray
        the indices of the centre of the peaks.
    peak_window: int, optional
        The number of indices the function looks ahead to see if the signal value decreases

    Returns
    --------

    peak_starts: list
        Indices of peak start.
    peak_ends: list
        Indices of peak ends.
    '''

    peak_starts = []
    for n in range(0,len(peaks_indices)):
        cursor = peaks_indices[n]-2
        while not cursor-1 < peak_window and signal[cursor] > np.min(signal[cursor-peak_window:cursor]):
            cursor -= 1
        while signal[cursor] > signal[cursor-1] and not cursor-1 != 0: #if the search is stopped because the peak window is nearing the border continue stepwise search
            cursor -= 1
        peak_starts.append(cursor)

    peak_ends = []

    for n in range(0,len(peaks_indices)):
        cursor = peaks_indices[n]+2
        while not cursor+1+peak_window > len(signal) and signal[cursor] > np.min(signal[cursor:cursor+peak_window]):
            cursor += 1  
        while signal[cursor] > 0 and not cursor + 1 < len(signal):
            cursor += 1
        peak_ends.append(cursor)

    return peak_starts, peak_ends

def find_peaks_scipy(signal, threshold=0.1, min_dist=1, max_inten = 1e100, prominence = 0.7, wlen = 1001, peak_window = 12):
    '''Peak finding function that relies on scipy.signal.find_peaks rather than the first derivative.
    This allows to specify peak prominence, making it more suited for a noisy data set.

    Parameters
    ----------
    time, signal: 1D arrays
        Signal from which peaks will be picked and its corresponding time axis.
    thres : float between [0., 1.]
        Normalized threshold. Only the peaks with amplitude higher than the
        threshold will be detected.
    min_dist : int
        Minimum distance between each detected peak. The peak with the highest
        amplitude is preferred to satisfy this constraint.
    max_inten: float
        peaks will not be detected above this threshold.
    min_inten: float
        peaks will not be detected below this threshold.
    Prominence: number or ndarray or sequence, optional
        Required prominence of peaks. Either a number, None, an array matching 
        x or a 2-element sequence of the former. The first element is always 
        interpreted as the minimal and the second, if supplied, as the maximal required prominence.
    wlen: int, optional
        The number of indices to search for when determining the peak promince. Recommend to be roughly equal to average peak width.
    peak_window: int
        Number of spaces the function will search through to find the start or end of a peak

    Returns
    -------
    ndarray
        Array containing the indexes of the peaks that were detected
    '''
    
    import numpy as np
    from scipy.signal import find_peaks

    smooth_signal = sig.savitzky_golay(signal, 7, 3, deriv=0, rate=1) #smooth out the function before starting the picking
    peaks_indices, properties = find_peaks(smooth_signal,  
                                    distance = min_dist, 
                                    height = [max(smooth_signal)*threshold,max_inten], 
                                    prominence = prominence, 
                                    wlen = wlen)
    peak_starts, peak_ends = find_peak_boundaries(smooth_signal, peaks_indices, peak_window=peak_window)

    return {'Peak_indices':peaks_indices, 'Peak_start_indices':peak_starts, 'Peak_end_indices':peak_ends}

def find_peaks(signal, thres=0.1, min_dist=1, min_inten = -1e100):

    '''
    Peak detection routine.

    Modified from PeakUtils:
    https://github.com/atjacobs/PeakUtils/tree/master/peakutils
    Finds the peaks in *y* by taking its first order difference. By using
    *thres* and *min_dist* parameters, it is possible to reduce the number of
    detected peaks.

    Parameters
    ----------
    y: ndarray
        1D amplitude data to search for peaks.
    thres: float between [0., 1.]
        Normalized threshold. Only the peaks with amplitude higher than the
        threshold will be detected.
    min_dist: int
        Minimum distance between each detected peak. The peak with the highest
        amplitude is preferred to satisfy this constraint.
    min_inten: float
        peaks will not be detected below this threshold.

    Returns
    -------
    dict
        dict containing the indexes of the peaks that were detected and their
        start and end indices.
    '''

    thres *= np.max(signal) - np.min(signal)

    # find the peaks by using the first order difference
    diff = np.diff(signal)

    smoothed_diff = sig.savitzky_golay(diff, 5, 3, deriv=0, rate=1)

    pre_peak_inds = np.where(
                              (np.hstack([smoothed_diff, 0.]) < 0.)
                            & (np.hstack([0., smoothed_diff]) > 0.)
                            & (signal > thres)
                            & (signal > min_inten))

    peaks_indices = pre_peak_inds[0]

    if peaks_indices.size > 1 and min_dist > 1:
        highest = peaks_indices[np.argsort(signal[peaks_indices])][::-1]
        rem = np.ones(signal.size, dtype=bool)
        rem[peaks_indices] = False

        for peak in highest:
            if not rem[peak]:
                sl = slice(max(0, peak - min_dist), peak + min_dist + 1)
                rem[sl] = True
                rem[peaks_indices] = False

        peaks_indices = np.arange(signal.size)[~rem]

    # find beginning of peaks
    peak_starts = []
    for n in range(0,len(peaks_indices)):
        cursor = peaks_indices[n]-2
        while diff[cursor] > 0:
            cursor -= 1
        peak_starts.append(cursor)

    # find end of peaks
    peak_ends = []
    for n in range(0,len(peaks_indices)):
        cursor = peaks_indices[n]+2
        if int(cursor) >= len(diff):
            cursor = len(diff)-1

        while diff[cursor] < 0:
            cursor += 1
            if int(cursor) >= len(diff):
                cursor = len(diff)-1
                break
        peak_ends.append(cursor)

    return {'Peak_indices':peaks_indices, 'Peak_start_indices':peak_starts, 'Peak_end_indices':peak_ends}

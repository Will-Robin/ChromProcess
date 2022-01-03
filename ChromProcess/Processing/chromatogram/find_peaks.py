from ChromProcess.Utils.utils import utils
from ChromProcess.Utils.peak_finding import pick_peaks as pfind

def find_peaks_in_region(chromatogram, start, end, threshold = 0.1):
    '''
    Find peaks within the chromatogram between start and end retention times.

    Parameters
    ----------
    chromatogram: ChromProcess Chromatogram object

    start: float
        Start of time region (< end).
    end: float
        End of time region (> start).
    threshold: float
        Peaks below this fraction of the highest intensity of the chromatogram
        will not be picked.

    Returns
    -------
    None
    '''
    if start > end:
        print(f'peak start ({start}) > peak end, ({end}) returning None')
        return None

    inds = utils.indices_from_boundary(chromatogram.time, start, end)

    time = chromatogram.time[inds]
    signal = chromatogram.signal[inds]

    picked_peaks = pfind.find_peaks(
                                    signal,
                                    thres = threshold
                                    )

    peak_features = []
    for x in range(0,len(picked_peaks['Peak_indices'])):

        rt_ind = picked_peaks['Peak_indices']
        start_ind = picked_peaks['Peak_start_indices']
        end_ind = picked_peaks['Peak_end_indices']

        retention_time = time[rt_ind][0]
        start = time[start_ind][0]
        end = time[end_ind][0]

        peak_params = [start, retention_time, end]

        peak_features.append(peak_params)

    return peak_features

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
    print("you are here")
    picked_peaks = pfind.find_peaks(
                                    signal,
                                    thres = threshold
                                    )
    peak_features = utils.peak_indices_to_times(time,picked_peaks)

    return peak_features


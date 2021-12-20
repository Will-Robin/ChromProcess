from ChromProcess.Utils.utils import utils
from ChromProcess.Utils.utils import peak_finding as pfind

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

    rts, peak_times = pfind.find_peaks(
                                        time,
                                        signal,
                                        threshold = threshold
                                        )
    peak_features = []
    for rt, times in zip(rts,peak_times):
        feature = [times[0], rt, times[-1]]
        peak_features.append(feature)

    return peak_features

from ChromProcess import Classes
from ChromProcess.Utils.utils import utils
from ChromProcess.Utils.peak_finding import pick_peaks as pfind


def find_peaks_in_region(chromatogram, start, end, threshold=0.1):
    """
    Find peaks within the chromatogram between start and end retention times.

    Parameters
    ----------
    chromatogram: ChromProcess Chromatogram object

    start: float
        Start of time region (< end).
    end: float End of time region (> start).
    threshold: float
        Peaks below this fraction of the highest intensity of the chromatogram
        will not be picked.

    Returns
    -------
    peaks: list of Peak objects
    """

    if start > end:
        print(f"peak start ({start}) > peak end, ({end}) returning None")
        return None

    inds = utils.indices_from_boundary(chromatogram.time, start, end)

    time = chromatogram.time[inds]
    signal = chromatogram.signal[inds]

    picked_peaks = pfind.find_peaks(signal, thres=threshold)

    peaks = []
    for x in range(0, len(picked_peaks["Peak_indices"])):

        pk_idx = picked_peaks["Peak_indices"][x]
        start_idx = picked_peaks["Peak_start_indices"][x]
        end_idx = picked_peaks["Peak_end_indices"][x]

        retention_time = time[pk_idx]
        start = time[start_idx]
        end = time[end_idx]

        peaks.append(Classes.Peak(retention_time, start, end, indices=[]))

    return peaks

import numpy as np
from ChromProcess.Classes import Peak
from ChromProcess.Classes import Chromatogram


def peak_from_chromatogram(
    chrom: Chromatogram, start: float, end: float
) -> Peak | None:
    """
    Create a peak using the boundaries defined within a chromatogram.

    Parameters
    ----------
    chrom: Chromatogram
    start: float
        Start of the peak
    end: float
        End of the peak

    Returns
    -------
    peak: Peak
        peak created from the bounds
    """

    if start > end:
        print(f"peak start ({start}) > peak end, ({end}) returning None")
        return None

    inds = np.where((chrom.time > start) & (chrom.time < end))[0].tolist()

    timeseg = chrom.time[inds]
    sigseg = chrom.signal[inds]

    peak_idx = np.argmax(sigseg)
    retention_time = timeseg[peak_idx]

    start = timeseg[0]
    end = timeseg[-1]
    peak = Peak(retention_time, start, end, indices=inds)
    peak.set_height(chrom.signal[peak_idx])

    return peak

import numpy as np
from ChromProcess import Classes

def peak_from_chromatogram(chrom, start, end):
    '''
    Create a peak using the boundaries defined within chrom

    chrom: Chromatogram object
    start: float
        Start of the peak
    end: float
        End of the peak
    
    peak: Peak object
        peak created from the bounds
    '''

    if start > end:
        print(f'peak start ({start}) > peak end, ({end}) returning None')
        return None

    inds = np.where(
            (chrom.time > start) & 
            (chrom.time < end)
            )[0]

    timeseg = chrom.time[inds]
    sigseg  = chrom.signal[inds]

    peak_idx = np.argmax(sigseg)
    retention_time = timeseg[peak_idx] 
 
    peak = Classes.Peak(retention_time, inds)

    return peak


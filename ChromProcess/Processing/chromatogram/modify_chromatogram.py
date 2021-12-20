from ChromProcess.Utils.utils import utils
from ChromProcess.Processing.chromatogram import find_peaks

def add_peaks_to_chromatogram(peak_times, chromatogram):
    '''
    Parameters
    ----------
    peak_times: list
        [start, peak, end] in units of the chromatogram's
        retention time axis.

    chromatogram: Chromatogram object

    Returns
    ------
    None
    '''

    time = chromatogram.time
    for p in peak_times:
        start, retention_time, end = p[0], p[1], p[2]

        idx = utils.indices_from_boundary(time, start, end)

        peak = Classes.Peak(retention_time, idx)
        chromatogram.peaks[retention_time] = peak

def integrate_chromatogram_peaks(chromatogram, baseline_subtract = False):
    '''
    Parameters
    ----------

    chromatogram: Chromatogram object

    Returns
    ------
    None
    '''

    for p in chromatogram.peaks:
        chromatogram.peaks[p].get_integral(
                                        chromatogram, 
                                        baseline_subtract = baseline_subtract
                                        )

def internal_ref_integral(chromatogram, internal_ref_region):
    '''
    Finds and adds internal reference information into a chromatogram.

    Parameters
    ----------
    series: Chromatogram_Series object
        Object containing chromatograms and associated series data which is
        modified by the function.

    Returns
    ------
    None
    '''

    start = internal_ref_region[0]
    end = internal_ref_region[1]
    if len(chromatogram.peaks) == 0:
        find_peaks.find_peaks_in_region(
                                        chromatogram, 
                                        start, 
                                        end,
                                        threshold = 0.1
                                        )
    else:
        pass

    for p in chromatogram.peaks:
        if p > start and p < end:
            rt = chromatogram.peaks[p].retention_time
            indices = chromatogram.peaks[p].indices
            chromatogram.internal_reference = Classes.Peak(rt, indices)
            chromatogram.internal_reference.get_integral(chromatogram)


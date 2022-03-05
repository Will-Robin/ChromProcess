import numpy as np
from ChromProcess.Processing.chromatogram import find_peaks


def add_peaks_to_chromatogram(peaks, chromatogram):
    """
    Add peaks to a chromatogram (modifies the chromatogram in place).

    Parameters
    ----------
    peaks: list of Peak objects

    chromatogram: Chromatogram object

    Returns
    ------ None
    """

    for peak in peaks:
        rt = peak.retention_time
        indices = np.where(
            (chromatogram.time >= peak.start) & (chromatogram.time <= peak.end)
        )[0]
        peak.indices = indices
        chromatogram.peaks[rt] = peak


def integrate_chromatogram_peaks(chromatogram, baseline_subtract=False):
    """
    Integrate all of the peaks in a chromatogram (modifies them in place).

    Parameters
    ----------
    chromatogram: Classes.Chromatogram object
        Chromatogram containing peaks.
    baseline_subtract: bool
        Whether to perform a local baseline subtraction on the peak.

    Returns
    ------
    None
    """

    for p in chromatogram.peaks:
        chromatogram.peaks[p].get_integral(
            chromatogram, baseline_subtract=baseline_subtract
        )


def internal_standard_integral(chromatogram, is_start, is_end):
    """
    Finds and adds internal standard information into a chromatogram.

    Parameters
    ----------
    series: Chromatogram_Series object
        Object containing chromatograms and associated series data which is
        modified by the function.
    is_start: float
        Start of the region of the chromatogram in which the internal standard
        is found.
    is_end: float
        End of the region of the chromatogram in which the internal standard
        is found.

    Returns
    ------
    None
    """

    peaks = find_peaks.find_peaks_in_region(
        chromatogram, is_start, is_end, threshold=0.1
    )

    peak = peaks[0]

    peak.indices = np.where(
        (chromatogram.time >= peak.start) & (chromatogram.time <= peak.end)
    )[0]

    peak.get_integral(chromatogram)

    chromatogram.internal_standard = peak

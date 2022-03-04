from ChromProcess import Classes
from ChromProcess.Utils.utils import utils
from ChromProcess.Processing.chromatogram import find_peaks


def add_peaks_to_chromatogram(peak_times, chromatogram):
    """
    Add peaks to a chromatogram (modifies the chromatogram in place).

    Parameters
    ----------
    peak_times: list
        [start, peak, end] in units of the chromatogram's
        retention time axis.

    chromatogram: Chromatogram object

    Returns
    ------
    None
    """

    time = chromatogram.time
    for p in peak_times:
        start, retention_time, end = p[0], p[1], p[2]

        idx = utils.indices_from_boundary(time, start, end)

        peak = Classes.Peak(retention_time, start, end, indices=idx)
        chromatogram.peaks[retention_time] = peak


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

    start, retention_time, end = peaks[0][0], peaks[0][1], peaks[0][2]

    time = chromatogram.time
    idx = utils.indices_from_boundary(time, start, end)

    peak = Classes.Peak(retention_time, start, end, indices=idx)
    peak.get_integral(chromatogram)

    chromatogram.internal_standard = peak

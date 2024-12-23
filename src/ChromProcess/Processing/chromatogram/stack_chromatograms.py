"""
Functions for plotting a heat map of
a chromatogram series.
"""

import numpy as np
from scipy import interpolate
from ChromProcess.Classes import Chromatogram


def get_chrom_time_min_max(chromatograms: list[Chromatogram]) -> tuple[float, float]:
    """
    Get the highest and lowest retention times from a set of chromatograms.

    Parameters
    ----------
    chromatograms: list[ChromProcess.Classes.Chromatogram]

    Returns
    -------
    max_time, min_time: float
    """
    max_time = 1e100
    min_time = 0

    for c in chromatograms:
        if c.time.min() > min_time:
            min_time = c.time.min()
        if c.time.max() < max_time:
            max_time = c.time.max()

    return min_time, max_time


def stack_chromatograms(
    chromatograms: list[Chromatogram],
) -> tuple[np.ndarray, np.ndarray]:
    """
    Create a stack of chromatogram signals in a numpy array.

    Parameters
    ----------
    chromatograms: list[Chromatogram]

    Returns
    -------
    time_axis: np.ndnumpy
    chrom_stack: 2d numpy array (len(time_axis),len(chromatograms))
    """

    min_time, max_time = get_chrom_time_min_max(chromatograms)

    interpolation_length = len(chromatograms[0].time)

    chrom_stack = np.empty((len(chromatograms), interpolation_length))
    time_axis = np.linspace(min_time, max_time, num=interpolation_length)

    for c, chrom in enumerate(chromatograms):
        interp_function = interpolate.interp1d(chrom.time, chrom.signal)
        interpolated_signal = interp_function(time_axis)
        chrom_stack[c] = interpolated_signal

    return time_axis, chrom_stack

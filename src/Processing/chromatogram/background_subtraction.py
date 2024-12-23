from copy import copy
import numpy as np
from scipy.interpolate import interp1d
from ChromProcess.Classes import Chromatogram


def background_subtraction(
    chromatogram: Chromatogram, control_chromatogram: Chromatogram
) -> Chromatogram:
    """
    Subtract the signal of the control_chromatogram from a chromatogram.

    Parameters
    ----------
    chromatogram: Chromatogram
    control_chromatogram: Chromatogram

    Returns
    -------
    baseline_subtracted: Chromatogram
    """
    # Determine a common range between the two chromatograms
    lower = max(control_chromatogram.time.min(), chromatogram.time.min())
    upper = min(control_chromatogram.time.max(), chromatogram.time.max())

    region_mask = (chromatogram.time >= lower) & (chromatogram.time <= upper)
    new_time_axis = chromatogram.time[region_mask]
    sample_signal_axis = chromatogram.signal[region_mask]

    # Interpolate control_chromatogram so that it has the same timepoints as
    # chromatogram
    interp_f = interp1d(control_chromatogram.time, control_chromatogram.signal)
    interpolated_control_signal = interp_f(new_time_axis)

    # Subtract the interpolated control_chromatogram from the chromatogram
    baseline_subtracted = sample_signal_axis - interpolated_control_signal

    # return the result
    baseline_substracted_chromatogram = copy(chromatogram)

    baseline_substracted_chromatogram.time = new_time_axis
    baseline_substracted_chromatogram.signal = baseline_subtracted

    return baseline_substracted_chromatogram


def ic_background_subtraction(
    chromatogram: Chromatogram, threshold: int = 500
) -> np.ndarray:
    """
    Gets the ion chromatograms of the analysis and reconstitutes the total ion
    chromatogram ommitting m/z signals which do not exceed a threshold.

    Parameters
    ----------
    chromatogram: Chromatogram
        Chromatogram to be processed.
    threshold: float
        Ion chromatograms which do not exceed this threshold will be removed
        before the total ion chromatogram is reconsituted.

    Returns
    -------
    1D numpy array.
        Original signal if no mass spectral information is present in the
        chromatogram, processed signal if ms info is present.
    """

    if len(chromatogram.mz_intensity) == 0:
        return chromatogram.signal
    else:
        new_chromatogram = np.zeros(len(chromatogram.time))
        for s in range(0, len(chromatogram.point_counts)):
            start = chromatogram.scan_indices[s]
            end = start + chromatogram.point_counts[s]
            inten = chromatogram.mz_intensity[start:end]
            new_chromatogram[s] = np.sum(inten[inten > threshold])

        return new_chromatogram

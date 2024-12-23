import numpy as np
from .ion_chromatogram_from_region import ion_chromatogram_from_region
from ChromProcess.Classes import Peak
from ChromProcess.Classes import Chromatogram


def ion_chromatogram_from_peak(
    peak: Peak, parent_chromatogram: Chromatogram, threshold: float = 0.1
) -> dict[float, np.ndarray]:
    """
    Create a dictionary of ion chromatograms using information from the peak
    and its parent chromatogram.

    Parameters
    ----------
    peak: Peak
        Peak containing information.
    parent_chromatogram: Chromatogram
       Chromatogram containing information.
    threshold: float
        Threshold for mass spectra extraction relative to the maximum signal
        in the region.

    Returns
    -------
    ion_chromatograms: dict[mz, np.ndarray]
        A dictionary of ion chromatograms keyed by m/z value.
    """

    # Get the indices of the peaks data points in the chromatogram
    inds = peak.indices

    # find the relevant chromatogram attribute sections
    time = parent_chromatogram.time[inds]

    if len(time) != 0:
        ion_chromatograms = ion_chromatogram_from_region(
            parent_chromatogram,
            time.min(),
            time.max(),
            threshold=threshold,
        )
    else:
        ion_chromatograms = dict()

    return ion_chromatograms

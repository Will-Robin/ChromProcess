import numpy as np
from ChromProcess.Classes import Chromatogram


def ic_background_subtraction(
    chromatogram: Chromatogram, threshold: int = 500
) -> np.ndarray:
    """
    Gets the ion chromatograms of the analysis and reconstitutes the total ion
    chromatogram ommitting m/z signals which do not exceed a threshold.

    Parameters
    ----------
    chromatogram: ChromProcess.Classes.Chromatogram
        Chromatogram to be processed.
    threshold: float
        Ion chromatograms which do not exceed this threshold will be removed
        before the total ion chromatogram is reconsituted.

    Returns
    -------
    1D numpy array.
        Original signal if not mass spectra information present.
        processed signal if ms info present.
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

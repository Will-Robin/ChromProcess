import numpy as np

def ic_background_subtraction(chromatogram, threshold = 500):
    '''
    Gets the ion chromatograms of the analysis and reconstitutes the total ion
    chromatogram ommitting mass signals which do not exceed a threshold.

    Parameters
    ----------
    chromatogram: ChromProcess Chromatogram object.
        Chromatogram to be processed.
    threshold: float
        Ion chromatograms which do not exceed this threshold will be removed
        before the total ion chromatogram is reconsituted.

    Returns: 1D numpy array.
        original signal if not mass spectra information present.
        processed signal if ms info present.
    '''

    if len(chromatogram.mass_intensity) == 0:
        return chromatogram.signal
    else:
        inds = chromatogram.mass_intensity < threshold
        mass_intens = np.copy(chromatogram.mass_intensity)
        mass_intens[inds] = 0.0

        new_chromatogram = np.zeros(len(chromatogram.time))
        for s in range(0,len(chromatogram.point_counts)):

            start = chromatogram.scan_indices[s]
            end = start + chromatogram.point_counts[s]

            inten = mass_intens[start:end]

            new_chromatogram[s] = np.sum(inten)

        return new_chromatogram

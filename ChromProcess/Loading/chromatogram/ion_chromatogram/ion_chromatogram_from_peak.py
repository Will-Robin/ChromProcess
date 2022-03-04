from .ion_chromatogram_from_region import ion_chromatogram_from_region


def ion_chromatogram_from_peak(peak, parent_chromatogram, threshold=0.1):
    """
    Create a dictionary of ion chromatograms using information from the peak
    and its parent chromatogram.

    Parameters
    ----------
    peak: ChromProcess Peak object
        Peak containing information.
    parent_chromatogram: ChromProcess Chromatogram object
       Chromatogram containing information.
    spectrum_filter: float
        m/z intensities which do not exceed this fraction its parent
        mass spectrum will be omitted from the ion chromatogram.
    threshold: float
        Threshold for mass spectra extraction relative to the maximum signal
        in the region.

    Returns
    -------
    ion_chromatograms: dict[mz]:array()
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
        ion_chromatograms = {}

    return ion_chromatograms

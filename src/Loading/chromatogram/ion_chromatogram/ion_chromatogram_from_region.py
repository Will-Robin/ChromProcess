import numpy as np
from ChromProcess import Utils
from ChromProcess.Classes import Chromatogram


def ion_chromatogram_from_region(
    chromatogram: Chromatogram, lower: float, upper: float, threshold: float = 0.1
) -> dict[float, np.ndarray]:
    """
    Get the ion chromatograms for a region of a chromatogram. Requires mass
    spectra information to be present in the chromatogram.

    Parameters
    ----------
    chromatogram: Chromatogram
        Chromatogram containing information
    lower: float
        Lower bound for retention time region in chromatogram
    upper: float
        Upper bound for retention time region in chromatogram
    threshold: float
        Threshold for mass spectra extraction relative to the maximum signal
        in the region.

    Returns
    -------
    ion_chromatograms: dict[float, np.ndarray]
        dictionary of ion chromatograms
        {m/z: intensities over time}
    """

    if len(chromatogram.mz_values) > 0:
        # Select scan indices and point counts in region
        inds = Utils.indices_from_boundary(chromatogram.time, lower, upper)
        scan_inds = chromatogram.scan_indices[inds]
        p_counts = chromatogram.point_counts[inds]

        ion_chromatograms = dict()

        # iterate over mass spectra recorded at each time point
        for s in range(0, len(scan_inds)):
            lower_idx = scan_inds[s]
            upper_idx = scan_inds[s] + p_counts[s]
            slice_indices = np.arange(lower_idx, upper_idx, 1)
            slice_indices = slice_indices[
                chromatogram.mz_values[slice_indices] > threshold
            ]

            mz_values = chromatogram.mz_values[slice_indices]
            inten = chromatogram.mz_intensity[slice_indices]

            # Add the intensity values into the appropriate m/z channel
            for m in range(0, len(mz_values)):
                mz = np.round(mz_values[m])
                if mz not in ion_chromatograms:
                    ion_chromatograms[mz] = np.zeros(len(scan_inds))
                ion_chromatograms[mz][s] += inten[m]

    return ion_chromatograms

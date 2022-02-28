import numpy as np
from ChromProcess import Classes


def mass_spectrum_from_peak(peak, chromatogram):
    """
    Get the mass spectrum at the apex of a peak. Inserts the mass spectrum into
    the Peak object and returns the mass spectrum. If mass spectra information
    is not present in the chromatogram, an ewmpty list is returned and the
    Peak's mass_spectra value remains unchanged (False by default).

    Parameters
    ----------
    peak: ChromProcess Peak object
        Peak to find mass spectrum for.
    chromatogram: ChromProcess Chromatogram object
        Parent chromatogram for the peak

    Returns
    -------
    mass_spectrum: Classes.MassSpectrum
    """

    if len(chromatogram.mz_values) > 0:
        idx = np.where(chromatogram.time == peak.retention_time)[0]

        start = chromatogram.scan_indices[idx][0]
        end = start + chromatogram.point_counts[idx][0]

        mz_values = np.round(chromatogram.mz_values[start:end], 2)
        intensities = chromatogram.mz_intensity[start:end]
        v_idx = np.where(intensities > 0.0)[0]

        return Classes.MassSpectrum(
            mz_values[v_idx], intensities[v_idx], pos=peak.retention_time
        )

    else:
        return Classes.MassSpectrum([], [])

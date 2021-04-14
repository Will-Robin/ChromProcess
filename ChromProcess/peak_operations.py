def peakIonChromatogram(peak, chromatogram):

    '''
    Get the mass spectrum at the apex of a peak. Inserts the mass spectrum into
    the Peak object and returns the mass spectrum.
    If mass spectra information is not present in the chromatogram, an ewmpty
    list is returned and the Peak's mass_spectra value remains unchanged
    (False by default).

    Parameters
    ----------
    peak: ChromProcess Peak object
        Peak to find mass spectrum for.
    chromatogram: ChromProcess Chromatogram object
        Parent chromatogram for the peak

    Returns
    -------
    mass_spectrum: list
        list; [m/z, intensity] as numpy arrays
    '''
    import numpy as np

    if chromatogram.mass_spectra:
        idx = np.where(chromatogram.time == peak.retention_time)[0]

        start = chromatogram.scan_indices[idx][0]
        end = start + chromatogram.point_counts[idx][0]

        masses = np.round(chromatogram.mass_values[start:end],2)
        intensities = chromatogram.mass_intensity[start:end]
        v_idx = np.where(intensities > 0.0)[0]
        mass_spectrum = [masses[v_idx],intensities[v_idx]]

        peak.mass_spectrum = mass_spectrum

        return mass_spectrum

    else:
        return [np.array([]),np.array([])]

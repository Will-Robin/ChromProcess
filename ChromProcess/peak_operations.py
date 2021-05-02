def peakMassSpectrum(peak, chromatogram):

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


def peakIonChromatograms(peak, parent_chromatogram, spectrum_filter = 0.1):

    import numpy as np
    from ChromProcess import mass_spectra as ms

    inds = peak.indices
    
    time = parent_chromatogram.time[inds]
    signal = parent_chromatogram.signal[inds]
    scan_inds = parent_chromatogram.scan_indices[inds]
    p_counts = parent_chromatogram.point_counts[inds]

    for s in range(0,len(time)):

        inten = parent_chromatogram.mass_intensity[scan_inds[s]:scan_inds[s]+p_counts[s]]
        masses = parent_chromatogram.mass_values[scan_inds[s]:scan_inds[s]+p_counts[s]]

        if len(inten) > 0:

            filt_inds = np.where(inten > spectrum_filter*np.amax(inten))[0]
            inten = inten[filt_inds]
            masses = masses[filt_inds]

            round = np.round(masses, 2)

            for m in range(0,len(round)):
                if round[m] in peak.ion_chromatograms:
                    peak.ion_chromatograms[round[m]][s] = inten[m]
                else:
                    peak.ion_chromatograms[round[m]] = np.zeros(len(time))
                    peak.ion_chromatograms[round[m]][s] = inten[m]

        else:
            pass

    ms.bin_ion_chromatograms(peak, stdev = 0.1)

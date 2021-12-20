'''
I would prefer to add these functions into Peak object methods.

However, they rely upon calling mass spectral information from the
parent Chromatogram object.

The fact that any chromatogram can be provided as a source of the mass
spectra is a weakness: if any chromatogram can be provided, then
the wrong mass spectrum for the peak can be obtained.
'''

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


def peakIonChromatograms(peak, parent_chromatogram,
                         ic_filter = 0.05, spectrum_filter = 0.1):
    '''
    peak: ChromProcess Peak object
    parent_chromatogram: ChromProcess Chromatogram object
    ic_filter: float
        ion chromatograms which do not exceed this fraction of the
        maximum total ion chromatogram peak height will not be included.
    spectrum_filter: float
        m/z intensities which do not exceed this fraction its parent
        mass spectrum will be omitted from the ion chromatogram.

    modifies peak.ion_chromatograms dictionary in Peak object
    '''
    import numpy as np
    from ChromProcess.Processing import mass_spectra as ms

    # wipe previously stored ion chromatograms
    peak.ion_chromatograms = {}

    # Get the indices of the peaks data points in the chromatogram
    inds = peak.indices

    # find the relevant chromatogram attribute sections
    time = parent_chromatogram.time[inds]
    signal = parent_chromatogram.signal[inds]
    scan_inds = parent_chromatogram.scan_indices[inds]
    p_counts = parent_chromatogram.point_counts[inds]

    # iterate over mass spectra recorded at each time point
    for s in range(0,len(time)):
        # get mass spectrum at time point
        inten = parent_chromatogram.mass_intensity[scan_inds[s]:scan_inds[s]+p_counts[s]]
        mz_values = parent_chromatogram.mass_values[scan_inds[s]:scan_inds[s]+p_counts[s]]

        # filter out low intensity m/z signals
        filt_inds = np.where(inten > spectrum_filter*np.amax(inten))[0]

        if len(filt_inds) > 0:

            inten = inten[filt_inds]
            mz_values = mz_values[filt_inds]

            # round m/z to 2 d.p.
            rounded_mz = np.round(mz_values, 2)

            # add the intensity values into the appropriate m/z channel
            for m in range(0,len(rounded_mz)):
                if rounded_mz[m] in peak.ion_chromatograms:
                    peak.ion_chromatograms[rounded_mz[m]][s] = inten[m]
                else:
                    peak.ion_chromatograms[rounded_mz[m]] = np.zeros(len(time))
                    peak.ion_chromatograms[rounded_mz[m]][s] = inten[m]

        else:
            pass

    # combine channels with close enough m/z value (given stdev)
    ms.bin_ion_chromatograms(peak, stdev = 0.1)

    # remove ion chromatograms whose maximum values do not exceed the
    # fraction of the total ion chormatogram defind by ic_filter
    threshold = ic_filter*signal.max()
    remove_ic = []
    for ic in peak.ion_chromatograms:
        if peak.ion_chromatograms[ic].max() < threshold:
            remove_ic.append(ic)

    for i in remove_ic:
        del peak.ion_chromatograms[i]

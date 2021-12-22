import numpy as np

def ion_chromatogram_from_region(chromatogram, lower, upper, threshold = 0.1):
    '''
    Get the ion chromatograms for a region of a chromatogram. Requires mass
    spectra information to be present in the chromatogram.

    Parameters
    ----------
    chromatogram: ChromProcess Chromatogram object
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

    ion_dict: dict
        dictionary of ion chromatograms
        {m/z: intensities over time}
    '''

    if not chromatogram.mass_spectra:
        return None
    else:
        ion_dict = {}
        inds = np.where((chromatogram.time > lower)&(chromatogram.time < upper))[0]

        time = chromatogram.time[inds]
        scan_inds = chromatogram.scan_indices[inds]
        p_counts = chromatogram.point_counts[inds]

        for s in range(0,len(time)):

            inten = chromatogram.mass_intensity[scan_inds[s]:scan_inds[s]+p_counts[s]]
            masses = chromatogram.mass_values[scan_inds[s]:scan_inds[s]+p_counts[s]]

            if len(inten) > 0:

                filt_inds = np.where(inten > threshold*np.amax(inten))[0]
                inten = inten[filt_inds]
                masses = masses[filt_inds]

                round = np.round(masses, 2)

                for m in range(0,len(round)):
                    if round[m] in ion_dict:
                        ion_dict[round[m]][s] = inten[m]
                    else:
                        ion_dict[round[m]] = np.zeros(len(time))
                        ion_dict[round[m]][s] = inten[m]
            else:
                pass

    return ion_dict


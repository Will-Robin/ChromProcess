import numpy as np

'''
Functions for dealing with mass spectra.
'''
def ion_chromatogram_from_region(chromatogram, lower, upper, threshold = 0.1):
    '''
    Get the ion chromatograms from a region of a chromatogram.

    Parameters
    ----------
    chromatogram: ChromProcess Chromatogram object
        Chromatogram containing mass spectral information.
    lower: float
        Lower retention time bound of the region.
    upper: float
        Upper retention time bound of the region

    threshold: float
        Intensities below this fraction of the total ion chromatogram peak
        singal will not be included in the output.
    
    Returns
    -------
    time: numpy array
        A time axis for the data
    ion_dict: dict
        A dictionary containing numpy arrays indexes my m/z values.
    '''

    time = np.array([])
    ion_dict = {}

    if len(chromatogram.scan_indices) == 0:
        pass
    else:
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

    return time, ion_dict
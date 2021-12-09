import os
import numpy as np

from ChromProcess import simple_functions as s_f

'''
Functions for dealing with mass spectra.
'''

def ion_chromatogram_region(chromatogram, lower, upper, threshold = 0.1):
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

    if len(chromatogram.scan_indices) == 0:
        pass
    else:
        ion_dict = {}
        inds = np.where((chromatogram.time > lower)&(chromatogram.time < upper))[0]

        time = chromatogram.time[inds]
        signal = chromatogram.signal[inds]
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

def get_peak_mass_spectra(series):
    for c in series.chromatograms:
        if len(c.scan_indices) == 0:
            pass
        else:
            for p in c.peaks:
                ind = np.where(c.time == c.peaks[p].retention_time)[0]

                start = c.scan_indices[ind][0]
                end = start + c.point_counts[ind][0]

                c.peaks[p].mass_spectrum = [np.round(c.mass_values[start:end],2), c.mass_intensity[start:end]]


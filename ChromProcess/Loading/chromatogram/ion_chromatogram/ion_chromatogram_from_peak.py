import numpy as np
from ChromProcess.Utils.utils import utils as ut

def ion_chromatogram_from_peak(
                                peak, 
                                parent_chromatogram,
                                ic_filter = 0.05, 
                                spectrum_filter = 0.1
                                ):
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
    ut.bin_dictionary(peak, stdev = 0.1)

    # remove ion chromatograms whose maximum values do not exceed the
    # fraction of the total ion chormatogram defind by ic_filter
    threshold = ic_filter*signal.max()
    remove_ic = []
    for ic in peak.ion_chromatograms:
        if peak.ion_chromatograms[ic].max() < threshold:
            remove_ic.append(ic)

    for i in remove_ic:
        del peak.ion_chromatograms[i]

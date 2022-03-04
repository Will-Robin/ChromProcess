import numpy as np
from ChromProcess.Utils.utils import utils


def ion_chromatogram_from_region(chromatogram, lower, upper, threshold=0.1):
    """
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
    spectrum_filter: float
        m/z intensities which do not exceed this fraction its parent
        mass spectrum will be omitted from the ion chromatogram.
    threshold: float
        Threshold for mass spectra extraction relative to the maximum signal
        in the region.

    Returns
    -------
    ion_chromatograms: dict
        dictionary of ion chromatograms
        {m/z: intensities over time}
    """

    mz_regions = {}
    if len(chromatogram.mz_values) > 0:
        inds = utils.indices_from_boundary(chromatogram.time, lower, upper)

        time = chromatogram.time[inds]
        scan_inds = chromatogram.scan_indices[inds]
        p_counts = chromatogram.point_counts[inds]

        # iterate over mass spectra recorded at each time point
        for s in range(0, len(time)):
            # get mass spectrum at time point
            inten = chromatogram.mz_intensity[scan_inds[s] : scan_inds[s] + p_counts[s]]
            mz_values = chromatogram.mz_values[
                scan_inds[s] : scan_inds[s] + p_counts[s]
            ]

            # filter out low intensity m/z signals
            filt_inds = np.where(inten > threshold * np.amax(inten))[0]

            if len(filt_inds) > 0:

                inten = inten[filt_inds]
                masses = mz_values[filt_inds]

                round = np.round(masses, 2)

                # add the intensity values into the appropriate m/z channel
                for m in range(0, len(round)):
                    if round[m] in mz_regions:
                        mz_regions[round[m]][s] = inten[m]
                    else:
                        mz_regions[round[m]] = np.zeros(len(time))
                        mz_regions[round[m]][s] = inten[m]
            else:
                pass

    # combine channels with close enough m/z value (given stdev)
    ion_chromatograms = utils.bin_dictionary(mz_regions, stdev=0.1)

    # remove ion chromatograms whose maximum values do not exceed the
    # fraction of the total ion chromatogram defind by threshold
    threshold = threshold * chromatogram.signal.max()
    remove_ic = []
    for ic in ion_chromatograms:
        if ion_chromatograms[ic].max() < threshold:
            remove_ic.append(ic)

    for i in remove_ic:
        del ion_chromatograms[i]

    return ion_chromatograms

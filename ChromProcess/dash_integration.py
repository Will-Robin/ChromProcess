def load_cdf_file(fname, import_ms = False):
    import numpy as np
    from ChromProcess import file_import
    data = {}
    data['time'] = file_import.get_data_cdf_GCMS(fname,
                            'scan_acquisition_time')/60 # converted to minutes
    data['signal'] = file_import.get_data_cdf_GCMS(fname,
                                                'total_intensity')

    data['mass_spectra'] = {}

    load_list = []
    if import_ms:
        data['mass_spectra']['intensities'] = file_import.get_data_cdf_GCMS(fname, "intensity_values")
        # Measured masses
        data['mass_spectra']['m/z'] = np.round(file_import.get_data_cdf_GCMS(fname, "mass_values"),3)
        # Scan index is starting index
        data['mass_spectra']['scan_indices'] = file_import.get_data_cdf_GCMS(fname, "scan_index")
        # Point count is number of elements to read
        data['mass_spectra']['point_counts'] = file_import.get_data_cdf_GCMS(fname, "point_count")

    return data

def pickPeaksRegion(time,signal, region, threshold = 0.1):
    '''
    Pick peaks in a region of a chromatogram.

    Parameters
    ----------
    chromatogram: ChromProcess Chromatogram object
        Chromatogram containing information
    region: list
        [lower bound, upper bound] for region

    threshold: float
        Threshold for peak picking relative to the maximum signal
        in the region.
    Returns
    -------
    None
    '''
    import numpy as np
    from ChromProcess import Classes
    from ChromProcess import chromatogram_operations as chrom_ops

    inds = np.where( (time > region[0]) &
                     (time < region[1]) )[0]

    reg_time = time[inds]
    reg_signal = signal[inds]

    rts, peak_times = chrom_ops.pickPeaks(reg_time, reg_signal,
                                          threshold = threshold)
    peaks_idx = {'indices':[]}
    for rt, times in zip(rts,peak_times):

        idx = np.where((time >= times[0]) &
                        (time <= times[-1]))[0]

        peaks_idx['indices'].append(idx)

    return peaks_idx

def mz_background_subtraction(data, threshold = 500):
    '''
    Gets the ion chromatograms of the analysis and reconstitutes the total ion
    chromatogram ommitting mass signals which do not exceed a threshold.

    Parameters
    ----------
    chromatogram: ChromProcess Chromatogram object.
        Chromatogram to be processed.
    threshold: float
        Ion chromatograms which do not exceed this threshold will be removed
        before the total ion chromatogram is reconsituted.

    Returns: None
        Modifies the chromatogram in-place.
    '''
    import numpy as np

    if 'mass_spectra' not in data:
        return data['signal']
    else:
        intensities = np.array(data['mass_spectra']['intensities'])
        inds = intensities < threshold
        intensities[inds] = 0.0

        new_chromatogram = np.zeros(len(data['time']))
        for s in range(0,len(data['mass_spectra']['point_counts'])):

            start = data['mass_spectra']['scan_indices'][s]
            end = start + data['mass_spectra']['point_counts'][s]

            inten = intensities[start:end]

            new_chromatogram[s] = np.sum(inten)
        print(new_chromatogram, 'haha')
        return new_chromatogram

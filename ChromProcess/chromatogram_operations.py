def pickPeaks(time, signal, threshold = 0.1):
    '''
    Finds peaks in signal and returns their times and indices

    Parameters
    ----------
    time, signal: 1D arrays
        Signal from which peaks will be picked and its corresponding time axis.
    threshold: float
        threshold for peak picking as a proportion of the maximum value of the
        signal.

    Returns
    -------
    retention_times: list
        list of peak retention times (floats)

    peak_windows:
        Time windows for peaks
    '''
    import numpy as np
    from ChromProcess import processing_functions as p_f

    peaks = p_f.Peak_finder(signal, thres = threshold, min_dist = 0.1,
                            max_inten = 1e100, min_inten = -1e100)

    retention_times = []
    peak_windows = []

    for p in range(0,len(peaks["Peak_indices"])):

        peak_idx = peaks["Peak_indices"][p]
        start_idx  = peaks["Peak_start_indices"][p]
        end_idx = peaks["Peak_end_indices"][p]

        rt = time[peak_idx]
        peak_inds = np.where( (time >= time[start_idx]) &
                              (time <= time[end_idx]) )[0]

        if len(peak_inds) > 0:
            retention_times.append(rt)
            peak_windows.append(time[peak_inds])


    return retention_times, peak_windows

def pickPeaksRegion(chromatogram, region, threshold = 0.1):
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

    inds = np.where( (chromatogram.time > region[0]) &
                     (chromatogram.time < region[1]) )[0]

    time = chromatogram.time[inds]
    signal = chromatogram.signal[inds]

    rts, peak_times = chrom_ops.pickPeaks(time, signal, threshold = threshold)

    for rt, times in zip(rts,peak_times):

        idx = np.where((chromatogram.time >= times[0]) &
                        (chromatogram.time <= times[-1]))[0]

        peak = Classes.Peak(rt, idx)
        peak.get_integral(chromatogram)
        chromatogram.peaks[rt] = peak

def mz_background_subtraction(chromatogram, threshold = 500):
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

    Returns: 1D numpy array.
        original signal if not mass spectra information present.
        processed signal if ms info present.
    '''
    import numpy as np

    if len(chromatogram.mass_intensity) == 0:
        return chromatogram.signal
    else:
        inds = chromatogram.mass_intensity < threshold
        mass_intens = np.copy(chromatogram.mass_intensity)
        mass_intens[inds] = 0.0

        new_chromatogram = np.zeros(len(chromatogram.time))
        for s in range(0,len(chromatogram.point_counts)):

            start = chromatogram.scan_indices[s]
            end = start + chromatogram.point_counts[s]

            inten = mass_intens[start:end]

            new_chromatogram[s] = np.sum(inten)

        return new_chromatogram

def get_mz_background(chromatogram, threshold = 500):
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

    if len(chromatogram.mass_intensity) == 0:
        return print("MS_intensity_threshold_chromatogram: no mass spectra information in chromatogram")
        return chromatogram.signal
    else:
        inds = chromatogram.mass_intensity > threshold
        mass_intens = np.copy(chromatogram.mass_intensity)
        mass_intens[inds] = 0.0

        new_chromatogram = np.zeros(len(chromatogram.time))
        for s in range(0,len(chromatogram.point_counts)):

            start = chromatogram.scan_indices[s]
            end = start + chromatogram.point_counts[s]

            inten = mass_intens[start:end]

            new_chromatogram[s] = np.sum(inten)

        return new_chromatogram

def internalRefIntegral(chromatogram, internal_ref_region):
    '''
    Get integrals of internal references for a Chromatogram_Series
    Parameters
    ----------
    series: Chromatogram_Series object
        Object containing chromatograms and associated series data which is
        modified by the function.
    '''
    from ChromProcess import Classes
    from ChromProcess import chromatogram_operations as chrom_ops
    if len(chromatogram.peaks) == 0:
        chrom_ops.pickPeaksRegion(chromatogram, internal_ref_region, threshold = 0.1)
    else:
        pass

    for p in chromatogram.peaks:
        if p > internal_ref_region[0] and p < internal_ref_region[1]:
            rt = chromatogram.peaks[p].retention_time
            indices = chromatogram.peaks[p].indices
            chromatogram.internal_reference = Classes.Peak(rt, indices)
            chromatogram.internal_reference.get_integral(chromatogram)

def getIonChromatogramsFromRegion(chromatogram, lower, upper, threshold = 0.1):
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

    import numpy as np

    if not chromatogram.mass_spectra:
        return None
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

    return ion_dict

def RegionSVD(chromatogram, region, ic_threshold = 0.1,
                          components = 5):

    '''
    Reduce a GC-MS chromatogram region dimensionality using singular value
    decomposition. The chromatogram must contain ion chromatogram information.

    Parameters
    ----------
    chromatogram: ChromProcess Chromatogram object
        Chromatogram to work on.
    region: list, tuple or array
        Two-element list of floats defining the lower and upper bounds of the
        region.

    Returns
    -------

    U: numpy array
        Unitary matrix ith left singular vectors

    S: numpy array
        array of singular values

    Vh: numpy array
        Unitary matrix with right singular vectors


    '''
    import numpy as np
    from ChromProcess import chromatogram_operations as chrom_ops
    from ChromProcess import simple_functions as simp_func

    ic_dict = chrom_ops.getIonChromatogramsFromRegion(chromatogram,
                                                      region[0], region[1],
                                                      threshold = ic_threshold)
    if ic_dict == None:
        return None
    else:
        # convert dict into nump arrays
        mass_ax = np.array([*ic_dict])
        ic_stack = np.array(list(ic_dict.values()))
        trans_IC = ic_stack.T

        samples, features = trans_IC.shape

        U, S, Vh = simp_func.runSVD(trans_IC)

        exp_variance = (S**2)/(samples - 1)

        # inverse transform
        #IT = np.dot(U*S, Vh)

        return U[:,:components], S[:components], Vh[:components]

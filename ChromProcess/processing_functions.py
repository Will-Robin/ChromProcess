import numpy as np
from ChromProcess import processing_functions as p_f

'''
Basic functions to be applied to chromatography data.
'''

def MS_intensity_threshold_chromatogram(chromatogram, threshold = 500):
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

    if len(chromatogram.mass_intensity) == 0:
        return print("MS_intensity_threshold_chromatogram: no mass spectra information in chromatogram")
    else:
        inds = chromatogram.mass_intensity < threshold
        chromatogram.mass_intensity[inds] = 0.0
        new_chromatogram = np.zeros(len(chromatogram.time))

        for s in range(0,len(chromatogram.point_counts)):

            start = chromatogram.scan_indices[s]
            end = start + chromatogram.point_counts[s]

            #mass = ch.mass_values[start:end]
            inten = chromatogram.mass_intensity[start:end]

            new_chromatogram[s] = np.sum(inten)

        chromatogram.signal = new_chromatogram

def Peak_finder(intensity, thres=0.1, min_dist=1, max_inten = 1e100, min_inten = -1e100):

    '''Peak detection routine.
    Modified from PeakUtils:
    https://github.com/atjacobs/PeakUtils/tree/master/peakutils
    Finds the peaks in *y* by taking its first order difference. By using
    *thres* and *min_dist* parameters, it is possible to reduce the number of
    detected peaks.

    Parameters
    ----------
    y : ndarray
        1D amplitude data to search for peaks.
    thres : float between [0., 1.]
        Normalized threshold. Only the peaks with amplitude higher than the
        threshold will be detected.
    min_dist : int
        Minimum distance between each detected peak. The peak with the highest
        amplitude is preferred to satisfy this constraint.
    max_inten: float
        peaks will not be detected above this threshold.
    min_inten: float
        peaks will not be detected below this threshold.

    Returns
    -------
    ndarray
        Array containing the indexes of the peaks that were detected
    '''

    test = np.where((intensity < max_inten))[0]
    thres *= np.max(intensity[test]) - np.min(intensity[test])

    # find the peaks by using the first order difference
    diff = p_f.savitzky_golay(np.diff(intensity), 7, 3, deriv=0, rate=1)
    
    peaks_indices = np.where((np.hstack([diff, 0.]) < 0.)
                     & (np.hstack([0., diff]) > 0.)
                     & (intensity > thres) & (intensity < max_inten) & (intensity > min_inten))[0]

    if peaks_indices.size > 1 and min_dist > 1:
        highest = peaks_indices[np.argsort(intensity[peaks_indices])][::-1]
        rem = np.ones(intensity.size, dtype=bool)
        rem[peaks_indices] = False

        for peak in highest:
            if not rem[peak]:
                sl = slice(max(0, peak - min_dist), peak + min_dist + 1)
                rem[sl] = True
                rem[peaks_indices] = False

        peaks_indices = np.arange(intensity.size)[~rem]

    # find beginning of peaks
    peak_starts = []
    for n in range(0,len(peaks_indices)):
        cursor = peaks_indices[n]-2
        while diff[cursor] > 0:
            cursor -= 1
        peak_starts.append(cursor)

    # find end of peaks
    peak_ends = []
    for n in range(0,len(peaks_indices)):
        cursor = peaks_indices[n]+2
        if int(cursor) >= len(diff):
            cursor = len(diff)-1

        while diff[cursor] < 0:
            cursor += 1
            if int(cursor) >= len(diff):
                cursor = len(diff)-1
                break
        peak_ends.append(cursor)

    return {'Peak_indices':peaks_indices, 'Peak_start_indices':peak_starts, 'Peak_end_indices':peak_ends}

def name_peak(peak_rt,bound_dict):
    """
    Takes a peak (retention times of peak) and assigns peak name based
    on a dictionary of boundaries. The assignment priority is based on the
    interation order of the dict (Python 3 dict)

    Parameters
    ----------
    peak_rt: float
        retention time of a peak
    bound_dict: dict
        dictionary of boundaries for peak assignments
    Returns
    -------
    peak_name: str
        name of the peak. If the peak has no assignment, the peak retention time
        is returned
    """

    # Assign peak name
    peak_name = str(round(peak_rt,3))
    for b in bound_dict:
        if bound_dict[b][0] < peak_rt < bound_dict[b][1]:
            peak_name = b
            break

    return peak_name

def int_test(x):
    '''
    Test if variable can be converted to an integer.

    Parameters
    ----------
    x: any type

    Returns
    -------
    bool
    '''
    try:
        int(x)
        return True
    except ValueError:
            return False

def savitzky_golay(y, window_size, order, deriv=0, rate=1):
    """
    From https://scipy-cookbook.readthedocs.io/items/SavitzkyGolay.html

    Smooth (and optionally differentiate) data with a Savitzky-Golay filter.
    From the Scipy Cookbook
    The Savitzky-Golay filter removes high frequency noise from data.
    It has the advantage of preserving the original shape and
    features of the signal better than other types of filtering
    approaches, such as moving averages techniques.
    Parameters
    ----------
    y : array_like, shape (N,)
        the values of the time history of the signal.
    window_size : int
        the length of the window. Must be an odd integer number.
    order : int
        the order of the polynomial used in the filtering.
        Must be less then `window_size` - 1.
    deriv: int
        the order of the derivative to compute (default = 0 means only smoothing)
    Returns
    -------
    ys : ndarray, shape (N)
        the smoothed signal (or it's n-th derivative).
    Notes
    -----
    The Savitzky-Golay is a type of low-pass filter, particularly
    suited for smoothing noisy data. The main idea behind this
    approach is to make for each point a least-square fit with a
    polynomial of high order over a odd-sized window centered at
    the point.
    Examples
    --------
    t = np.linspace(-4, 4, 500)
    y = np.exp( -t**2 ) + np.random.normal(0, 0.05, t.shape)
    ysg = savitzky_golay(y, window_size=31, order=4)
    import matplotlib.pyplot as plt
    plt.plot(t, y, label='Noisy signal')
    plt.plot(t, np.exp(-t**2), 'k', lw=1.5, label='Original signal')
    plt.plot(t, ysg, 'r', label='Filtered signal')
    plt.legend()
    plt.show()
    References
    ----------
    .. [1] A. Savitzky, M. J. E. Golay, Smoothing and Differentiation of
       Data by Simplified Least Squares Procedures. Analytical
       Chemistry, 1964, 36 (8), pp 1627-1639.
    .. [2] Numerical Recipes 3rd Edition: The Art of Scientific Computing
       W.H. Press, S.A. Teukolsky, W.T. Vetterling, B.P. Flannery
       Cambridge University Press ISBN-13: 9780521880688
    """
    import numpy as np
    from math import factorial

    try:
        window_size = np.abs(np.int(window_size))
        order = np.abs(np.int(order))
    except ValueError:
        raise ValueError("window_size and order have to be of type int")
    if window_size % 2 != 1 or window_size < 1:
        raise TypeError("window_size size must be a positive odd number")
    if window_size < order + 2:
        raise TypeError("window_size is too small for the polynomials order")
    order_range = range(order+1)
    half_window = (window_size -1) // 2
    # precompute coefficients
    b = np.mat([[k**i for i in order_range] for k in range(-half_window, half_window+1)])
    m = np.linalg.pinv(b).A[deriv] * rate**deriv * factorial(deriv)
    # pad the signal at the extremes with
    # values taken from the signal itself
    firstvals = y[0] - np.abs( y[1:half_window+1][::-1] - y[0] )
    lastvals = y[-1] + np.abs(y[-half_window-1:-1][::-1] - y[-1])
    y = np.concatenate((firstvals, y, lastvals))

    return np.convolve( m[::-1], y, mode='valid')

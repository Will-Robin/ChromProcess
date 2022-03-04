import numpy as np
from scipy.optimize import curve_fit
from ChromProcess.Utils.signal_processing import deconvolution as d_c


def _1gaussian(x, amp1, cen1, sigma1):
    """
    A single gaussian function

    Parameters
    ----------
    x: array
        x axis data
    amp1: float
        amplitude of the function
    cen1: float
        centre of the function (mean)
    sigma1: float
        width of the function (standard deviation)

    Returns
    -------
    function: numpy array
        y values for the function
    """
    return (
        amp1
        * (1 / (sigma1 * (np.sqrt(2 * np.pi))))
        * (np.exp(-((x - cen1) ** 2) / ((2 * sigma1) ** 2)))
    )


def _2gaussian(x, amp1, cen1, sigma1, amp2, cen2, sigma2):
    """
    A double gaussian function

    Parameters
    ----------
    x: array
        x axis data
    ampn: float
        amplitude of a component gaussian function
    cenn: float
        centre of a component gaussian function (mean)
    sigman: float
        width of a component gaussian function (standard deviation)

    Returns
    -------
    function: numpy array
        y values for the function
    """
    return d_c._1gaussian(x, amp1, cen1, sigma1) + d_c._1gaussian(x, amp2, cen2, sigma2)


def _3gaussian(x, amp1, cen1, sigma1, amp2, cen2, sigma2, amp3, cen3, sigma3):
    return d_c._1gaussian(x, amp1, cen1, sigma1) + d_c._2gaussian(
        x, amp2, cen2, sigma2, amp3, cen3, sigma3
    )


def fit_gaussian_peaks(
    time,
    sig,
    peaks,
    initial_guess=[10000, 1, 0.005],
    lowerbounds=[0, 0, 0.0],
    upperbounds=[1e100, 1, 0.025],
):
    """
    TODO: This kind of function could be useful, but a better adapted function
        for peak deconvolution should be written. The function could take
        similar concepts to this one, but with a different concept for its
        implementation.

    Fitting sums of gaussian peaks to data using supplied peak indices.

    Parameters
    ----------
    time: array
        time values
    sig: array
        signal to be fit to
    peaks: list of peak positions
        list peak positions in the time

    initial_guess: list
        Initial guess for the peak amplitude, position and width
        e.g. see _1gaussian() function arguments.

    lowerbounds: list
        Lower bounds for the peak amplitude, position and width
        e.g. see _1gaussian() function arguments.
    upperbounds: list
        Upper bounds for the peak amplitude, position and width
        e.g. see _1gaussian() function arguments.

    Returns
    -------
    popt: ndarray
        list of fitted values [[amplitude, centre, width],]
    pcov: array
        correlation matrix
    """

    guess = []  # amp, cen, sig
    lbds = []
    ubds = []

    for p in range(0, len(peaks)):  # extend the boundary
        initial_guess[0] = np.amax(sig)
        initial_guess[1] = peaks[p]
        lowerbounds[1] = time[0]
        upperbounds[1] = time[-1]
        guess.extend(initial_guess)
        lbds.extend(lowerbounds)
        ubds.extend(upperbounds)

    boundarr = [lbds, ubds]
    if len(peaks) == 1:
        popt, pcov = curve_fit(d_c._1gaussian, time, sig, p0=guess, bounds=boundarr)
    elif len(peaks) == 2:
        popt, pcov = curve_fit(d_c._2gaussian, time, sig, p0=guess, bounds=boundarr)
    elif len(peaks) == 3:
        popt, pcov = curve_fit(d_c._3gaussian, time, sig, p0=guess, bounds=boundarr)
    else:
        print("Error fitting peaks")
        popt, pcov = [0, 0, 0], 0

    return popt, pcov


def deconvolute_region(chromatogram, region, num_peaks=1):

    """
    TODO: Combine the ideas in this function with fit_gaussian_peaks()

    Parameters
    ----------
    chromatogram: ChromProcess Chromatogram object
        Chromatogram
    region: list
        region of chromatogram under operation [lower bound, upper bound]

    Returns
    -------
    popt: ndarray
        list of fitted values [[amplitude, centre, width],]
    pcov: array
        correlation matrix
    """

    upper = region[1]
    lower = region[0]

    inds = np.where((chromatogram.time > lower) & (chromatogram.time < upper))[0]

    time = chromatogram.time[inds]
    signal = chromatogram.signal[inds]

    signal = signal - np.average(signal[-5:-1])

    peak_list = np.array([*chromatogram.peaks])
    peak_inds = np.where((peak_list > lower) & (peak_list < upper))[0]

    peaks = peak_list[peak_inds]

    while len(peaks) < num_peaks:
        peaks = np.append(peaks, np.average(peaks))

    if len(peaks) > num_peaks:
        peaks = peaks[:num_peaks]

    return d_c.fit_gaussian_peaks(time, signal, peaks)


def deconvolute_peak(peak, chromatogram, num_peaks=2):

    """
    TODO: this function is quite similar in scope to deconvolute_region().
    Refactor with the other two deconvolution macros.

    Parameters
    ----------
    chromatogram: ChromProcess Chromatogram object
        Chromatogram
    region: list
        region of chromatogram under operation [lower bound, upper bound]

    Returns
    -------
    popt: ndarray
        list of fitted values [[amplitude, centre, width],]
    pcov: array
        correlation matrix
    """

    peaks = [peak.retention_time for _ in range(num_peaks)]
    time = chromatogram.time[peak.indices]
    signal = chromatogram.signal[peak.indices]
    baseline = np.interp(time, [time[0], time[-1]], [signal[0], signal[-1]])

    signal = signal - baseline

    return d_c.fit_gaussian_peaks(time, signal, peaks)

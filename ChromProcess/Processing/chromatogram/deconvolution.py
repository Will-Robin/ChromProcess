import numpy as np
from scipy.stats import norm
from ChromProcess.Classes import Peak
from ChromProcess.Classes import Chromatogram
from ChromProcess.Classes import Deconvolution
from ChromProcess.Utils import deconvolution


def deconvolute_region(
    chromatogram: Chromatogram,
    region: list[float],
    num_peaks: int = 1,
    baseline_subtract: bool = False,
) -> list[Peak]:
    """
    Deconvolute a region of a chromatogram.

    The fitting procedure assumes that the baselinne in the selected region of
    the chromatogram is constant.

    Parameters
    ----------
    chromatogram: Chromatogram
        Chromatogram
    region: list
        region of chromatogram under operation [lower bound, upper bound]
    num_peaks: int
        Number of peaks expected in the region.
    baseline_subtract: bool
        Whether to subtract the fitted baseline from the deconvoluted peaks

    Returns
    -------
    new_peaks: list[Peak]
        List of deconvoluted peaks.
    """

    ## STAGE 1 preprocessing
    lower, upper = (region[0], region[1])

    # Cut out portion of signal
    idx = np.where((chromatogram.time > lower) & (chromatogram.time < upper))[0]

    time = chromatogram.time[idx]
    signal = chromatogram.signal[idx]

    # Find any peaks that have already been picked in the region.
    relevant_peaks = [
        chromatogram.peaks[p] for p in chromatogram.peaks if lower <= p <= upper
    ]

    peaks = [p.retention_time for p in relevant_peaks]

    if len(peaks) > num_peaks:
        print("Deconvolution warning: deconvoluting more peaks than specified.")

    height_init = np.array([p.height for p in relevant_peaks])

    width_init = np.array(
        [np.std(chromatogram.time[p.indices], ddof=1) for p in relevant_peaks]
    )

    # Pad out any extra peaks if more peaks are specified as arguments than are
    # found in the peak picking operation
    if len(peaks) < num_peaks:
        pad_len = num_peaks - len(peaks)

        peaks = np.hstack((peaks, np.linspace(time.min(), time.max(), pad_len)))

        width_init = np.pad(
            width_init,
            (0, pad_len),
            mode="constant",
            constant_values=(0.0, width_init.mean()),
        )

        height_init = np.pad(
            height_init,
            (0, pad_len),
            mode="constant",
            constant_values=(0.0, signal.mean()),
        )

    # Scale heights
    height_init /= len(height_init)

    baseline = signal.min()

    ## STAGE 2: perform deconvolution.
    popt, pcov = deconvolution.fit_pdf(
        time,
        signal,
        peaks,
        height_init,
        width_init,
        baseline,
    )

    # Format output into Peaks
    # -> [[magnitude, position, width, baseline],]
    result = np.vstack((np.reshape(popt[:-1], (3, -1)), np.full(len(peaks), popt[-1])))

    result = result.T

    # Create peak objects from deconvolution results
    new_peaks = []
    for x in range(0, len(result)):
        peak_features = result[x]
        magnitude = peak_features[0]
        retention_time = peak_features[1]
        width = peak_features[2]
        baseline = peak_features[3]

        start = retention_time - 4 * width  # 4σ
        end = retention_time + 4 * width  # 4σ
        indices = np.where((chromatogram.time >= start) & (chromatogram.time <= end))[0]

        # Calculate integral numerically (same as 'real' peaks).
        time_axis = chromatogram.time[indices]

        deconvoluted_signal = (
            magnitude * norm.pdf(time_axis, loc=retention_time, scale=width) + baseline
        )

        if baseline_subtract:
            deconvoluted_signal -= baseline

        integral = np.trapz(deconvoluted_signal, x=time_axis)

        deconv_params = Deconvolution(
            {
                "magnitude": magnitude,
                "loc": retention_time,
                "scale": width,
                "baseline": baseline,
            },
            region,
            "scipy.norm.pdf() + baseline",
        )

        new_peaks.append(
            Peak(
                retention_time,
                start,
                end,
                indices=indices,
                height=magnitude,
                integral=integral,
                deconvolution_params=deconv_params,
            )
        )

    return new_peaks

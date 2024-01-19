from ChromProcess.Classes import Peak
from ChromProcess.Classes import Chromatogram
from ChromProcess.Utils.utils import utils
from ChromProcess.Processing.peak import operations
from ChromProcess.Utils.peak_finding import pick_peaks


def find_peaks_in_region(
    chromatogram: Chromatogram,
    start: float,
    end: float,
    threshold: float = 0.1,
    smooth_width: float = 10,
    wlen: int = 10,
    distance: int = 1,
    prominence: float = 0.1,
    peak_local_baseline: bool = False,
) -> list[Peak]:
    """
    Find peaks within the chromatogram between start and end retention times.

    Parameters
    ----------
    chromatogram: Chromatogram
    start: float
        Start of time region (< end).
    end: float
        End of time region (> start).
    threshold: float
        Peaks below this fraction of the highest intensity of the chromatogram
        will not be picked.
    smooth_width: int
        Width for smoothing window.
    wlen: int
        Size of peak picking window in indices.
    distance: int
    prominence: float
    peak_local_baseline: bool

    Returns
    -------
    peaks: list[Peak]
        list of Peak objects
    """

    if start > end:
        print(f"peak start ({start}) > peak end, ({end}) returning None")
        return []

    inds = utils.indices_from_boundary(chromatogram.time, start, end)

    time = chromatogram.time[inds]
    signal = chromatogram.signal[inds]

    picked_peaks = pick_peaks(
        signal,
        smooth_width=smooth_width,
        distance=distance,
        threshold=threshold,
        prominence=prominence,
        wlen=wlen,
    )

    peaks = []
    for x in range(0, len(picked_peaks["Peak_indices"])):
        pk_idx = picked_peaks["Peak_indices"][x]
        start_idx = picked_peaks["Peak_start_indices"][x]
        end_idx = picked_peaks["Peak_end_indices"][x]

        retention_time = time[pk_idx]
        start = time[start_idx]
        end = time[end_idx]

        p_inds = utils.indices_from_boundary(chromatogram.time, start, end)

        new_peak = Peak(retention_time, start, end, indices=p_inds)
        new_peak.set_height(
            operations.get_peak_height(
                new_peak, chromatogram, baseline_subtract=peak_local_baseline
            )
        )
        new_peak.set_integral(
            operations.get_peak_integral(
                new_peak, chromatogram, baseline_subtract=peak_local_baseline
            )
        )

        peaks.append(new_peak)

    return peaks

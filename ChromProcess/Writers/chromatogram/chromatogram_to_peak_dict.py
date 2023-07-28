from ChromProcess.Classes import Chromatogram


def chromatogram_to_peak_dict(chromatogram: Chromatogram) -> dict[str, list[float]]:
    """
    Create a spreadsheet-like dictionary of the peaks in a chromatogram.

    Parameters
    ----------
    chromatogram: Chromatogram

    Returns
    -------
    peaks: dict[str, list[float]]
    """

    peaks: dict[str, list[float]] = {
        "retention_time": [],
        "integral": [],
        "height": [],
        "start": [],
        "end": [],
        "is_retention_time": [],
        "is_integral": [],
        "deconvoluted": [],
    }

    complete_peaks = list(chromatogram.peaks.values()) + list(
        chromatogram.deconvoluted_peaks.values()
    )

    for peak in complete_peaks:
        integral = peak.integral
        retention_time = peak.retention_time
        height = peak.height
        start = chromatogram.time[peak.indices[0]]
        end = chromatogram.time[peak.indices[-1]]

        peaks["retention_time"].append(retention_time)
        peaks["integral"].append(integral)
        peaks["height"].append(height)
        peaks["start"].append(start)
        peaks["end"].append(end)

        if peak.deconvolution_params:
            peaks["deconvoluted"].append(True)
        else:
            peaks["deconvoluted"].append(False)

    is_rt = chromatogram.internal_standard.retention_time
    peaks["is_retention_time"] = [is_rt for p in complete_peaks]

    is_integral = chromatogram.internal_standard.integral
    peaks["is_integral"] = [is_integral for p in complete_peaks]

    return peaks

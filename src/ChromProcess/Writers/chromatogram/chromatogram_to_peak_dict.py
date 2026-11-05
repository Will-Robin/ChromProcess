from typing import Union
from ChromProcess.Classes import Peak
from ChromProcess.Classes import Chromatogram


def chromatogram_to_peak_dict(
    chromatogram: Chromatogram,
) -> dict[str, list[Union[float, str]]]:
    """
    Create a spreadsheet-like dictionary of the peaks in a chromatogram.

    Parameters
    ----------
    chromatogram: Chromatogram

    Returns
    -------
    peaks: dict[str, list[Union[float, str]]]
    """

    complete_peaks = (
        chromatogram.peaks
        + chromatogram.deconvoluted_peaks
    )

    peaks: dict[str, list[Union[float, str]]] = {
        "filename": [],
        "retention_time": [],
        "integral": [],
        "height": [],
        "start": [],
        "end": [],
        "is_retention_time": [],
        "is_integral": [],
        "deconvoluted": [],
    }

    for peak in complete_peaks:
        integral = peak.integral
        retention_time = peak.retention_time
        height = peak.height
        start = peak.start
        end = peak.end

        peaks["retention_time"].append(retention_time)
        peaks["integral"].append(integral)
        peaks["height"].append(height)
        peaks["start"].append(start)
        peaks["end"].append(end)

        if peak.deconvolution_params:
            peaks["deconvoluted"].append(True)
        else:
            peaks["deconvoluted"].append(False)

    if chromatogram.internal_standard is not None:
        is_rt = chromatogram.internal_standard.retention_time
        is_integral = chromatogram.internal_standard.integral
    else:
        is_rt = None
        is_integral = None

    peaks["is_retention_time"] = [is_rt for p in complete_peaks]
    peaks["is_integral"] = [is_integral for p in complete_peaks]

    peaks["filename"] = [chromatogram.filename for _ in complete_peaks]

    return peaks

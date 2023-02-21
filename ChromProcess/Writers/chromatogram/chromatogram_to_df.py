import pandas as pd
from ChromProcess.Classes import Chromatogram


def chromatogram_to_df(chromatogram: Chromatogram) -> pd.DataFrame:
    """
    Create a Pandas data frame object from the peaks stored in a chromatogram.

    Parameters
    ----------
    chromatogram: Chromatogram

    Returns
    -------
    df: pd.DataFrame
    """

    peaks = {
        "retention_time": [],
        "integral": [],
        "start": [],
        "end": [],
        "is_retention_time": [],
        "is_integral": [],
    }
    for peak in chromatogram.peaks.values():
        integral = peak.integral
        rt = peak.retention_time
        start = chromatogram.time[peak.indices[0]]
        end = chromatogram.time[peak.indices[-1]]
        peaks["retention_time"].append(rt)
        peaks["integral"].append(integral)
        peaks["start"].append(start)
        peaks["end"].append(end)

    peaks["is_retention_time"] = [
        chromatogram.internal_standard.retention_time for p in chromatogram.peaks
    ]
    peaks["is_integral"] = [
        chromatogram.internal_standard.integral for p in chromatogram.peaks
    ]

    df = pd.DataFrame(peaks)

    return df

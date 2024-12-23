import pandas as pd
from typing import Union
from ChromProcess.Classes import Chromatogram
from .chromatogram_to_peak_dict import chromatogram_to_peak_dict


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

    peaks: dict[str, list[Union[float, str]]] = chromatogram_to_peak_dict(chromatogram)

    df = pd.DataFrame(peaks)

    return df

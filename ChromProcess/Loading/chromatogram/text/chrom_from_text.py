import numpy as np
from ChromProcess.Classes import Chromatogram


def chrom_from_text(
    x_values: list[float],
    y_values: list[float],
    x_unit: str,
    y_unit: str,
    filename: str,
) -> Chromatogram:
    """
    Create a chromatogram object and insert time and signal information into
    it.

    Parameters
    ----------
    x_values: list
    y_values: list
    x_unit: str
    y_unit: str
    filename: str

    Returns
    -------
    chrom: Chromatogram
        Chromatogram derived from information in the file.

    """

    assert isinstance(x_values, list), "x_values arg should be a list"
    assert isinstance(y_values, list), "y_values arg should be list"
    assert isinstance(x_unit, str), "x_unit arg should be str"
    assert isinstance(y_unit, str), "y_unit arg should be str"

    chrom = Chromatogram()
    chrom.filename = filename

    chrom.x_unit = x_unit
    chrom.y_unit = y_unit

    x_values = [float(x) for x in x_values]
    y_values = [float(y) for y in y_values]

    chrom.time = np.array(x_values)
    chrom.signal = np.array(y_values)

    return chrom

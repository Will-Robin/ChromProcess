from pathlib import Path
from ChromProcess import Classes

from ChromProcess.Loading.chromatogram.cdf.instruments import JEOL
from ChromProcess.Loading.chromatogram.cdf.cdf_loading import load_from_cdf


def chrom_from_cdf(filename, instrument=JEOL, load_ms=False):
    """
    Load a chromatogram from .cdf file.

    Parameters
    ----------
    filename: str or pathlib Path

    Returns
    -------
    chrom: ChromProcess Chromatogram
    """

    if isinstance(filename, str):
        fname = Path(filename)
    else:
        fname = filename

    data_keys = [instrument.TIME_KEY, instrument.TIC_KEY]

    if load_ms is True:
        ms_extensions = [
            instrument.MS_INTENSITY_KEY,
            instrument.MZ_KEY,
            instrument.MS_INTENSITY_KEY,
            instrument.SCAN_INDEX_KEY,
            instrument.POINT_COUNT_KEY,
        ]

        data_keys.extend(ms_extensions)

    data_container = load_from_cdf(fname, data_keys)

    chrom = Classes.Chromatogram()

    chrom.time = data_container[instrument.TIME_KEY] / instrument.TIME_CONVERSION
    chrom.signal = data_container[instrument.TIC_KEY]

    chrom.filename = fname.name

    chrom.x_unit = instrument.X_UNIT
    chrom.y_unit = instrument.Y_UNIT

    if load_ms is True:
        chrom.mz_values = data_container[instrument.MZ_KEY]
        chrom.mz_intensity = data_container[instrument.MS_INTENSITY_KEY]
        chrom.scan_indices = data_container[instrument.SCAN_INDEX_KEY]
        chrom.point_counts = data_container[instrument.POINT_COUNT_KEY]

    return chrom

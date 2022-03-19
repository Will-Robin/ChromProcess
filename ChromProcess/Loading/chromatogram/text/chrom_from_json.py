import json
import numpy as np
from pathlib import Path

from ChromProcess import Classes


def chrom_from_json(filename):
    """
    A specific parser for .txt files exported from Shimadzu LabSolutions
    software.

    Extracts data from the chromatogram file into a dictionary using string
    manipulation and regex parsing (not all information in the file is
    scraped).

    Parameters
    ----------
    filename: str or pathlib Path
        Name of chromatogram file (including the path to the file).

    Returns
    -------
    data_container: dict
        Dictionary containing data scraper from the file.
    """

    if isinstance(filename, str):
        fname = Path(filename)
    else:
        fname = filename

    assert isinstance(fname, Path), "filename should be string or pathlib Path"

    with open(fname, "r") as file:
        text = file.read()

    chrom_dict = json.loads(text)

    chrom = Classes.Chromatogram()

    chrom.time = np.array(chrom_dict["time"])
    chrom.signal = np.array(chrom_dict["signal"])

    for peak in chrom_dict["peaks"]:

        retention_time = peak["retention_time"]
        start = peak["start"]
        end = peak["end"]

        inds = np.where((chrom.time >= start) & (chrom.time <= end))[0]
        chrom.peaks[retention_time] = Classes.Peak(
            retention_time, start, end, indices=inds, integral=peak["integral"]
        )

    return chrom

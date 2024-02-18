import json
import numpy as np
from pathlib import Path

from ChromProcess.Classes import Peak
from ChromProcess.Classes import Chromatogram


def chrom_from_json(filename: str) -> Chromatogram:
    """
    Load a chromatogram from a .json file:

    Expected structure:

    ```json
    {
        "time": [],
        "signal": [],
        "peaks": [
                    {
                        "retention_time": float,
                        "start": float,
                        "end": float,
                        "integral": float,
                    },
        ]

    }
    ```

    Parameters
    ----------
    filename: str or pathlib Path
        Name of chromatogram file (including the path to the file).

    Returns
    -------
    chrom: Chromatogram
        Chromatogram derived from information in the file.
    """

    if isinstance(filename, str):
        fname = Path(filename)
    else:
        fname = filename

    assert isinstance(fname, Path), "filename should be string or pathlib Path"

    with open(fname, "r") as file:
        text = file.read()

    chrom_dict = json.loads(text)

    chrom = Chromatogram()

    chrom.time = np.array(chrom_dict["time"])
    chrom.signal = np.array(chrom_dict["signal"])

    for peak in chrom_dict["peaks"]:
        retention_time = peak["retention_time"]
        start = peak["start"]
        end = peak["end"]

        inds = np.where((chrom.time >= start) & (chrom.time <= end))[0].tolist()
        chrom.peaks[retention_time] = Peak(
            retention_time, start, end, indices=inds, integral=peak["integral"]
        )

    return chrom

import numpy as np
from pathlib import Path
from ChromProcess.Classes import Chromatogram


def chrom_from_csv(file_path: str | Path) -> Chromatogram:
    """
    Load a chromatogram from a .csv file.

    Expected format:

    ```csv
    time, signal\n
    float,float\n
    float,float\n
    float,float\n
    ```

    Parameters
    ----------
    file_path: str or pathlib Path
        Name of chromatogram file (including the path to the file).

    Returns
    -------
    chrom: Chromatogram
        Chromatogram derived from information in the file.

    """

    chromatogram = Chromatogram()

    if isinstance(file_path, Path):
        chromatogram.filename = file_path.name
    else:
        chromatogram.filename = Path(file_path).name

    x_vals = []
    y_vals = []
    with open(file_path, "r") as file:
        for c, line in enumerate(file):
            if c == 0:
                header = line.strip("\n").split(",")
                chromatogram.x_unit = header[0]
                chromatogram.y_unit = header[1]
            else:
                time, signal = map(float, line.split(","))
                x_vals.append(time)
                y_vals.append(signal)

    chromatogram.time = np.array(x_vals)
    chromatogram.signal = np.array(y_vals)

    return chromatogram

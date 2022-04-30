from pathlib import Path

from ChromProcess.Loading.parsers import parsers
from ChromProcess.Loading.chromatogram.text import chrom_from_text


def chrom_from_csv(filename):
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
    filename: str or pathlib Path
        Name of chromatogram file (including the path to the file).

    Returns
    -------
    chrom: ChromProcess.Classes.Chromatogram
        Chromatogram derived from information in the file.
 
    """

    if isinstance(filename, str):
        fname = Path(filename)
    else:
        fname = filename

    assert isinstance(fname, Path), "filename should be string or pathlib Path"

    with open(fname, "r") as file:
        text = file.read()

    data = parsers.parse_text_columns(text, "\n", ",")

    x_values = [float(x) for x in data[0][1:]]
    y_values = [float(y) for y in data[1][1:]]

    x_name = data[0][0]
    y_name = data[1][0]

    chrom = chrom_from_text.chrom_from_text(
        x_values, y_values, x_name, y_name, fname.name
    )

    return chrom

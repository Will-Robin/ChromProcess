def write_chromatogram_csv_text(chromatogram):
    """
    Write the text for a chromatogram in a column-wise csv format.

    Parameters
    ----------
    chromatogram: Classes.Chromatogram
        Chromatogram to be converted.

    Returns
    -------
    csv_string: str
        Chromatogram time and signal in column-wise csv format.
    """

    time = chromatogram.time
    signal = chromatogram.signal
    x_unit = chromatogram.x_unit
    y_unit = chromatogram.y_unit

    csv_string = ""

    csv_string += f"{x_unit},{y_unit}\n"
    for x in range(0, len(time)):
        csv_string += f"{time[x]},{signal[x]}\n"

    return csv_string


def chromatogram_to_csv(chromatogram, filename=""):
    """
    Write chromatogram to a .csv file.

    Parameters
    ----------
    chromatogram: Classes.Chromatogram
        Chromatogram to be converted.
    filename: str
        Name for the file.

    Returns
    -------
    None
    """

    if filename == "":
        filename = chromatogram.filename.split(".")[0] + ".csv"

    csv_text = write_chromatogram_csv_text(chromatogram)

    with open(filename, "w") as f:
        f.write(csv_text)

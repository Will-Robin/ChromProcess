from ChromProcess.Writers.peak.peak_to_entry_text import peak_to_entry_text


def write_peak_collection_text(chromatogram, header_text=""):
    """
    Create the text for a peak collection based on the Peak objects in the
    chromatogram.

    Parameters
    ----------
    chromatogram: Classes.Chromatogram
    header_text: string

    Returns
    -------
    peak_collection_string: str
    """

    peak_collection_string = ""

    if header_text != "":
        peak_collection_string += header_text

    peak_collection_string += "IS_retention_time/ min,"
    peak_collection_string += "IS_integral,IS_peak start/ min,"
    peak_collection_string += "IS_peak end/ min\n"

    i_s = chromatogram.internal_standard
    IS_entry = peak_to_entry_text(i_s, chromatogram)

    peak_collection_string += IS_entry

    peak_collection_string += "Retention_time/ min,"
    peak_collection_string += "integral,peak start/ min,"
    peak_collection_string += "peak end/ min\n"

    for p in chromatogram.peaks:
        peak = chromatogram.peaks[p]
        peak_collection_string += peak_to_entry_text(peak, chromatogram)

    return peak_collection_string


def chromatogram_to_peak_collection(
    chromatogram, filename="peak_collection.csv", header_text=""
):
    """
    For writing peak integrals from a chromatogram to a .csv file.

    Parameters
    ----------
    filename: str or pathlib Path
        Name for the file
    header_text: str
       Text to place at the top of the peak table.

    Returns
    -------
    None
    """

    output_text = write_peak_collection_text(chromatogram, header_text=header_text)
    with open(filename, "w") as f:
        f.write(output_text)

from ChromProcess.Classes import Peak
from ChromProcess.Classes import Chromatogram

def peak_to_entry_text(peak: Peak, chromatogram: Chromatogram) -> str:
    """
    For writing the peak's row in a peak table.

    Parameters
    ----------
    peak: Peak
    chromatogram: Chromatogram

    Returns
    -------
    entry: str
        Row for the peak table output
    """

    # Default values
    st_ind = 0
    end_ind = 0
    peak_start = 0.0
    peak_end = 0.0
    peak_height = 0.0

    if len(peak.indices) > 0:
        st_ind = peak.indices[0]
        end_ind = peak.indices[-1]

        peak_start = chromatogram.time[st_ind]
        peak_end = chromatogram.time[end_ind]
        peak_height = peak.height

    rtn_time = peak.retention_time
    integral = peak.integral

    entry_list = [rtn_time, integral, peak_start, peak_end, peak_height, "\n"]

    entry = ",".join([str(x) for x in entry_list])

    return entry

def write_peak_collection_text(
    chromatogram: Chromatogram, header_text: str = ""
) -> str:
    """
    Create the text for a peak collection based on the Peak objects in the
    chromatogram.

    Parameters
    ----------
    chromatogram: Chromatogram
    header_text: str

    Returns
    -------
    peak_collection_string: str
    """

    peak_collection_string = ""

    if header_text != "":
        peak_collection_string += header_text

    peak_collection_string += "IS_retention_time/ min,"
    peak_collection_string += "IS_integral,IS_peak start/ min,"
    peak_collection_string += "IS_peak end/ min,"
    peak_collection_string += f"peak height/ {chromatogram.y_unit}"
    peak_collection_string += "\n"

    i_s = chromatogram.internal_standard
    IS_entry = peak_to_entry_text(i_s, chromatogram)

    peak_collection_string += IS_entry

    peak_collection_string += "Retention_time/ min,"
    peak_collection_string += "integral,peak start/ min,"
    peak_collection_string += "peak end/ min,"
    peak_collection_string += f"peak height/ {chromatogram.y_unit}"
    peak_collection_string += "\n"

    for p in chromatogram.peaks:
        peak = chromatogram.peaks[p]
        peak_collection_string += peak_to_entry_text(peak, chromatogram)

    return peak_collection_string


def chromatogram_to_peak_table(
    chromatogram: Chromatogram,
    filename: str = "peak_collection.csv",
    header_text: str = "",
) -> None:
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

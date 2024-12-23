from ChromProcess.Classes import Chromatogram
from .chromatogram_to_peak_dict import chromatogram_to_peak_dict


def chromatogram_to_peak_table_text(
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
    output_text: str
    """

    peak_dict = chromatogram_to_peak_dict(chromatogram)

    headings = list(peak_dict.keys())
    output_text = ",".join(headings)
    output_text += "\n"

    for c, peak in enumerate(peak_dict["retention_time"]):
        entry = ",".join([str(peak_dict[h][c]) for h in headings])
        output_text += f"{entry}\n"

    return output_text


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

    output_text = chromatogram_to_peak_table_text(chromatogram, header_text=header_text)

    with open(filename, "w") as f:
        f.write(output_text)

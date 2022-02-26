import re
from pathlib import Path

from ChromProcess import Classes

from ChromProcess.Loading.chromatogram.text import chrom_from_text
from ChromProcess.Loading.parsers import parsers


def chrom_from_labsolutions_ascii(filename, data_key="Detector A-Ch1"):
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

    blocks = text.split("\n\n")

    item_regex = r"(?:[A-Z][a-z]*\()(.*)(?:\)\])"
    data_regex = r"(?:Intensity\n)([\s\S]*)"
    x_units_regex = r"Intensity\sUnits\s(.+)"
    y_units_regex = r"R.Time\s\((.+)\)"

    data_container = {}
    for b in blocks:
        name_segment = re.findall(item_regex, b)
        data_segment = re.findall(data_regex, b)
        x_units_segment = re.findall(x_units_regex, b)
        y_units_segment = re.findall(y_units_regex, b)

        # the following conditions are met if the block contains a
        # data trace.
        if len(data_segment) > 0 and len(name_segment) > 0:
            name = name_segment[0]
            data = parsers.parse_text_columns(data_segment[0], "\n", "\t")
            data_container[name] = {
                "data": data,
                "x_unit": x_units_segment[0],
                "y_unit": y_units_segment[0],
            }

    if "data_key" in data_container:
        chrom = chrom_from_text.chrom_from_text(
            data_container[data_key]["data"][0],
            data_container[data_key]["data"][1],
            data_container[data_key]["x_unit"],
            data_container[data_key]["y_unit"],
            fname.name,
        )

        return chrom

    else:
        print("Scraping data from file failed.")

    return Classes.Chromatogram()

import json


def write_chromatogram_json_text(chromatogram):
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

    chrom_dict = {}

    chrom_dict["time"] = chromatogram.time.tolist()
    chrom_dict["signal"] = chromatogram.signal.tolist()
    chrom_dict["x_unit"] = chromatogram.x_unit
    chrom_dict["y_unit"] = chromatogram.y_unit

    chrom_dict["peaks"] = []
    for p_k in chromatogram.peaks:
        peak = chromatogram.peaks[p_k]
        chrom_dict["peaks"].append(
            {
                "retention_time": peak.retention_time,
                "start": peak.start,
                "end": peak.end,
                "integral": peak.integral,
            }
        )

    json_string = json.dumps(chrom_dict)

    return json_string


def chromatogram_to_json(chromatogram, filename=""):
    """
    Write chromatogram to a .json file.

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

    json_text = write_chromatogram_json_text(chromatogram)

    with open(filename, "w") as f:
        f.write(json_text)

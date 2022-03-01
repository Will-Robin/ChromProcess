from pathlib import Path

def write_peak_collection_text(peak_collection):
    """
    Write a PeakCollection to a csv string.

    Parameters
    ----------
    peak_collection: Classes.PeakCollection
    directory: str or pathlib Path

    Returns
    -------
    peak_collection_text: string
    """

    unit = peak_collection.series_unit
    value = peak_collection.series_value

    IS_rt = peak_collection.internal_standard.retention_time
    IS_integ = peak_collection.internal_standard.integral
    strt = peak_collection.internal_standard.start
    end = peak_collection.internal_standard.end

    peak_collection_text = f"{unit},{value}\n"
    peak_collection_text += "IS_retention_time/ min,"
    peak_collection_text += "IS_integral,"
    peak_collection_text += "IS_peak start/ min,"
    peak_collection_text += "IS_peak end/ min\n"
    peak_collection_text += f"{IS_rt},{IS_integ},{strt},{end}\n"

    peak_collection_text += "Retention_time/ min,"
    peak_collection_text += "integral,"
    peak_collection_text += "peak start/ min,"
    peak_collection_text += "peak end/ min\n"

    for p in peak_collection.peaks:
        rt = p.retention_time
        integ = p.integral
        start = p.start
        end = p.end
        peak_collection_text += f"{rt},{integ},{start},{end}\n"

    return peak_collection_text

def peak_collection_to_csv(peak_collection, directory=""):
    """
    Write a PeakCollection to a csv file.

    Parameters
    ----------
    peak_collection: ChromProcess.Classes.PeakCollection
    directory: str or pathlib Path

    Returns
    -------
    None
    """

    fname = Path(peak_collection.filename)
    if isinstance(directory, str):
        if directory == "":
            fname = Path(f"{peak_collection.filename}")
    elif isinstance(directory, Path):
        fname = directory / f"{peak_collection.filename}"

    text = write_peak_collection_text(peak_collection)

    with open(fname, "w") as f:
        f.write(text)


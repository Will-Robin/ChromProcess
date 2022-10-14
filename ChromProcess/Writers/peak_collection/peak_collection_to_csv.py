from pathlib import Path

IS_ENTRY_HEADER = ",".join(
    [
        "IS_retention_time/ min,",
        "IS_integral,",
        "IS_peak start/ min,",
        "IS_peak end/ min\n",
    ]
)

PEAK_ENTRY_HEADER = ",".join(
    ["Retention_time/ min,", "integral,", "peak start/ min,", "peak end/ min\n"]
)


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
    is_line = ",".join([IS_rt, IS_integ, strt, end, "\n"])

    peak_collection_text = f"{unit},{value}\n"
    peak_collection_text += IS_ENTRY_HEADER
    peak_collection_text += is_line
    peak_collection_text += PEAK_ENTRY_HEADER

    for p in peak_collection.peaks:

        peak_line = ",".join([p.retention_time, p.integral, p.start, p.end, "\n"])

        peak_collection_text += peak_line

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

from pathlib import Path


def peak_collection_to_csv(peak_collection, directory=""):
    """
    Write a PeakCollection to a csv file.

    Parameters
    ----------
    peak_collection: Classes.PeakCollection
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

    IS_header = ""
    IS_header += "IS_retention_time/ min,"
    IS_header += "IS_integral,"
    IS_header += "IS_peak start/ min,"
    IS_header += "IS_peak end/ min\n"

    pk_header = ""
    pk_header += "Retention_time/ min,"
    pk_header += "integral,"
    pk_header += "peak start/ min,"
    pk_header += "peak end/ min\n"

    unit = peak_collection.series_unit
    value = peak_collection.series_value
    with open(fname, "w") as f:
        f.write(f"{unit},{value}\n")
        f.write(IS_header)

        IS_rt = peak_collection.internal_standard.retention_time
        IS_integ = peak_collection.internal_standard.integral
        strt = peak_collection.internal_standard.start
        end = peak_collection.internal_standard.end

        f.write(f"{IS_rt},{IS_integ},{strt},{end}\n")
        f.write(pk_header)
        for p in peak_collection.peaks:
            rt = p.retention_time
            integ = p.integral
            start = p.start
            end = p.end
            f.write(f"{rt},{integ},{start},{end}\n")

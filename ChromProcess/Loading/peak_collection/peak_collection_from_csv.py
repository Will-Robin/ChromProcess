import sys
from pathlib import Path
from ChromProcess import Classes


def peak_collection_from_csv(filename, round_digits=3):
    """
    Create a PeakCollection object from a formatted file.

    Parameters
    ----------
    filename: str or pathlib Path
    round_digits: int

    Returns
    -------
    peak_collection: Classes.PeakCollection
    """

    peak_collection = Classes.PeakCollection()

    fname = Path("")
    if isinstance(filename, str):
        fname = Path(filename)
    elif isinstance(filename, Path):
        fname = filename
    else:
        sys.exit(
            """PeakCollection requires file kwarg to be string or
        pathlib Path."""
        )

    peak_collection.filename = fname.name

    read_line = lambda line: [float(x) for x in line.strip("\n").split(",") if x != ""]
    peaks = []

    IS_line_num = -1

    IS = Classes.Peak(0.0, 0.0, 0.0, integral=1.0)

    value = 0.0
    variable = ""
    with open(fname, "r") as f:

        parent_filename = peak_collection.filename.split(".")[0]

        for c, line in enumerate(f):
            if "None" in line:
                pass

            elif c == 0:
                read = [x for x in line.strip("\n").split(",") if x != ""]
                variable = read[0]
                value = float(read[1])

            elif "IS_" in line:
                IS_line_num = c + 1

            elif c == IS_line_num:
                if "None" in line:
                    pass

                else:
                    read = read_line(line)

                    IS = Classes.Peak(
                        round(read[0], round_digits),
                        round(read[2], round_digits),
                        round(read[3], round_digits),
                        integral=read[1],
                        parent=parent_filename,
                    )

            elif c < 4:
                pass

            else:
                rd = read_line(line)
                peaks.append(
                    Classes.Peak(
                        round(rd[0], round_digits),
                        round(rd[2], round_digits),
                        round(rd[3], round_digits),
                        integral=rd[1],
                        parent=parent_filename,
                    )
                )

    peak_collection.series_value = value
    peak_collection.series_unit = variable
    peak_collection.internal_standard = IS
    peak_collection.peaks = peaks
    peak_collection.mass_spectra = []
    peak_collection.initial_IS_pos = IS.retention_time

    return peak_collection

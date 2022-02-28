from pathlib import Path

from ChromProcess.Utils.utils import utils
from ChromProcess.Writers.general import write_header


def peak_collection_series_to_data_report(
    peak_collection_series, filename, information
):
    """
    Write a peak collection series as a formatted data report file.

    Parameters
    ----------
    peak_collection_series: Classes.PeakCollectionSeries
    filename: name for file including path
    information: ChromProcess Analysis_Information object

    Returns
    -------
    None
    """

    if isinstance(filename, str):
        filename = filename
    elif isinstance(filename, Path):
        filename = str(filename)

    analysis_type = information.analysis_type

    conc_fname = f"{filename}_{analysis_type}_concentration_report.csv"
    integral_fname = f"{filename}_{analysis_type}_integral_report.csv"

    # create output dictionaries
    conc_dict, err_dict, integral_dict = peak_collection_series.series_traces_as_dict()
    
    # create spreadsheet-like output
    conc_header, conc_grid = utils.peak_dict_to_spreadsheet(
        conc_dict,
        peak_collection_series.series_values,
        peak_collection_series.series_unit,
    )

    peak_integral_header, integ_grid = utils.peak_dict_to_spreadsheet(
        integral_dict,
        peak_collection_series.series_values,
        peak_collection_series.series_unit,
    )

    peak_err_header, err_grid = utils.peak_dict_to_spreadsheet(
        err_dict,
        peak_collection_series.series_values,
        peak_collection_series.series_unit,
    )

    header_text = write_header.write_conditions_header(
        peak_collection_series.name, peak_collection_series.conditions, information
    )

    # Write concentration report to file
    with open(conc_fname, "w") as outfile:

        outfile.write(header_text)

        outfile.write("start_data\n")

        [outfile.write("{},".format(x)) for x in conc_header]

        outfile.write("\n")

        for x in range(0, len(conc_grid)):
            for y in range(0, len(conc_grid[x])):
                val = conc_grid[x][y]
                outfile.write(f"{val},")
            outfile.write("\n")

        outfile.write("end_data\n")

        outfile.write("start_errors\n")

        [outfile.write("{},".format(x)) for x in peak_err_header]

        outfile.write("\n")

        for x in range(0, len(err_grid)):
            for y in range(0, len(err_grid[x])):
                val = err_grid[x][y]
                outfile.write(f"{val},")
            outfile.write("\n")

        outfile.write("end_errors\n")

    # Write integral report to file
    with open(integral_fname, "w") as outfile:

        outfile.write(header_text)

        outfile.write("start_data\n")

        [outfile.write(f"{x},") for x in peak_integral_header]

        outfile.write("\n")
        for x in range(0, len(integ_grid)):
            for y in range(0, len(integ_grid[x])):
                val = integ_grid[x][y]
                outfile.write(f"{val},")
            outfile.write("\n")

        outfile.write("end_data\n")

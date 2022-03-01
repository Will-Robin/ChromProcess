from pathlib import Path

from ChromProcess.Utils.utils import utils
from ChromProcess.Writers.general import write_header

def peak_collection_series_conc_report_text(peak_collection_series, information):
    """
    Write a peak collection series as formatted data report text.

    Parameters
    ----------
    peak_collection_series: Classes.PeakCollectionSeries
    information: ChromProcess Analysis_Information object

    Returns
    -------
    conc_report_text: str
    """

    # create output dictionaries
    conc_dict, err_dict, _ = peak_collection_series.series_traces_as_dict()

    # create spreadsheet-like output
    conc_header, conc_grid = utils.peak_dict_to_spreadsheet(
        conc_dict,
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

    conc_report_text = header_text

    conc_report_text += "start_data\n"

    for x in conc_header:
        conc_report_text += f"{x},"

    conc_report_text += "\n"

    for x in range(0, len(conc_grid)):
        for y in range(0, len(conc_grid[x])):
            val = conc_grid[x][y]
            conc_report_text += f"{val},"
        conc_report_text += "\n"

    conc_report_text += "end_data\n"

    conc_report_text += "start_errors\n"

    for x in peak_err_header:
        conc_report_text += f"{x},"

    conc_report_text += "\n"

    for x in range(0, len(err_grid)):
        for y in range(0, len(err_grid[x])):
            val = err_grid[x][y]
            conc_report_text += f"{val},"
        conc_report_text += "\n"

    conc_report_text += "end_errors\n"

    return conc_report_text

def peak_collection_series_integral_report_text(peak_collection_series,
        information):
    """
    Write a peak collection series as formatted data report text.

    Parameters
    ----------
    peak_collection_series: Classes.PeakCollectionSeries
    information: ChromProcess Analysis_Information object

    Returns
    -------
    integral_report_text: str
    """

    _, _, integral_dict = peak_collection_series.series_traces_as_dict()

    header_text = write_header.write_conditions_header(
        peak_collection_series.name, peak_collection_series.conditions, information
    )

    peak_integral_header, integ_grid = utils.peak_dict_to_spreadsheet(
        integral_dict,
        peak_collection_series.series_values,
        peak_collection_series.series_unit,
    )

    integral_report_text = header_text

    integral_report_text += "start_data\n"

    for x in peak_integral_header:
        integral_report_text += f"{x},"
    integral_report_text += "\n"

    for x in range(0, len(integ_grid)):
        for y in range(0, len(integ_grid[x])):
            val = integ_grid[x][y]
            integral_report_text += f"{val},"
        integral_report_text += "\n"

    integral_report_text += "end_data\n"

    return integral_report_text

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
    conc_rep_text = peak_collection_series_conc_report_text(peak_collection_series, information)
    with open(conc_fname, "w") as outfile:
        outfile.write(conc_rep_text)

    # Write integral report to file
    i_data_rep_text = peak_collection_series_integral_report_text(peak_collection_series, information)
    with open(integral_fname, "w") as outfile:
        outfile.write(i_data_rep_text)


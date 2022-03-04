from ChromProcess.Utils.utils import utils
from ChromProcess.Writers.general import write_header


def write_data_report_text(data_report):
    """
    Write a formatted text version of the data report.

    Parameters
    ----------
    data_report: Classes.DataReport
        DataReport to convert to text.

    Returns
    -------
    data_report_text: str
        Text output for the data report.
    """

    header_text = write_header.write_conditions_header(
        data_report.name, data_report.conditions, data_report.analysis
    )

    data_header, data_grid = utils.peak_dict_to_spreadsheet(
        data_report.data, data_report.series_values, data_report.series_unit
    )

    peak_err_header, err_grid = utils.peak_dict_to_spreadsheet(
        data_report.errors, data_report.series_values, data_report.series_unit
    )

    data_report_text = ""

    data_report_text += header_text

    data_report_text += "start_data\n"

    for x in data_header:
        data_report_text += f"{x},"

    data_report_text += "\n"

    for x in range(0, len(data_grid)):
        for y in range(0, len(data_grid[x])):
            val = data_grid[x][y]
            data_report_text += f"{val},"
        data_report_text += "\n"

    data_report_text += "end_data\n"

    data_report_text += "start_errors\n"

    for x in peak_err_header:
        data_report_text += f"{x},"

    data_report_text += "\n"

    for x in range(0, len(err_grid)):
        for y in range(0, len(err_grid[x])):
            val = err_grid[x][y]
            data_report_text += f"{val},"
        data_report_text += "\n"

    data_report_text += "end_errors\n"

    return data_report_text


def data_report_to_csv(data_report, filename=""):
    """
    Write a data report to a csv file.

    Parameters
    ----------
    data_report: Classes.DataReport
        DataReport to write to file.
    filename: str
        name for file
    path: pathlib Path object
        Path to folder for file storage.

    Returns
    -------
    None
    """

    if filename == "":
        filename = data_report.filename

    data_report_as_string = write_data_report_text(data_report)

    with open(filename, "w") as file:
        file.write(data_report_as_string)

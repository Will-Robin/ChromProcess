import numpy as np
from pathlib import Path

from ChromProcess import Classes
from ChromProcess.Loading import parsers

# TODO: clean up this code (and use regex instead?)


def data_report_from_csv(file):
    """
    Load a data report object from a formatted .csv file.

    Parameters
    ----------
    file: pathlib Path or str
        Path to file.

    Returns
    -------
    data_report: Classes.DataReport
    """

    data_report = Classes.DataReport()

    spl_lin = lambda x: [e for e in x.strip("\n").split(",") if e != ""]

    if type(file) == str:
        file = Path(file)

    data_report.filename = file.name

    with open(file, "r", encoding="latin-1") as f:
        for line in f:
            ins = spl_lin(line)
            if "Dataset" in line:
                data_report.experiment_code = ins[1]

    # parse conditions section
    condset = parsers.import_file_section(file, "start_conditions", "end_conditions")

    for c in condset:
        entry = [float(x) for x in c[1:]]
        if len(entry) == 1:
            data_report.conditions[c[0]] = entry[0]
        else:
            data_report.conditions[c[0]] = np.array(entry)

    # parse data section
    dataset = parsers.import_file_section(file, "start_data", "end_data")

    transposed_datalines = [list(i) for i in zip(*dataset)]
    d_out = {}
    for s in transposed_datalines:
        d_out[s[0]] = np.array([0 if x == "nan" else float(x) for x in s[1:]])

    data_report.series_unit = dataset[0][0]
    data_report.series_values = d_out[data_report.series_unit]

    del d_out[data_report.series_unit]

    data_report.data = d_out

    # parse errors section
    errors = parsers.import_file_section(file, "start_errors", "end_errors")

    if len(errors) == 0:
        data_report.errors = {
            d: np.zeros(len(data_report.series_values)) for d in data_report.data
        }
    else:
        transposed_error_lines = [list(i) for i in zip(*errors)]
        errors_out = {}
        for s in transposed_error_lines:
            errors_out[s[0]] = np.array([0 if x == "nan" else float(x) for x in s[1:]])
        del errors_out[data_report.series_unit]
        data_report.errors = errors_out

    # Parse analysis details section
    analysis = parsers.import_file_section(
        file, "start_analysis_details", "end_analysis_details"
    )
    for a in analysis:
        data_report.analysis_details[a[0]] = [x for x in a[1:]]

    return data_report

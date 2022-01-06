from pathlib import Path

from ChromProcess.Utils.utils import utils
from ChromProcess.Writers.general import write_header

def data_report_to_csv(data_report, filename = '', path = None):
    '''
    Parameters
    ----------
    data_report: Classes.DataReport
        DataReport to write to file.
    filename: str
        name for file
    path: pathlib Path object
        Path to folder for file storage.
    '''

    if isinstance(path, str):
        path = Path(path)

    if filename == '':
        filename = data_report.filename

    if path is None:
        fname = filename
    else:
        fname = path/filename

    header_text = write_header.write_conditions_header(
                                            data_report.name, 
                                            data_report.conditions, 
                                            data_report.analysis
                                            )
    conc_header, conc_grid = utils.peak_dict_to_spreadsheet(
                                        data_report.data, 
                                        data_report.series_values,
                                        data_report.series_unit
                                        )
    peak_err_header, err_grid = utils.peak_dict_to_spreadsheet(
                                        data_report.errors, 
                                        data_report.series_values,
                                        data_report.series_unit
                                        )

    with open(fname, 'w') as outfile:

        outfile.write(header_text)

        outfile.write("start_data\n")

        [outfile.write("{},".format(x)) for x in conc_header]

        outfile.write("\n")

        for x in range(0,len(conc_grid)):
            for y in range(0,len(conc_grid[x])):
                val = conc_grid[x][y]
                outfile.write(f"{val},")
            outfile.write("\n")

        outfile.write("end_data\n")

        outfile.write("start_errors\n")

        [outfile.write("{},".format(x)) for x in peak_err_header]

        outfile.write("\n")

        for x in range(0,len(err_grid)):
            for y in range(0,len(err_grid[x])):
                val = err_grid[x][y]
                outfile.write(f"{val},")
            outfile.write("\n")

        outfile.write("end_errors\n")


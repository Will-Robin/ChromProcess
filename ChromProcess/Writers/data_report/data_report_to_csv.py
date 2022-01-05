import numpy as np
from pathlib import Path
from ChromProcess.Writers import write_conditions_text as w_ct

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

    header = w_ct.write_conditions_text(data_report)

    with open(fname, 'w') as outfile:
        # writing experiment conditions to file
        outfile.write(header)
        # writing data
        sorted_keys = sorted([*data_report.data], key = lambda x:x.count('C'))

        outfile.write("start_data\n")

        p_header = [data_report.series_unit]
        out = np.array([data_report.series_values])

        for s in sorted_keys:
            p_header.append(s)
            out = np.vstack((out,data_report.data[s]))

        out = out.T

        [outfile.write("{},".format(x)) for x in p_header]

        outfile.write("\n")

        for x in range(0,len(out)):
            for y in range(0,len(out[x])):
                outfile.write("{},".format(out[x,y]))
            outfile.write("\n")

        outfile.write("end_data\n")

        outfile.write('start_errors\n')

        if len(data_report.errors)> 0:
            err_out = np.array([data_report.series_values])

            for s in sorted_keys:
                trace = data_report.errors[s]
                err_out = np.vstack((err_out,trace))

            err_out = err_out.T
            [outfile.write("{},".format(x)) for x in p_header]
            outfile.write("\n")
            for x in range(0,len(err_out)):
                for y in range(0,len(err_out[x])):
                    outfile.write("{},".format(err_out[x,y]))
                outfile.write("\n")

        outfile.write('end_errors\n')


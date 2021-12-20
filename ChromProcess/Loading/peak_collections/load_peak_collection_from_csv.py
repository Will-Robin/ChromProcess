from pathlib import Path
from ChromProcess import Classes

def read_from_file(filename):
    '''
    Read in information to the object from a file.

    filename: str or pathlib Path
    '''

    peak_collection = Classes.PeakCollection()

    fname = Path('')
    if isinstance(filename, str):
        fname = Path(file)
    elif isinstance(filename, Path):
        fname = file
    else:
        sys.exit('''PeakCollection requires file kwarg to be string or
        pathlib Path.''')

    peak_collection.filename = fname.name

    read_line = lambda line: [float(x) for x in line.strip('\n').split(",") if x != '']
    peaks = []

    IS_line_num = -1
    IS = Classes.PeakCollectionElement(0.0, 1, 0.0, 0.0)

    value = 0.0
    variable = ''
    with open(file, "r") as f:
        for c,line in enumerate(f):
            if 'None' in line:
                pass
            elif c == 0:
                read = [x for x in line.strip('\n').split(',') if x != '']
                variable = read[0]
                value = float(read[1])
            elif 'IS_' in line:
                IS_line_num = c + 1

            elif c == IS_line_num:
                if 'None' in line:
                    pass
                else:
                    read = read_line(line)

                    IS = Classes.PeakCollectionElement(round(read[0],3),
                                   read[1],
                                   round(read[2],3),
                                   round(read[3],3),
                                   parent = peak_collection.filename.split('.')[0])
            elif c < 4:
                pass
            else:
                rd = read_line(line)
                peaks.append(
                Classes.PeakCollectionElement(round(rd[0],3), rd[1],
                                              round(rd[2],3),
                                              round(rd[3],3),
                                              parent = peak_collection.filename.split('.')[0]))

    peak_collection.series_value = value
    peak_collection.series_unit = variable
    peak_collection.internal_standard = IS
    peak_collection.peaks = peaks
    peak_collection.assignment = 'not specified'
    peak_collection.mass_spectra = []
    peak_collection.initial_IS_pos = IS.retention_time


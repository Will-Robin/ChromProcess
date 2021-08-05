from ChromProcess import file_import as f_i
from ChromProcess import Classes
from netCDF4 import Dataset
import numpy as np
import os

'''
Functions for importing chromatography data.
'''

def read_local_assignments(file):
    import os

    modified_bounds = {}

    if not os.path.exists(file):
        with open(file, 'w') as f:
            f.write('')
    else:
        with open(file, 'r') as f:
            for line in f:
                spl = line.strip('\n').split(',')
                spl = [x for x in spl if x != '']
                if len(spl) == 0:
                    pass
                elif 'set_IS_pos' in line:
                    modified_bounds[spl[0]] = [float(spl[1])]
                else:
                    modified_bounds[spl[0]] = [float(spl[1]), float(spl[2])]

    return modified_bounds

def get_calibration_file_allocations(file, exp_name):
    calib_file = False
    IS_pos = False
    try:
        with open(file, 'r') as f:
            for line in f:
                if exp_name in line:
                    spl = line.strip('\n').split(',')

                    calib_file = spl[1]
                    IS_pos = float(spl[2])
    except:
        pass

    if calib_file == False and IS_pos == False:
        return None, None
    else:
        return calib_file, IS_pos

def read_mass_spectra_report(file):
    from ChromProcess import simple_functions as s_f

    mass_spectra = []

    with open(file, 'r') as f:
        fltln = lambda x: [float(e) for e in x.strip('\n').split(',')[1:]
                           if e != '']

        make_ms = False
        for line in f:
            if 'Peak retention time' in line:
                read = line.strip('\n').split(',')
                rt = float(read[1])
            if 'm/z' in line:
                if s_f.isfloat(line.split(',')[1]):
                    read = fltln(line)
                    mz = np.array(read)
                else:
                    mz = np.array([0.0])

            if 'relative abundance' in line:
                if s_f.isfloat(line.split(',')[1]):
                    read = fltln(line)
                    r_a = np.array(read)
                else:
                    r_a = np.array([0.0])

                make_ms = True

            if make_ms:
                mass_spectra.append(Classes.MassSpectrum(file, mz, r_a, pos = round(rt,3)))
                make_ms = False

    return mass_spectra

def read_point_removals_file(fname):
    '''
    For reading a file containing datasets and point indices to remove in the
    data sets.
    Parameters
    ----------
    fname: str
        Path to the file

    Returns
    -------
    point_removals: dict
        dictionary of point removal indices. keys are experiment codes,
        list of ints are items.
    '''
    point_removals = {}
    with open(fname, 'r') as f:
        for c,line in enumerate(f):
            if c == 0:
                pass
            else:
                ins = line.strip('\n').split(',')
                point_removals[ins[0]] = [int(x) for x in ins[1:] if x!='']
    return point_removals

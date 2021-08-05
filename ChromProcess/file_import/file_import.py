from ChromProcess import file_import as f_i
from ChromProcess import Classes
from netCDF4 import Dataset
import numpy as np
import os

'''
Functions for importing chromatography data.
'''
def load_cdf_from_directory(directory, dir_return = False, ms = False, limit = 1e100):
    '''
    Parameters
    ----------
    directory: str or pathlib Path
        Directory containing .cdf files.

    dir_return: bool
        Whether to return to the original working directory (True) or to stay in
        directory (False).
    ms: bool
        Whether to load mass spectral information from .cdf files or not.
    limit: int
        A limit to number of files to be read.

    Returns
    -------
    chroms: list of ChromProcess Chromatogram objects
        Chromatograms obtained from the directory.
    cond_file: str or None
        File name of the conditions file. If there is no file with 'conditions'
        in the name returns None.
    '''
    cwd = os.getcwd()

    os.chdir(directory) # changes the directory to where files are stored
    filelist = os.listdir() # get list of files in directory
    filelist.sort() # sort the files

    chroms = []
    cond_file = None
    for c,f in enumerate(filelist):
        if f.endswith(".cdf") and not f.startswith("._"):
            print("loading", f)
            chroms.append(Classes.Chromatogram(f, mass_spec = ms))
        if c > limit:
            break

    for c,f in enumerate(filelist):
        if "conditions" in f and not f.startswith("._"):
            cond_file = f
            break

    if dir_return:
        os.chdir(cwd)
    else:
        print("####################################################################")
        print("Now working in {}".format(directory))
        print("####################################################################")

    return chroms, cond_file

def directoryLoadShimadzuASCII(directory, dir_return = False, limit = 1e100):

    cwd = os.getcwd()

    os.chdir(directory) # changes the directory to where files are stored
    filelist = os.listdir() # get list of files in directory
    filelist.sort() # sort the files

    chroms = []
    for f in filelist:
        if f.endswith(".txt"): # this is how the programs finds the data files: .txt is assumed to be a HPLC file, .cdf is assumed to be GCMS
            print(f)
            chroms.append(Classes.Chromatogram(f,channel_select = '215nm'))
        if "conditions" in f:
            cond_file = f

    if dir_return:
        os.chdir(cwd)
    else:
        print("####################################################################")
        print("Now working in {}".format(directory))
        print("####################################################################")

    return chroms, cond_file

def get_info_Shimadzu_HPLC(file):
    '''
    Extracts key information from the exported chromatography file (see dat_dict).

    Parameters
    ----------
    file: str
         name of chromatogram file (ASCII exported from Shimadzu LC software)
    Returns
    -------
    dat_dict: dict
        dictionary of information extracted from file

    '''

    dat_dict = {"Data set start"   : [],
                "Start Time"       : [],
                "End Time"         : [],
                "Intensity Units"  : [],
                "Time"             : [],
                "Signal"           : [],
                "Wavelength"       : []}

    with open(file, 'r') as f:
        for c,line in enumerate(f):
            if 'Start Time'in line:
                ex =  line.split('\t')
                ex = [e.replace(',','.') for e in ex]
                dat_dict['Start Time'].append(float(ex[1].strip('\n')))
            if 'R.Time (min)' in line:
                dat_dict["Data set start"].append(c+1)
            if 'End Time'in line:
                ex =  line.split('\t')
                ex = [e.replace(',','.') for e in ex]
                dat_dict['End Time'].append(float(ex[1].strip('\n')))
            if 'Intensity Units'in line:
                ex =  line.split('\t')
                ex = [e.replace(',','.') for e in ex]
                dat_dict['Intensity Units'].append(ex[1].strip('\n'))
            if 'Wavelength(nm)' in line:
                ex =  line.split('\t')
                ex = [e.replace(',','.') for e in ex]
                dat_dict['Wavelength'].append(ex[1].strip('\n'))

    return dat_dict

def get_data_Shimadzu_HPLC(file,dat_dict):
    '''
    Extracts data from the chromatogram into a dictionary.

    Parameters
    ----------
    file: str
        name of chromatogram file (ASCII exported from Shimadzu LC software)
    dat_dict: dict
        dictionary of information from chromatogram file.

    Returns
    -------
    Updated dat_dict with chromatogram data added from the file

    '''

    for b in range(0,len(dat_dict['Data set start'])):
        with open(file, 'r') as f:
            time = []
            signal = []
            for c,line in enumerate(f):
                if c > dat_dict['Data set start'][b]:
                    ex = line.strip("\n")
                    inp = ex.split("\t")
                    inp = [e.replace(',','.') for e in inp]
                    time.append(float(inp[0]))
                    signal.append(float(inp[1]))
                    if round(float(inp[0]),3) == dat_dict['End Time'][b]:
                        dat_dict['Time'].append(np.array(time))
                        dat_dict['Signal'].append(np.array(signal))
                        break
    return dat_dict

def get_info_cdf_GCMS(f):
    '''
    Opens a .cdf file and prints its contents.
    Uses the Dataset function from netCDF4 library.
    Parameters
    ----------
    f: str
        file name of the GCMS .cdf file
    '''
    f = Dataset(f, "r")
    keys = [*f.variables]
    for n in range(0,len(keys)):
        output = f.variables[keys[n]][:]
        length = output.shape
        print()
        print(keys[n], "contents: ", output,"length: ", length)
        print()
    f.close()

    return 0

def get_data_cdf_GCMS(f, key):
    '''
    Extracts data from a .cdf file using the Dataset function
    from the netCDF4 library.
    Parameters
    ----------
    f: str
        file name of the GCMS .cdf file
    key: str
        key to a set of data in the .cdf file

    Returns
    -------
    output: numpy array
        array filled with values of selected variable

    Notes
    -----
    For a list of contents of the .cdf file, use cdf_info function.

    '''
    f = Dataset(f, "r")
    f.set_auto_mask(False)
    output = f.variables[key][:]
    f.close()
    return output

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

def importCalibrationFileAllocations(filename):
    from ChromProcess import Classes
    return Classes.CalibrationAllocations(filename)

def importDataInformationFile(filename):

    '''
    Importing directories of data from


    '''
    from ChromProcess import Classes
    return Classes.DataPaths(filename)

def read_peak_table(file):
    from ChromProcess import Classes
    read_line = lambda line: [float(x) for x in line.strip('\n').split(",") if x != '']
    peaks = []
    with open(file, "r") as f:
        for c,line in enumerate(f):
            if c == 0:
                read = [x for x in line.strip('\n').split(',') if x != '']
                variable = read[0]
                value = float(read[1])
            elif c == 2:
                read = read_line(line)
                IS_pos = read[0]
                IS_integral = read[1]
                IS_bound = read[2:]
            elif c < 4:
                pass
            else:
                rd = read_line(line)
                peaks.append(Classes.PeakCollectionElement(round(rd[0],3), rd[1],
                                                           round(rd[2],3),
                                                           round(rd[3],3)))

    IS = Classes.PeakCollectionElement(round(IS_pos,3), IS_integral,
                                       round(IS_bound[0],3),
                                       round(IS_bound[1],3))

    return Classes.PeakCollection(file,variable, value, IS, peaks)

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

def load_chromatogram_csv(file):
    '''
    For loading a chromatogram from a .csv file
    Parameters
    ----------
    file: str
        Path to file
    Returns
    -------
    time: numpy array
        time axis
    signal: numpy array
        signal axis
    '''
    data = []
    with open(file, 'r') as f:

        for c,line in enumerate(f):
            if c == 0:
                pass
            else:
                data.append([float(x) for x in line.strip('\n').split(',')])

    data = np.array(data)

    time = data[:,0]
    signal = data[:,1]

    return time, signal

def read_data_from_report(file):
    dataset = []
    readstate = False
    with open(file, 'r',encoding='latin1') as f:
        for n,line in enumerate(f,0):
            if "start_data" in line:
                readstate = True
                line = next(f)
            if "end_data" in line:
                readstate = False
            if readstate:
                newline = line.strip("\n")
                dataset.append([x for x in newline.split(",") if x != ""])

    e = [list(i) for i in zip(*dataset)]
    d_out = {}
    for s in e:
        d_out[s[0]] = np.array([0 if x == 'nan' else float(x) for x in s[1:]])

    return d_out

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

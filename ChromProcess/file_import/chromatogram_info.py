from netCDF4 import Dataset

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

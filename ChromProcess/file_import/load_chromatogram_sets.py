import os
from ChromProcess import Classes

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

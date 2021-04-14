import sys
sys.path.append('/Users/williamrobinson/documents/nijmegen/packages')
from ChromProcess import Classes, os, plt, series_builder, plotting, file_output, np
from ChromProcess import curve_fit
from shutil import copyfile
from pathlib import Path

calib_folder = Path('/Users/williamrobinson/Documents/Nijmegen/Data/2020_calibrations')
os.chdir(calib_folder)

paths = {}
with open("Calibration_paths.csv", "r") as f:
    next(f)
    for line in f:
        l = line.strip("\n")
        l = l.split(",")

        paths[l[0]] = [x for x in l[1:] if x != ""]
        for c,x in enumerate(paths[l[0]]):
            paths[l[0]][c] = x.replace('\\','/')

calibration_parameters = {"Sugar":["A","B","C"]}

for p in paths:
    print(p)

    os.makedirs(p,exist_ok = True)
    datasets = []
    for c,address in enumerate(paths[p]):
        os.chdir(calib_folder/"{}/Data_reports".format(address))
        f_lst = os.listdir()
        for f in f_lst:
            if "all_integrals" in f and not "._" in f:
                #create a dataset
                print(address)
                copyfile(f, r"/Users/williamrobinson/Documents/Nijmegen/Data/2020_calibrations/{}/{}{}.csv".format(p,p,c))

    os.chdir(r"/Users/williamrobinson/Documents/Nijmegen/Data/2020_calibrations")

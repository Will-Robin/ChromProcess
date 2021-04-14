import sys
sys.path.append(r'/Users/williamrobinson/documents/nijmegen/2020_03_13_Packages')
from ChromProcess import Classes, os, plt, series_builder, plotting, file_output, np
from ChromProcess import curve_fit
from NetFit import Classes as NClass
from NetFit import file_in_out as N_io

os.chdir(r"/Users/williamrobinson/Documents/Nijmegen/2020_calibrations/")

paths = {}
with open("Calibration_paths.csv", "r") as f:
    next(f)
    for line in f:
        l = line.strip("\n")
        l = l.split(",")

        paths[l[0]] = [x for x in l[1:] if x != ""]

calibration_parameters = {"Sugar":["A","B","C"]}

for p in paths:

    datasets = []
    for address in paths[p]:
        os.chdir(r"{}/Calibration_output".format(address))

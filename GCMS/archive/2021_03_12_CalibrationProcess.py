import sys
sys.path.append(r'/Users/williamrobinson/documents/nijmegen/packages')
import os
from pathlib import Path

base_folder = Path('/Users/williamrobinson/Documents/Nijmegen/Data')
calib_folder = base_folder/'GCMS/2020_calibrations'
os.chdir(calib_folder)

paths = {}
with open(calib_folder/"Calibration_paths.csv", "r") as f:
    for c,line in enumerate(f):
        if c == 0:
            pass
        else:
            l = line.strip("\n").split(",")
            paths[l[0]] = [x for x in l[1:] if x != ""]

for p in paths:

    concentrations = []
    peak_1 = []
    peak_2 = []
    readstate = False
    for pa in paths[p]:
        print('here',pa)
        for file in os.listdir(pa):
            if 'Combined_Results' in file: continue
            if file.endswith('.png'): continue
            if file.endswith('DS_Store'): continue

            dataset = []
            with open(file, 'r', encoding='latin1') as f:
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
                insert = [0 if x == 'nan' else float(x) for x in s[1:]]
                d_out[s[0]] = np.array(insert)

            i = np.argsort(d_out['[sample]/[internal standard]'])
            for d in d_out:
                d_out[d] = d_out[d][i]

import sys
sys.path.append(r'C:\Users\willi\Documents\Packages')
from ChromProcess import Classes, file_import
from ChromProcess import curve_fit
import matplotlib.pyplot as plt
import os
import numpy as np

noise_reg = [9.44,9.58]

paths = {}
with open(r"C:\Users\willi\Documents\Data\GCMS\2020_calibrations\Calibration_paths.csv", "r") as f:
    next(f)
    for line in f:
        l = line.strip("\n")
        l = l.split(",")

        paths[l[0]] = [x for x in l[1:] if x != ""]

calibration_parameters = {"Sugar":["A","B","C"]}

all_chroms = []
sig_regions = []
for p in paths:
    for x in paths[p]:
        chroms, cond_file = file_import.load_cdf_from_directory(x, ms = False) # get experimental data and conditions
        all_chroms.extend(chroms)
        c_f = Classes.Information(cond_file)
        [sig_regions.append(c_f.regions) for r in chroms]

noises = np.array([])
sigs = np.array([])
for c,a in enumerate(all_chroms):
    sig_inds = np.array([])
    for R in sig_regions[c]:
        sig_inds = np.hstack((sig_inds,np.where((a.time > R[0])&(a.time < R[1]))[0]))

    noise_inds = np.where((a.time > noise_reg[0])&(a.time < noise_reg[1]))[0]
    sig_inds = sig_inds.astype(int)

    noise_integral = np.std(a.signal[noise_inds], ddof = 1)
    signal_integral = np.std(a.signal[sig_inds], ddof = 1)

    noises = np.hstack((noises,noise_integral))
    sigs = np.hstack((sigs,signal_integral))

SN = sigs/noises
inds = np.where(SN < 10)[0]

for i in inds:
    plt.plot(all_chroms[i].time,all_chroms[i].signal,label = SN[i])
plt.legend()
plt.show()

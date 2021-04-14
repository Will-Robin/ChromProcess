import sys
sys.path.append(r'C:\Users\willi\Documents\Packages')
from ChromProcess import Classes, os, plt, series_builder, plotting, file_output, file_import
from ChromProcess import mass_spectra, calibration_functions, info_params
from ChromProcess import processing_functions
from ChromProcess import series_operations as s_o
import pickle

directory = r"C:\Users\willi\Documents\Data\GCMS\test_data"

calib_file = r"C:\Users\willi\Documents\Packages\info_files\2020_03_16_Combined_Species_Info.csv"

thold = 0.01 # threshold for peak detection as a fraction of the signal maximum
min_inten = 1e4 # to cut off noise if needed.
use_mass_spectra = False
information = Classes.Calibration_File(calib_file, type = "GCMS") # get calibration and assignment information form file.

chroms, cond_file = file_import.load_cdf_from_directory(directory, ms = use_mass_spectra) # get experimental data and conditions

import numpy as np


min_length = 1e100
for c in series.chromatograms:
    if len(c.time) < min_length:
        min_length = len(c.time)
        time = c.time

ser = np.array(series.x_series)
time_region = np.where((ser > 0)&(ser < 12000))[0]

heat_stack = np.empty((len(series.chromatograms),len(series.chromatograms[0].signal[:min_length])))

for c in range(0,len(series.chromatograms)):
    heat_stack[c] = series.chromatograms[c].signal[:min_length]

pca = PCA(n_components=3)

pca.fit(heat_stack)

count = 1
for length, vector in zip(pca.explained_variance_, pca.components_):
    v = vector * 3 * np.sqrt(length)
    plt.plot(time,pca.mean_ + v, label = "pca {} {}%".format(count, round(100*pca.explained_variance_ratio_[count-1],2)), linewidth = 3)
    count +=1

import sys
sys.path.append(r'C:\Users\willi\Documents\Packages')
from ChromProcess import Classes, os, plt, series_builder, plotting, file_output, file_import
from ChromProcess import mass_spectra, calibration_functions, info_params
from ChromProcess import processing_functions
from ChromProcess import series_operations as s_o
import pickle
import numpy as np

directory = r"C:\Users\willi\Documents\Data\GCMS\FRN\FRN094"

calib_file = r"C:\Users\willi\Documents\Packages\info_files\2020_03_16_Combined_Species_Info.csv"

thold = 0.1 # threshold for peak detection as a fraction of the signal maximum
min_inten = 40000 # to cut off noise if needed.
use_mass_spectra = True
information = Classes.Calibration_File(calib_file, type = "GCMS") # get calibration and assignment information form file.

chroms, cond_file = file_import.load_cdf_from_directory(directory, ms = use_mass_spectra, limit = 1) # get experimental data and conditions

ch = chroms[0]





plt.plot(ch.time, ch.signal)
processing_functions.MS_intensity_threshold_chromatogram(ch, threshold = 500)
plt.plot(ch.time, ch.signal)
plt.show()

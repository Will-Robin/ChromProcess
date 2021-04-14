import os
from pathlib import Path
import numpy as np
from ChromProcess import Classes
from ChromProcess import simple_functions as s_f
from ChromProcess import calibration_operations as cal_ops
from ChromProcess import file_import
from NorthNet import Global_formatting
from NorthNet import info_params
import matplotlib.pyplot as plt

# Get paths to calibration data
exp_paths = Classes.DataPaths(r'C:\Users\willi\Documents\Data\2021_calibration_paths.csv')
storage_stem = Path(r'C:\Users\willi\Documents\Data')

calib_output_folder = Path(r'C:\Users\willi\Documents\Data\Analysis_information\GCMS')
calib_dict = {}

IS_positions = []
peak_tables = []
# Import calibration data
for e in exp_paths.paths:
    exp_name = e.experiment_code
    exp_folder = storage_stem/e.data_type/'FRN'/exp_name
    data_folder =exp_folder/'PeakTables'

    for file in os.listdir(data_folder):
        peak_tables.append(Classes.PeakCollection(file = data_folder/file))


for p in peak_tables:
    IS_positions.append(p.internal_standard.retention_time)

print('averages', sum(IS_positions)/len(IS_positions))
plt.hist(IS_positions, bins = len(IS_positions))
plt.show()
quit()
all_peaks = []
for p in peak_tables:
    for pk in p.peaks:
        all_peaks.append(pk.retention_time)

plt.hist(all_peaks, bins  = 500)
plt.show()

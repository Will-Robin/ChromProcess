import os
import numpy as np
from pathlib import Path
from ChromProcess import Classes
from ChromProcess import simple_functions as s_f
from NorthNet import Global_formatting
from ChromProcess import calibration_operations as cal_ops
from ChromProcess import file_import

from NorthNet import info_params
import matplotlib.pyplot as plt

# Get paths to calibration data
exp_paths = Classes.DataPaths(r'C:\Users\willi\Documents\PrebioticDatabase\Data_Information\2021_calibration_paths.csv')
storage_stem = Path(r'C:\Users\willi\Documents\PrebioticDatabase\Data')

peak_positions = {}
peak_weights = {}
# Import calibration data
for e in exp_paths.paths:


    exp_name = e.experiment_code
    sugar = exp_name.split('_')[0]

    data_folder = storage_stem/'GCMS'/'FRN'/exp_name/'PeakTables'

    pos = []
    wt = []
    # read in peak table
    peak_tables = []
    for file in os.listdir(data_folder):
        peak_tables.append(Classes.PeakCollection(file = data_folder/file))

    for p in peak_tables:
        for pk in p.peaks:
            if pk.retention_time == p.internal_standard.retention_time:
                pass
            else:
                pos.append(pk.retention_time)
                wt.append(pk.integral/ p.internal_standard.integral)

    if sugar in peak_positions:
        peak_positions[sugar].extend(pos)
        peak_weights[sugar].extend(wt)
    else:
        peak_positions[sugar] = pos
        peak_weights[sugar] = wt

fig,ax = plt.subplots(nrows = len(peak_positions), sharex = True)
time = np.linspace(5,20, num = 1000)
for c,p in enumerate(peak_positions):

    positions = np.array(peak_positions[p])
    weights = np.array(peak_weights[p])

    nm = info_params.canonical_SMILES[p]
    clr = info_params.colour_assignments[nm]

    cen1 = np.mean(positions)
    sigma1 = s_f.weighted_stdev(positions, weights)
    idx = np.where(weights == np.amax(weights))[0]
    amp1 = weights[idx]
    y_ax = s_f._1gaussian(time, amp1, cen1, sigma1)

    ax[c].hist(positions, bins = len(positions), color = clr, weights = weights)
    ax[c].plot(time, y_ax, color = clr)

plt.show()

import os
from pathlib import Path
from ChromProcess import Classes
import matplotlib.pyplot as plt
from ChromProcess import file_import as f_i
import numpy as np
from ChromProcess import info_params

# Choose the experiment code (see Data_information.csv)
exp_name = 'glucose_1_2021'

# importing information
Path_file = r'C:\Users\willi\Documents\PrebioticDatabase\Data_Information\Data_information.csv'
# path to calibration folder
calib_path = Path(r'C:\Users\willi\Documents\PrebioticDatabase\Analysis_information')
exp_paths = Classes.DataPaths(Path_file)
storage_stem = Path(r'C:\Users\willi\Documents\PrebioticDatabase\Data')
experiment = exp_paths.exp_code_path[exp_name]
data_type = experiment.data_type
exp_path = experiment.path
# State directory in which to store results
store_folder = Path(storage_stem/data_type/'FRN'/exp_name)
analysis_path = store_folder/'{}_analysis_details.csv'.format(exp_name)
chrom_folder = store_folder/'Chromatograms'
# Get calibration files to use
alloc_file = calib_path/'Calibration_file_allocations.csv'
cal_alloc = f_i.importCalibrationFileAllocations(alloc_file)

Calib_file = cal_alloc.GCMS_allocations[exp_name]
calib_file_path = calib_path/'GCMS'/Calib_file
calib = Classes.Instrument_Calibration(file = calib_file_path)
analysis = Classes.Analysis_Information(information_file = analysis_path)

ir_0 = analysis.internal_ref_region[0]
ir_1 = analysis.internal_ref_region[1]
n_plots = 1+len(analysis.regions)
IS_pos = 7.53
fig,ax = plt.subplots(nrows = n_plots, figsize = (10,2*n_plots))
annotated_peaks = []
for file in os.listdir(chrom_folder):
    time, signal = f_i.load_chromatogram_csv(chrom_folder/file)
    idx = np.where((time>ir_0)&(time<ir_1))[0]

    pk = np.where(signal[idx] == np.amax(signal[idx]))[0]
    Ipos = time[idx][pk[0]]

    time = time - Ipos + IS_pos

    ax[0].plot(time[idx], signal[idx], label = file)

    for c,r in enumerate(analysis.regions,1):
        idx2 = np.where((time>r[0])&(time<r[1]))[0]
        ax[c].plot(time[idx2], signal[idx2], c= 'k', alpha = 0.5, zorder = 0)
        for b in calib.boundaries:
            bl = calib.boundaries[b][0]
            bh = calib.boundaries[b][1]
            idx3 = np.where((time>bl)&(time<bh))[0]
            if r[0] < bl and r[1] > bh and len(idx3) > 0:

                sm = info_params.canonical_SMILES[b]
                ax[c].plot(time[idx3], signal[idx3],
                c= info_params.colour_assignments[sm],
                zorder = 1)

                if b not in annotated_peaks:
                    ax[c].annotate(b, xy = (np.average(calib.boundaries[b]),np.amax(signal[idx3])),
                                    xycoords='data',
                                    xytext=(-90, 50), textcoords='offset points',
                                    arrowprops=dict(arrowstyle="->",
                                    connectionstyle="arc,angleA=0,armA=50,rad=10"))
                    annotated_peaks.append(b)



plt.legend()
plt.show()

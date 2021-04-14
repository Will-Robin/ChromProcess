from sklearn.decomposition import PCA

import sys
sys.path.append(r'C:\Users\willi\Documents\Packages')
from ChromProcess import Classes, os, plt, series_builder, plotting, file_output, file_import
from ChromProcess import mass_spectra, calibration_functions, info_params
from ChromProcess import processing_functions
import pickle

from ChromProcess import deconvolution
from ChromProcess import np
from ChromProcess import series_operations


directory = r"C:\Users\willi\Documents\Data\GCMS\2020_calibrations\2020_02_04\xylulose"

calib_file = r"C:\Users\willi\Documents\Packages\info_files\2020_03_16_Combined_Species_Info.csv"

thold = 0.05 # threshold for peak detection as a fraction of the signal maximum
min_inten = 5e4 # to cut off noise if needed.
use_mass_spectra = False
information = Classes.Calibration_File(calib_file, type = "GCMS") # get calibration and assignment information form file.

chroms, cond_file = file_import.load_cdf_from_directory(directory, ms = use_mass_spectra) # get experimental data and conditions
'''
for c in chroms: # plotting chromatograms for inspection.
    plt.plot(c.time,c.signal)
plt.show()
'''

os.chdir("deconvolution")

for f in os.listdir():
    if "conditions" in f:
        cond_file = f

series = Classes.Chromatogram_Series(chroms, cond_file) # Convert chromatogram list and conditions into a Series object

series_operations.pick_peaks(series, threshold = 0.1, max_intensity = 1e100, min_intensity = min_inten)
os.makedirs("plots",exist_ok = True)
os.chdir("plots")

series_operations.get_internal_ref_integrals(series)


series_integrals = []
stacking = np.array([])

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

reg = series.regions[0]
inds = np.where((time > reg[0])&(time < reg[1]))[0]
time = time[inds]

heat_stack =  heat_stack[time_region]
stacking = heat_stack[:,inds]

pca = PCA(n_components=3)

pca.fit(stacking)

for x in range(0,len(stacking)):
    if x == 0:
        plt.plot(time, stacking[x], label = "signal", c = "k", alpha = 0.5)
    else:
        plt.plot(time, stacking[x], c = "k", alpha = 0.5)

count = 1
for length, vector in zip(pca.explained_variance_, pca.components_):
    v = vector * 3 * np.sqrt(length)
    plt.plot(time,pca.mean_ + v, label = "pca {} {}%".format(count, round(100*pca.explained_variance_ratio_[count-1],2)), linewidth = 3)
    count +=1

v1 = np.dot( pca.components_[0] * 3 * np.sqrt(pca.explained_variance_[0]),  pca.components_[2] * 3 * np.sqrt(pca.explained_variance_[2]) )
plt.plot(time,  pca.mean_ + v1, label = "{}".format("LC 1, 3"), linewidth = 3, c = "r")

plt.legend()
plt.savefig("PCA_vectors_{}_to_{}.png".format(reg[0],reg[1]))


plotting.region_heatmap(series, information)

import sys
sys.path.append(r'C:\Users\willi\Documents\Packages')
from ChromProcess import Classes, os, plt
from ChromProcess import mass_spectra, info_params
from ChromProcess import processing_functions
from ChromProcess import series_operations
import pickle

directory = r"C:\Users\willi\Documents\Data\GCMS\FRN\FRN088\B"
chromatogram_number = 4

calib_file = r"C:\Users\willi\Documents\Packages\info_files\2020_03_16_Combined_Species_Info.csv" # 3/16 vs 2/18
IS_pos = 6.74 # 6.925/ 6.74

thold = 0.01 # threshold for peak detection as a fraction of the signal maximum
min_inten = 1000 # to cut off noise if needed.
use_mass_spectra = True
information = Classes.Calibration_File(calib_file, type = "GCMS") # get calibration and assignment information form file.
os.chdir(directory)
f_list = os.listdir()
c_names = []
for f in f_list:
    if ".cdf" in f:
        c_names.append(f)
    if "conditions" in f:
        cond_file = Classes.Information(f)

c_names.sort()
chrom = Classes.Chromatogram(c_names[chromatogram_number-1], mass_spec = use_mass_spectra)

import numpy as np
processing_functions.MS_intensity_threshold_chromatogram(chrom, threshold = 500)

inds = np.where((chrom.time > cond_file.internal_ref_region[0]) & (chrom.time < cond_file.internal_ref_region[1]))[0]
max = np.amax(chrom.signal[inds])
inds = chrom.signal == max
time_base = chrom.time[inds]

plt.plot(chrom.time, chrom.signal, alpha = 0.5, c = "k")
chrom.time = chrom.time - time_base
chrom.time = chrom.time + IS_pos
plt.plot(chrom.time, chrom.signal, c= "k")
plt.show()

os.makedirs("inspect_regions", exist_ok = True)
os.chdir("inspect_regions")

for b in cond_file.regions:
    fig = plt.figure(figsize = (info_params.across + 4,info_params.up+5))
    ax = plt.subplot(111)
    data_max = 0
    data_min = 1e100
    time, ion_dict = mass_spectra.ion_chromatogram_region(chrom, b[0], b[1])
    series_operations.bin_dictionary(ion_dict, bound = 0.1, rounding = 2)

    for i in ion_dict:
        if np.amax(ion_dict[i]) < 500:
            pass
        else:
            if int(i) in info_params.frag_colours:
                ax.plot(time, ion_dict[i], c= info_params.frag_colours[int(i)])
            else:
                ax.plot(time, ion_dict[i])

        if np.amax(ion_dict[i]) > data_max:
            data_max = np.amax(ion_dict[i])
        if np.amin(ion_dict[i]) < data_min:
            data_min =np.amin(ion_dict[i])

    ylim = [data_min,data_max]
    for g in information.boundaries:
        region_inds = np.where((chrom.time > b[0]) & (chrom.time < b[1]))[0]

        if information.boundaries[g][1] < b[1] and information.boundaries[g][0] > b[0]:

            lower_bound_y = np.linspace(ylim[0],ylim[1], num = 100)
            lower_bound_x = np.full(len(lower_bound_y),information.boundaries[g][0])

            upper_bound_y = np.linspace(ylim[0],ylim[1], num = 100)
            upper_bound_x = np.full(len(upper_bound_y),information.boundaries[g][1])

            ax.scatter(lower_bound_x, lower_bound_y, marker ="|", c= "g", alpha = 0.5, edgecolors = None, s = 10)
            ax.scatter(upper_bound_x, upper_bound_y, marker ="|", c= "r", alpha = 0.5, edgecolors = None, s = 10)
            ax.text(np.average(information.boundaries[g]),ylim[1],g, fontsize = 10, rotation = 90, rotation_mode = "anchor", horizontalalignment = "center", verticalalignment = "center")


    plt.savefig("{}_{}.png".format(chrom.filename, np.average(b)))
    plt.close()


os.chdir("..")

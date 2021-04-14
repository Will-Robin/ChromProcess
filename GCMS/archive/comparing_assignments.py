import sys
sys.path.append(r'C:\Users\willi\Documents\Packages')
from ChromProcess import Classes, os, plt, series_builder, plotting, file_output, file_import
from ChromProcess import mass_spectra, calibration_functions, info_params
from ChromProcess import processing_functions
import pickle

directory = r"C:\Users\willi\Documents\Data\GCMS\Collecting_chroms"

calib_file = r"C:\Users\willi\Documents\Packages\info_files\2020_03_16_Combined_Species_Info.csv"

thold = 0.01 # threshold for peak detection as a fraction of the signal maximum
min_inten = 1e4 # to cut off noise if needed.
use_mass_spectra = False
information = Classes.Calibration_File(calib_file, type = "GCMS") # get calibration and assignment information form file.

os.chdir(directory) # changes the directory to where files are stored
filelist = os.listdir() # get list of files in directory
filelist.sort() # sort the files

chroms = []

names = ["FRN071_039.cdf", # 0
         "FRN073_050.cdf", # 1
         "FRN077_050.cdf", # 2
         "FRN087_050.cdf", # 3
         "FRN088_050.cdf", # 4
         "FRN089_030.cdf", # 5
         "FRN090_030.cdf", # 6
         "FRN091_010.cdf", # 7
         "FRN092_078.cdf", # 8
         "FRN093_057.cdf"] # 9

filelist = [names[4], names[5]]
time_add = [0, -0.01]

for f in names:
    if f.endswith(".cdf"):
        print("loading", f)
        chroms.append(Classes.Chromatogram(f, mass_spec = use_mass_spectra))

cond_file = "collection_conditions.csv"

fig = plt.figure()
ax = fig.gca(projection='3d')
import numpy as np
for count,c in enumerate(chroms): # plotting chromatograms for inspection.
    inds = np.where((c.time > 9.5)&(c.time < 11.5))[0]
    ax.plot(c.time[inds],np.full(len(c.signal[inds]), count),c.signal[inds]/1000000)
ax.xaxis.set_pane_color((1.0, 1.0, 1.0, 0.0))
ax.yaxis.set_pane_color((1.0, 1.0, 1.0, 0.0))
ax.zaxis.set_pane_color((1.0, 1.0, 1.0, 0.0))
ax.xaxis._axinfo["grid"]['color'] =  (1,1,1,0)
ax.yaxis._axinfo["grid"]['color'] =  (1,1,1,0)
ax.zaxis._axinfo["grid"]['color'] =  (1,1,1,0)
ax.set_xlabel("retention time/ min.")
ax.set_ylabel("chrom_index/ n")
ax.set_zlabel("signal/ 10$^6$.counts")
plt.show()

quit()

# Call a series of functions on the Series object.
series_builder.series_from_chromatogram_set(series, information,
                                            threshold = thold,
                                            min_intensity = min_inten,
                                            plot_mass_spectra = use_mass_spectra,
                                            label_assignments = False,
                                            combine_assignments = False,
                                            internal_ref_correct = False)

os.makedirs("ion_chromatograms", exist_ok = True)
os.chdir("ion_chromatograms")
plotting.plot_peak_ion_chromatogram_from_series(series)
plotting.write_peak_ion_chromatogram_from_series(series)
#plotting.plot_ion_channels(series)
mass_spectra.ion_chromatograms_relative_to_peak(series)
print("Finished.")

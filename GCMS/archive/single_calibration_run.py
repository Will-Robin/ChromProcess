import sys
sys.path.append(r'C:\Users\willi\Documents\Packages')
from ChromProcess import Classes, os, plt, series_builder, plotting, file_output, file_import
from ChromProcess import mass_spectra, calibration_functions, info_params
from ChromProcess.series_operations import cluster
from ChromProcess import processing_functions
import numpy as np

path_list = r"C:\Users\willi\Documents\Data\GCMS\2020_calibrations\Calibration_paths.csv"

path_dict = {}
wd= os.getcwd()
with open(path_list, "r") as f:
    for c,line in enumerate(f):
        if c > 0:
            x = line.strip("\n").split(",")
            path_dict[x[0]]= [z for z in x[1:] if z != ""]

calib_file = r"C:\Users\willi\Documents\Packages\info_files\2020_03_16_Combined_Species_Info.csv"
information = Classes.Calibration_File(calib_file, type = "GCMS") # get calibration and assignment information from file.
information.calibrations = {"none":[0,1]}
information.boundaries = {"none":[0,1]}

thold = 0.1 # threshold for peak detection as a fraction of the signal maximum
min_inten = 10000 # to cut off noise if needed.
use_mass_spectra = True
#path_dict = {"X": "/Users/williamrobinson/Documents/Nijmegen/2020_calibrations/2020_03_12/threose/C"}
directory = r"C:\Users\willi\Documents\Data\GCMS\2020_calibrations\2020_03_12\threose\C"
chroms, cond_file = file_import.load_cdf_from_directory(directory, ms = use_mass_spectra) # get experimental data and conditions

for c in chroms: # plotting chromatograms for inspection.
    processing_functions.MS_intensity_threshold_chromatogram(c, threshold = 500)
    plt.plot(c.time,c.signal)
plt.show()

series = Classes.Chromatogram_Series(chroms, cond_file) # Convert chromatogram list and conditions into a Series object

calibration_functions.calibration_series(series) # Convert the series into a calibration series (since experiment info is slightly different)

# Call a series of functions on the Series object.
series_builder.series_from_chromatogram_set(series, information,
                                            threshold = thold,
                                            min_intensity = min_inten,
                                            plot_mass_spectra = use_mass_spectra,
                                            label_assignments = False,
                                            combine_assignments = False)

os.makedirs("Calibration_output", exist_ok = True)
os.chdir("Calibration_output")

cabs_1 = calibration_functions.fit_calibration_curve_TIC(series) # Get calibration factors from the total ion chromatogram.
cabs_2 = calibration_functions.fit_calibration_curve_IC(series) # Get calibration factors from ion chromatograms.        print(series.ion_series)

# Converting integral series from ion chromatogram to series
# corresponding to the intergal of m/z = 73 ion chromatograms
for peak_pos in series.ion_series:
    for i in series.ion_series[peak_pos]:
        if int(i) == 73:
            series.integral_series[peak_pos] = series.ion_series[peak_pos][i]
        else:
            pass
os.chdir("..")
os.makedirs("ion_chromatograms", exist_ok = True)
os.chdir("ion_chromatograms")
plotting.plot_peak_ion_chromatogram_from_series(series)
plotting.write_peak_ion_chromatogram_from_series(series)
#plotting.plot_ion_channels(series)
mass_spectra.ion_chromatograms_relative_to_peak(series)
print("Finished.")

import sys
sys.path.append(r'C:\Users\willi\Documents\Packages')
from ChromProcess import Classes, os, plt, series_builder, plotting, file_output, file_import
from ChromProcess import mass_spectra, calibration_functions, info_params
from ChromProcess import processing_functions
import pickle

os.chdir(r"C:\Users\willi\Documents\Data\GCMS\FRN\FRN089")
calib_file = r"C:\Users\willi\Documents\Packages\info_files\2020_03_16_Combined_Species_Info.csv"

with open("FRN089", "rb") as f:
    series = pickle.load(f)

information = Classes.Calibration_File(calib_file, type = "GCMS") # get calibration and assignment information form file.

# Call a series of functions on the Series object.
series_builder.series_from_chromatogram_set(series, information,
                                            threshold = thold,
                                            min_intensity = min_inten,
                                            plot_mass_spectra = use_mass_spectra,
                                            label_assignments = False)



plotting.plot_peak_ion_chromatogram_from_series(series)
#plotting.plot_ion_channels(series)
mass_spectra.ion_chromatograms_relative_to_peak(series)
print("Finished.")

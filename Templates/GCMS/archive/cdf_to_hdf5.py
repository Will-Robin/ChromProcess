import sys
sys.path.append(r'C:\Users\willi\Documents\Packages')
from ChromProcess import Classes, os, plt, series_builder, plotting, file_output, file_import
from ChromProcess import mass_spectra, calibration_functions, info_params
from ChromProcess import processing_functions
import pickle



directory = r"C:\Users\willi\Documents\Data\GCMS\FRN\FRN090"

calib_file = r"C:\Users\willi\Documents\Packages\info_files\2020_03_16_Combined_Species_Info.csv"

thold = 0.1 # threshold for peak detection as a fraction of the signal maximum
min_inten = 3e4 # to cut off noise if needed.
use_mass_spectra = True
information = Classes.Calibration_File(calib_file, type = "GCMS") # get calibration and assignment information form file.

chroms, cond_file = file_import.load_cdf_from_directory(directory, ms = use_mass_spectra) # get experimental data and conditions

series = Classes.Chromatogram_Series(chroms, cond_file) # Convert chromatogram list and conditions into a Series object

for c in chroms: # plotting chromatograms for inspection.
    #c.signal = processing_functions.savitzky_golay(c.signal,3,1)
    #plt.plot(c.time,c.signal)
    file_output.store_chromatogram(c)

    cr = Classes.Chromatogram("{}.hdf5".format(c.filename))

    plt.plot(cr.time, cr.signal)

plt.show()

# Call a series of functions on the Series object.
series_builder.series_from_chromatogram_set(series, information,
                                            threshold = thold,
                                            min_intensity = min_inten,
                                            plot_mass_spectra = use_mass_spectra,
                                            label_assignments = False,
                                            combine_assignments = True)

plotting.plot_peak_ion_chromatogram_from_series(series)
plotting.write_peak_ion_chromatogram_from_series(series)
#plotting.plot_ion_channels(series)
mass_spectra.ion_chromatograms_relative_to_peak(series)
print("Finished.")

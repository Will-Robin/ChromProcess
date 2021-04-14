import sys
sys.path.append(r'C:\Users\willi\Documents\Packages')
from ChromProcess import Classes, os, plt, series_builder, plotting, file_output, file_import
from ChromProcess import mass_spectra, calibration_functions, info_params
from ChromProcess import processing_functions
from ChromProcess import series_operations
import pickle

directory = r"C:\Users\willi\Documents\Data\GCMS\FRN\FRN073"

calib_file = r"C:\Users\willi\Documents\Packages\info_files\2020_03_16_Combined_Species_Info.csv"

thold = 0.01 # threshold for peak detection as a fraction of the signal maximum
min_inten = 500 # to cut off noise if needed.
max_intensity = 1e100
use_mass_spectra = True
information = Classes.Calibration_File(calib_file, type = "GCMS") # get calibration and assignment information form file.

chroms, cond_file = file_import.load_cdf_from_directory(directory, ms = use_mass_spectra) # get experimental data and conditions

series = Classes.Chromatogram_Series(chroms, cond_file) # Convert chromatogram list and conditions into a Series object

for c in chroms: # plotting chromatograms for inspection.
    processing_functions.MS_intensity_threshold_chromatogram(c, threshold = 500)
    plt.plot(c.time,c.signal)
plt.show()

series_operations.get_internal_ref_integrals(series)

# Find peaks in each chromatogram, region by region
series_operations.pick_peaks(series, threshold = thold, max_intensity = max_intensity, min_intensity = min_inten)

series_operations.get_integrals(series) # Adds intergral information into the picked peaks

os.makedirs(r"Peak_tables", exist_ok = True)
os.chdir(r"Peak_tables")
def write_peak_table(chromatogram, value = 0.0, series_unit = "None"):
    with open("{}_peak_table.csv".format(chromatogram.filename), "w") as f:
        f.write("{},{}\n".format(series_unit,value))
        f.write("Retention_time/ min, I/Istand, peak start/ min, peak end/ min\n")
        st_ind = chromatogram.internal_reference.indices[0]
        end_ind = chromatogram.internal_reference.indices[-1]
        f.write("{},{},{},{}\n".format(chromatogram.internal_reference.retention_time, chromatogram.internal_reference.integral, chromatogram.time[st_ind], chromatogram.time[end_ind]))
        for p in chromatogram.peaks:
            st_ind = chromatogram.peaks[p].indices[0]
            end_ind = chromatogram.peaks[p].indices[-1]
            f.write("{},{},{},{}\n".format(chromatogram.peaks[p].retention_time, chromatogram.peaks[p].integral, chromatogram.time[st_ind], chromatogram.time[end_ind]))

for count,c in enumerate(series.chromatograms):
    write_peak_table(c, value = series.x_series[count],series_unit = "Time point/ s")

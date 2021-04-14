import sys
sys.path.append(r'C:\Users\willi\Documents\Packages')
from ChromProcess import Classes, os, plt, series_builder, plotting, file_output, file_import
from ChromProcess import mass_spectra, calibration_functions, info_params
from ChromProcess import processing_functions
import pickle

directory = r"C:\Users\willi\Documents\Data\GCMS\FRN\FRN103"

calib_file = r"C:\Users\willi\Documents\Packages\info_files\2020_03_16_Combined_Species_Info.csv"

thold = 0.1 # threshold for peak detection as a fraction of the signal maximum
min_inten = 4000 # to cut off noise if needed.
use_mass_spectra = True
information = Classes.Calibration_File(calib_file, type = "GCMS") # get calibration and assignment information form file.

chroms, cond_file = file_import.load_cdf_from_directory(directory, ms = use_mass_spectra) # get experimental data and conditions

series = Classes.Chromatogram_Series(chroms, cond_file) # Convert chromatogram list and conditions into a Series object
import numpy as np

for c in chroms: # plotting chromatograms for inspection.

    processing_functions.MS_intensity_threshold_chromatogram(c, threshold = 500)

    inds = np.where((c.time > series.internal_ref_region[0]) & (c.time < series.internal_ref_region[1]))[0]
    max = np.amax(c.signal[inds])
    inds = c.signal == max
    time_base = c.time[inds]

    c.time = c.time - time_base
    c.time = c.time + 6.74

    plt.plot(c.time,c.signal)

plt.show()

os.makedirs("ion_timecourses", exist_ok = True)
os.chdir("ion_timecourses")

region = [10.6,10.8]

for c in chroms:

    t, ion_dict = mass_spectra.ion_chromatogram_region(c, region[0], region[1])

    for i in ion_dict:
        if int(i) == 73:
            pass
        elif np.amax(ion_dict[i]) < 4000:
            pass
        elif int(i) in info_params.frag_colours:
            plt.plot(t,ion_dict[i], c=  info_params.frag_colours[int(i)])
            plt.annotate(i, xy = (10.6,np.amax(ion_dict[i])))
        else:
            plt.plot(t,ion_dict[i])
            plt.annotate(i, xy = (10.6,np.amax(ion_dict[i])))

    plt.savefig("{}.png".format(c.filename))
    plt.close()
print("Finished.")

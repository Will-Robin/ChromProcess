import sys
sys.path.append(r'/Users/williamrobinson/documents/nijmegen/packages')

import numpy as np
from ChromProcess import Classes, os, plt, series_builder, plotting, file_output, file_import
from ChromProcess import mass_spectra, calibration_functions, info_params
from ChromProcess import processing_functions
import pickle

directory = r"/Users/williamrobinson/Documents/Nijmegen/Data/GCMS/FRN107/B"

calib_file, IS_pos = file_import.get_calibration_file_allocations(r'/Users/williamrobinson/Documents/Nijmegen/Dynamic_environment_project/Calibration_file_allocations.csv', directory.split('/')[-2])

if calib_file != None:
    calib_file = r"/Users/williamrobinson/Documents/Nijmegen/Packages/info_files/{}".format(calib_file)
else:
    calib_file =  r"/Users/williamrobinson/Documents/Nijmegen/Packages/info_files/{}".format('2020_03_16_Combined_Species_Info.csv')
#IS_pos = 6.975 # 6.956/6.925/ 6.74

thold = 0.01 # threshold for peak detection as a fraction of the signal maximum
min_inten = 1000 # to cut off noise if needed.
use_mass_spectra = True
read_modified_bounds = True
chroms, cond_file = file_import.load_cdf_from_directory(directory, ms = use_mass_spectra) # get experimental data and conditions
series = Classes.Chromatogram_Series(chroms, cond_file) # Convert chromatogram list and conditions into a Series object

if read_modified_bounds:
    try:
        modified_bounds = file_import.read_local_assignments('{}_local_assignments.csv'.format(series.set_name))
    except:
        modified_bounds = {'' : [0,1],
                           "set_IS_pos" : [IS_pos]}
else:
    modified_bounds = {'' : [0,1],
                       "set_IS_pos" : [IS_pos]}

information = Classes.Calibration_File(calib_file, type = "GCMS") # get calibration and assignment information form file.

temp_dict = {}
for m in modified_bounds:
    if m == "set_IS_pos":
        pass
    else:
        temp_dict[m] = modified_bounds[m]

for i in information.boundaries:
    if i in temp_dict:
        pass
    else:
         temp_dict[i] = information.boundaries[i]

information.boundaries = {}

for t in temp_dict:
    information.boundaries[t] = temp_dict[t]

for c in chroms: # plotting chromatograms for inspection.

    processing_functions.MS_intensity_threshold_chromatogram(c, threshold = 500)
    if modified_bounds['set_IS_pos'] == None:
        pass
    else:
        inds = np.where((c.time > series.internal_ref_region[0]) & (c.time < series.internal_ref_region[1]))[0]
        max = np.amax(c.signal[inds])
        inds = c.signal == max
        time_base = c.time[inds]

        c.time = c.time - time_base
        c.time = c.time + modified_bounds['set_IS_pos']

'''    plt.plot(c.time,c.signal)

plt.show()'''

# Call a series of functions on the Series object.
series_builder.series_from_chromatogram_set(series, information,
                                            threshold = thold,
                                            min_intensity = min_inten,
                                            plot_mass_spectra = use_mass_spectra,
                                            label_assignments = False,
                                            combine_assignments = False,
                                            cluster_dev = 0.02)

os.makedirs("heatmap_ion_chromatograms", exist_ok = True)
os.chdir("heatmap_ion_chromatograms")
plotting.region_heatmap_with_ion_chroms(series, information, cutoff = 0.05)
os.chdir('..')
'''
os.makedirs("ion_chromatograms", exist_ok = True)
os.chdir("ion_chromatograms")
plotting.plot_peak_ion_chromatogram_from_series(series, information)
#plotting.write_peak_ion_chromatogram_from_series(series)
#plotting.plot_ion_channels(series)
mass_spectra.ion_chromatograms_relative_to_peak(series)
os.chdir("..")
'''
file_output.write_modified_calibration_params(modified_bounds, series)

print("Finished.")

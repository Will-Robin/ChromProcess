import sys
sys.path.append(r'C:\Users\willi\Documents\Packages')
from ChromProcess import Classes, os, plt, series_builder, plotting, file_output, file_import
from ChromProcess import mass_spectra, calibration_functions, info_params
from ChromProcess import processing_functions, series_operations
import pickle
import numpy as np

directory = r"C:\Users\willi\Documents\Data\GCMS\FRN\FRN103"

chroms, cond_file = file_import.load_cdf_from_directory(directory, ms = True, limit = 1) # get experimental data and conditions

for c in chroms:
    time, ion_dict = mass_spectra.ion_chromatogram_region(c, 9, 15)

    sorted_masses = sorted([*ion_dict])
    clusters = []
    for clust in series_operations.cluster(sorted_masses, bound = 0.5):
        clusters.append(clust)

    out_log = {}
    for clust in clusters:

        position = round(np.average(clust),2)

        out_log[position] = []

        for o in range(0,len(sorted_masses)):
            if sorted_masses[o] in clust:
                if len(out_log[position]) == 0:
                    out_log[position] = ion_dict[ sorted_masses[o] ]
                else:
                    for count,p in enumerate(out_log[position]):
                        if p == 0:
                            out_log[position][count] = ion_dict[ sorted_masses[o] ][count]
                        else:
                            pass

    X = time
    Y = [*out_log]

    fig, ax = plt.subplots(figsize=(15,10))

    for count, o in enumerate(out_log):
        ax.plot(time, out_log[o]+count*500)

    plt.show()
    plt.close()

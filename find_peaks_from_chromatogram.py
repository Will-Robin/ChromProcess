#from pathlib import Path
import os
import matplotlib.pyplot as plt
from ChromProcess.Loading import chrom_from_csv
from ChromProcess.Loading import analysis_from_csv
from ChromProcess.Loading import conditions_from_csv
experiment_number = 'FRN140'
experiment_folder = r"C:\Users\thijs\Documents\PhD\Data\FRN140"
from ChromProcess.Utils.peak_finding import find_peaks_scipy
from ChromProcess.Utils import indices_from_boundary, peak_indices_to_times

from ChromProcess.Processing import add_peaks_to_chromatogram
from ChromProcess.Processing import integrate_chromatogram_peaks
from ChromProcess.Processing import internal_standard_integral

chromatogram_directory = f'{experiment_folder}\ChromatogramCSV'
conditions_file = f'{experiment_folder}\{experiment_number}_conditions.csv'
analysis_file = f'{experiment_folder}\{experiment_number}_analysis_details.csv'
peak_collection_directory = f'{experiment_folder}\PeakCollections'


conditions = conditions_from_csv(conditions_file)
analysis = analysis_from_csv(analysis_file)

os.makedirs(peak_collection_directory, exist_ok = True)
chromatogram_files = os.listdir(chromatogram_directory)
chromatogram_files.sort()
chroms = []
for f in chromatogram_files:
    chroms.append(chrom_from_csv(f'{chromatogram_directory}/{f}'))

fig, ax = plt.subplots()
for c in chroms:
    ax.plot(c.time, c.signal, label = c.filename)

#plt.show()

is_start = analysis.internal_standard_region[0]
is_end = analysis.internal_standard_region[1]
for c in chroms:
    internal_standard_integral(c, is_start, is_end)

threshold = analysis.peak_pick_threshold
if type(threshold) == float:
    threshold = [threshold for r in analysis.regions]
for chrom in chroms:
    for reg,thres in zip(analysis.regions,threshold):
        inds = indices_from_boundary(chrom.time, reg[0], reg[1])
        time = chrom.time[inds]
        signal = chrom.signal[inds]
        picked_peaks = find_peaks_scipy(signal, 
                        threshold=thres, 
                        min_dist=1, 
                        max_inten = 1e100, 
                        prominence = 5000, 
                        wlen = 1001, 
                        look_ahead = 12
                        )
        peak_features = peak_indices_to_times(time,picked_peaks)
        add_peaks_to_chromatogram(peaks, chrom)
        integrate_chromatogram_peaks(chrom)


for c,v in zip(chroms, conditions.series_values):
    c.write_peak_collection(filename = f'{peak_collection_directory}/{c.filename}',
                        header_text = f"{conditions.series_unit},{v}\n",
                        )
import time
import timeit
from ChromProcess import chromate
from ChromProcess.Loading import chrom_from_csv
from ChromProcess.Utils.signal_processing import signal_processing as sig
from scipy.signal import find_peaks, peak_widths
import matplotlib.pyplot as plt

file = "Tutorials/example_data/Example/ExperimentalData/ChromatogramCSV/chrom_001.csv"

chrom = chrom_from_csv(file)

smoothed = sig.adjacent_average(chrom.signal, 10)

start = time.time()
peaks, _ = find_peaks(smoothed, height=chrom.signal.max()*0.1)
results_half = peak_widths(smoothed, peaks, rel_height=0.5)
end = time.time()
print(f"Scipy benchmark {(end-start)*1000} ms")

print("next")
sig_list = smoothed.tolist()
start = time.time()
peaks, starts, ends = chromate.find_peaks(sig_list, chrom.signal.max()*0.1)
end = time.time()
print(f"Took {(end-start)*1000} ms")

fig, ax = plt.subplots()

ax.plot(chrom.time, chrom.signal)
ax.scatter(chrom.time[peaks], chrom.signal[peaks], c= "k")
ax.scatter(chrom.time[starts], chrom.signal[starts], c= "g")
ax.scatter(chrom.time[ends], chrom.signal[ends], c= "r")
plt.show()


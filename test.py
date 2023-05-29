import time
from chromate import find_peaks_wrapper
from ChromProcess.Loading import chrom_from_csv
from ChromProcess.Utils.signal_processing import signal_processing as sig
import matplotlib.pyplot as plt

file = "Tutorials/example_data/Example/ExperimentalData/ChromatogramCSV/chrom_001.csv"

chrom = chrom_from_csv(file)

smoothed = sig.adjacent_average(chrom.signal, 5)
smoothed = sig.adjacent_average(chrom.signal, 10)

start = time.time()
peaks, starts, ends = find_peaks_wrapper(smoothed, chrom.signal.max()*0.1)
end = time.time()
print(f"Took {(end-start)*1000} ms")

fig, ax = plt.subplots()

ax.plot(chrom.time, chrom.signal)
ax.scatter(chrom.time[peaks], chrom.signal[peaks], c= "k")
ax.scatter(chrom.time[starts], chrom.signal[starts], c= "g")
ax.scatter(chrom.time[ends], chrom.signal[ends], c= "r")
plt.show()


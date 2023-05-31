import time
from ChromProcess import chromate
from ChromProcess.chromate import find_peaks
import numpy as np
from ChromProcess.Loading import chrom_from_csv
from ChromProcess.Utils.signal_processing import signal_processing as sig
import matplotlib.pyplot as plt

file = "Tutorials/example_data/Example/ExperimentalData/ChromatogramCSV/chrom_001.csv"
chrom = chrom_from_csv(file)

smoothed = sig.adjacent_average(chrom.signal, 10)

peaks, starts, ends = find_peaks(smoothed, chrom.signal.max()*0.1)

for (p,s,e) in zip(peaks, starts, ends):
    time = chrom.time[s:e]
    sig = chrom.signal[s:e]

    np_result = np.trapz(sig, x=time)
    chromate_result = chromate.trapz(sig, time)

    print(f"numpy result:    {np_result}")
    print(f"chromate result: {chromate_result}")
    print()

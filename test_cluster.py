import time
from ChromProcess import chromate
from ChromProcess.chromate import find_peaks
from ChromProcess.Loading import chrom_from_csv
from ChromProcess.Utils.signal_processing import signal_processing as sig
import matplotlib.pyplot as plt

file = "Tutorials/example_data/Example/ExperimentalData/ChromatogramCSV/chrom_001.csv"
chrom = chrom_from_csv(file)
smoothed = sig.adjacent_average(chrom.signal, 10)
peaks, starts, ends = find_peaks(smoothed, chrom.signal.max() * 0.001)

peak_list = [chrom.time[x] for x in peaks]

chromate_clusters = [c for c in chromate.cluster(peak_list, 0.1)]
chromate_clusters_idx = [c for c in chromate.cluster_indices(peak_list, 0.1)]


for b,c in zip(chromate_clusters, chromate_clusters_idx):
    print(f"chromate_clusters:     {b}")
    print(f"chromate_clusters idx: {c}")
    print()

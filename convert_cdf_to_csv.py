#from pathlib import Path
import os
from ChromProcess.Loading import chrom_from_cdf
import matplotlib.pyplot as plt
from ChromProcess.Loading import analysis_from_csv
from ChromProcess.Processing import ic_background_subtraction

experiment_folder = r"C:\Users\thijs\Documents\PhD\Data\FRN141"

chromatogram_directory = f'{experiment_folder}\cdf'

chromatogram_output_folder = f'{experiment_folder}/ChromatogramCSV'
os.makedirs(chromatogram_output_folder,exist_ok=True)


chromatogram_files = os.listdir(chromatogram_directory)
chromatogram_files.sort()
chroms = []
for f in chromatogram_files:
    chroms.append(chrom_from_cdf(f'{chromatogram_directory}/{f}',load_ms=True))


for count,c in enumerate(chroms):
    chroms[count].signal = ic_background_subtraction(c, threshold = 500)


fig, ax = plt.subplots()
for c in chroms:
    ax.plot(c.time, c.signal, label = c.filename)
plt.show()

for c in chroms:
    c.write_to_csv(filename = f'{chromatogram_output_folder}/{c.filename}')


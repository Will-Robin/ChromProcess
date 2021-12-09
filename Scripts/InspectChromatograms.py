'''
Loading and plotting chromatograms.
'''

import os
import matplotlib.pyplot as plt

from ChromProcess import Classes

chromatogram_dir = '../Examples/Example/ExperimentalData/ExampleChromatograms'

chromatogram_file_list = os.listdir(chromatogram_dir)
chromatogram_file_list.sort()

chroms = []

for f in chromatogram_file_list:
    chroms.append(Classes.Chromatogram(file = f'{chromatogram_dir}/{f}'))

fig, ax = plt.subplots()
for c in chroms:
    ax.plot(c.time, c.signal, label = c.filename)

plt.show()

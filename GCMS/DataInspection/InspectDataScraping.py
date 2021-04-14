from ChromProcess import file_import as f_i
import matplotlib.pyplot as plt
from pathlib import Path
import os
import numpy as np
from ChromProcess import Classes

exp_name = 'allose_1_2021'
method = 'GCMS'

rt_fldr = r'C:\Users\willi\Documents\Data\{}\FRN\{}'.format(method,exp_name)
root_data_folder = Path(rt_fldr)
chrom_folder = root_data_folder/'Chromatograms'
drep_folder = root_data_folder/'PeakTables'
chroms = []
for file in os.listdir(chrom_folder):
    if not file.startswith('.'):
        print(file)
        chroms.append(f_i.load_chromatogram_csv(chrom_folder/file))

# read in peak table
peak_tables = []
for file in os.listdir(drep_folder):
    print(file)
    peak_tables.append(Classes.PeakCollection(file = drep_folder/file))

analysis_path = root_data_folder/'{}_analysis_details.csv'.format(exp_name)
analysis = Classes.Analysis_Information(information_file = analysis_path)

for c in chroms:
    plt.plot(c[0],c[1])

p_scatter_x = []
begin_scatter_x = []
end_scatter_x = []
for c,p in enumerate(peak_tables):
    for pk in p.peaks:
        idx = np.where(chroms[c][0]==pk.retention_time)[0]
        plt.scatter(chroms[c][0][idx], chroms[c][1][idx])
        plt.scatter(pk.retention_time, 200000)

for r in analysis.regions:
    plt.scatter(r[0], 200000,marker = '|', s = 10, c = 'g')
    plt.scatter(r[1], 200000, marker = '|', s = 10, c= 'r')
plt.show()

import os
from pathlib import Path
from ChromProcess import Classes
import matplotlib.pyplot as plt
from ChromProcess import file_import as f_i

# Choose the experiment code (see Data_information.csv)
exp_name = 'glucose_1_2021'

# importing information
Path_file = r'C:\Users\willi\Documents\Data\Data_information.csv'
exp_paths = Classes.DataPaths(Path_file)
storage_stem = Path(r'C:\Users\willi\Documents\Data')
experiment = exp_paths.exp_code_path[exp_name]
data_type = experiment.data_type
exp_path = experiment.path
# State directory in which to store results
store_folder = Path(storage_stem/data_type/'FRN'/exp_name/'Chromatograms')

for file in os.listdir(store_folder):
    time, signal = f_i.load_chromatogram_csv(store_folder/file)
    plt.plot(time, signal, label = file)
plt.legend()
plt.show()

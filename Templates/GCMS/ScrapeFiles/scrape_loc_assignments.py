from ChromProcess import Classes
import shutil
from pathlib import Path
import os

GCMS_root = Path(r'C:\Users\willi\Documents\Data\GCMS\FRN')
exp_paths = Classes.DataPaths(r'C:\Users\willi\Documents\Data\Data_information.csv')
exp_names = [e.experiment_code for e in exp_paths.paths]
print(exp_names)
storage_stem = Path(r'C:\Users\willi\Documents\Data')

file_dict = {}
for folder in os.listdir(GCMS_root):
    for (root,dirs,files) in os.walk(GCMS_root/folder):

        for f in files:
            if 'local_assignments' in f:
                file_dict[f] = root

store_folder = Path(storage_stem/'GCMS'/'FRN')

for f in file_dict:
    stem = Path(file_dict[f])

    exp_name = f.split('_')[0]

    destination = store_folder/exp_name

    if not os.path.exists(destination/f) and exp_name in exp_names:
        shutil.copyfile(stem/f, destination/f)

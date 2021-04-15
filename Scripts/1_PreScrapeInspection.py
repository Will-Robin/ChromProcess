import os
from pathlib import Path
from ChromProcess import Classes
from ChromProcess import file_import
from ChromProcess import processing_functions as p_f
import matplotlib.pyplot as plt

exp_name = 'xylulose_1_2021'
MS_baseline_subtract = True
ms_thres = 500
Path_file = r'C:\Users\willi\Documents\Data\Data_information.csv'

'''
A pre-inspection process for data scraping.

Checks for the presence of conditions and analysis files in the source and
target storage folders. Their presence is either True or Flase as indicated
via an on-screen printout

Plots of the chromatograms are made to aid with selecting the correct regions
for peak integration.
'''


exp_paths = Classes.DataPaths(Path_file)
storage_stem = Path(r'C:\Users\willi\Documents\Data')

data_type = exp_paths.exp_code_path[exp_name].data_type
exp_path = exp_paths.exp_code_path[exp_name].path
# Point the program to where the data are
source_folder = Path(exp_path)
# State directory in which to store results
store_folder = Path(storage_stem/data_type/'FRN'/exp_name)

checks = {}
# try to find conditions file
source_cond_file = source_folder/'{}_conditions.csv'.format(exp_name)
if not os.path.exists(source_cond_file):
    checks['Source conditions file'] = False
else:
    checks['Source conditions file'] = True

local_cond_file = store_folder/'{}_conditions.csv'.format(exp_name)
if not os.path.exists(local_cond_file):
    checks['Local conditions file'] = False
else:
    checks['Local conditions file'] = True


source_analysis_det = source_folder/'{}_analysis_details.csv'.format(exp_name)
if not os.path.exists(source_analysis_det):
    checks['Source analysis file'] = False
else:
    checks['Source analysis file'] = True

local_analysis_det = store_folder/'{}_analysis_details.csv'.format(exp_name)
if not os.path.exists(local_analysis_det):
    checks['Local analysis file'] = False
else:
    checks['Local analysis file'] = True

# Read in the data files
if data_type == 'GCMS':
    # load .cdf files from GCMS analysis
    chroms, _ = file_import.load_cdf_from_directory(source_folder,
                                                    ms = MS_baseline_subtract)
for c in chroms:
    if data_type == 'GCMS' and MS_baseline_subtract:
        # remove low intensity ion chromatograms
        # reconstitue total ion chromatogram
        p_f.MS_intensity_threshold_chromatogram(c, threshold = ms_thres)

for c in checks:
    print(c,':',checks[c])

for c in chroms:
    plt.plot(c.time, c.signal)
plt.show()

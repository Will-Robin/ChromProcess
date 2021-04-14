import os
from pathlib import Path
from ChromProcess import file_import
from ChromProcess import Classes
from ChromProcess import info_params

exp_paths = Classes.DataPaths(r'C:\Users\willi\Documents\Data\Data_information.csv')
storage_stem = Path(r'C:\Users\willi\Documents\Data')

# path to calibration folder
calib_path = Path(r'C:\Users\willi\Documents\Data\Analysis_information')
# Get calibration files to use
alloc_file = calib_path/'Calibration_file_allocations.csv'
cal_alloc = file_import.importCalibrationFileAllocations(alloc_file)

path_list = [e for e in exp_paths.paths]

start_token = ''
end_token = ''
idx_s = 0
idx_e = len(path_list)
for c,p in enumerate(exp_paths.experiment_codes):
    if p == start_token:
        idx_s = c
    if p == end_token:
        idx_e = c + 1

path_list = path_list[idx_s:idx_e]

for e in path_list:
    print(e.experiment_code)

    exp_name = e.experiment_code
    method = e.data_type

    # point to data folder
    rt_fldr = r'C:\Users\willi\Documents\Data\{}\FRN\{}'.format(method,exp_name)
    root_data_folder = Path(rt_fldr)

    # point program to peak table
    data_folder = root_data_folder/'PeakTables'

    # read in peak table
    peak_tables = []
    for file in os.listdir(data_folder):
        peak_tables.append(Classes.PeakCollection(file = data_folder/file))

    for p in peak_tables:
        p.write_to_file(directory = data_folder)

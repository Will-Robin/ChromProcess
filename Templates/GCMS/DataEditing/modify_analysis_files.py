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
    exp_name = e.experiment_code
    method = e.data_type
    # Get analysis file
    rt_fldr = r'C:\Users\willi\Documents\Data\{}\FRN\{}'.format(method,exp_name)
    root_data_folder = Path(rt_fldr)
    analysis_path = root_data_folder/'{}_analysis_details.csv'.format(exp_name)
    analysis = Classes.Analysis_Information(information_file = analysis_path)

    # overwrite some information
    analysis.dilution_factor_error = 4.47E-02
    analysis.internal_ref_concentration_error = 9.89E-6

    # re-write to folder
    analysis.write_to_file(directory = root_data_folder)

from pathlib import Path
from ChromProcess import Classes

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
        idx_e = c

path_list = path_list[idx_s:idx_e]

for e in path_list:
    exp_name = e.experiment_code
    method = e.data_type

    # point to data folder
    rt_fldr = r'C:\Users\willi\Documents\Data\{}\FRN\{}'.format(method,exp_name)
    root_data_folder = Path(rt_fldr)

    file_path = root_data_folder/'DataReports'

    # Load in integral data report
    report = Classes.DataReport(file = file_path/'FRN073E_GCMS_integral_report.csv')

    # Remove duplicate entries, taking the entry with the higher values
    report.remove_repeat_entries()

    # remove low-concentration components below a bound
    report.remove_entries_below_threshold(0.1)

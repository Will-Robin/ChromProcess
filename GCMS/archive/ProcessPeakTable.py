import os
from pathlib import Path
from ChromProcess import file_import
from ChromProcess import Classes
from ChromProcess import info_params

# Choose experiment
exp_name = 'FRN107A'
method = 'GCMS'

# point to data folder
root_data_folder = Path(r'C:\Users\willi\Documents\Data\{}\FRN\{}'.format(method,exp_name))
# point program to peak table
data_folder = root_data_folder/'PeakTables'
# Point program to MS folder
MS_folder = root_data_folder/'PeakMassSpectra'
# path to calibration folder
calib_path = Path(r'C:\Users\willi\Documents\Data\Analysis_information')
# Get calibration files to use
alloc_file = calib_path/'Calibration_file_allocations.csv'
cal_alloc = file_import.importCalibrationFileAllocations(alloc_file)
if method == 'GCMS':
    Calib_file = cal_alloc.GCMS_allocations[exp_name]
elif method == 'HPLC':
    Calib_file = cal_alloc.HPLC_allocations[exp_name]
# read in modified boundaries
modified_bounds = file_import.read_local_assignments(root_data_folder/'{}_local_assignments.csv'.format(exp_name))

calib = Classes.Instrument_Calibration(file = calib_path/e.data_type/Calib_file)
calib.modify_boundaries(modified_bounds)

# Get conditions file
for file in os.listdir(root_data_folder):
    if 'conditions' in file:
        series_info = file_import.ImportExperimentConditionsFile(root_data_folder/file)

IS_pos = calib.internal_standard_position
IS_conc = series_info.internal_ref_concentration
dil = series_info.dilution

# read in peak table
peak_tables = []
for file in os.listdir(data_folder):
    peak_tables.append(file_import.read_peak_table(data_folder/file))
# Read in mass spectra
mass_spectra = []
for file in os.listdir(MS_folder):
    mass_spectra.append(file_import.read_mass_spectra_report(MS_folder/file))
# Combine mass spectra with peak_tables
for pt,ms in zip(peak_tables, mass_spectra):
    pt.mass_spectra = ms
# Create series of peak collections
series = Classes.PeakCollectionSeries(peak_tables, name = exp_name,
                                      conditions = series_info.conditions)

# align the peaks using the internal standard
series.align_peaks_to_IS(IS_pos)
# reference integrals to internal standard
series.reference_integrals_to_IS()
# cluster peaks
series.get_peak_clusters(bound = 0.025)
# assign peaks
series.assign_peaks(calib.boundaries)
# assign clusters
series.assign_clusters(calib.boundaries)
# Apply calibrations to assigned peaks
series.apply_calibrations(series_info, calib)
# Correct for sample dilution
series.apply_peak_dilution_factors(dil)
# Create arrays for the series data
series.make_integral_series()
series.make_concentration_series()

# Output data report
os.makedirs(root_data_folder/'DataReports', exist_ok = True)
dest_dir = root_data_folder/'DataReports'
series.write_concentrations_to_file(dest_dir/exp_name,calib)
series.write_integrals_to_file(dest_dir/exp_name,calib)

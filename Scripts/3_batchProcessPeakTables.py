import os
from pathlib import Path
from ChromProcess import file_import
from ChromProcess import Classes
from ChromProcess import info_params

start_token = '1A' # experiment code from which to begin
end_token = '5C' # Experiment code to end at

# peaks with integral values below peak_removal_limit will be ignored
peak_removal_limit = 0.1
# peak_agglomeration_boundary is the expected standard deviation for the
# peak position. This value is used to group peaks with similar retention times
# together (between chromatograms)
peak_agglomeration_boundary = 0.025

# importing information
system_root = Path('/Users/williamrobinson/Documents/Nijmegen')
storage_stem = system_root/'PrebioticDatabase'
data_folder = storage_stem/'Data/GCMS/FRN'
Path_file = storage_stem/'Data_information/Data_information.csv'
# path to calibration folder
calib_path = storage_stem/'Analysis_information'
# calibration allocations file
alloc_file = calib_path/'Calibration_file_allocations.csv'

exp_paths = Classes.DataPaths(Path_file)
cal_alloc = Classes.CalibrationAllocations(alloc_file)

path_list = [e for e in exp_paths.paths]

# build a list of experiments to process
idx_s = 0
idx_e = len(path_list)
for c,p in enumerate(exp_paths.experiment_codes):
    if p == start_token:
        idx_s = c
    if p == end_token:
        idx_e = c + 1

path_list = path_list[idx_s:idx_e]

# iterate through the peak tables and collate the data into data reports
for e in path_list:
    print(e.experiment_code)

    exp_name = e.experiment_code
    method = e.data_type

    # point to data folder
    rt_fldr = data_folder/f'{exp_name}'
    root_data_folder = Path(rt_fldr)

    # point program to peak table
    peak_table_folder = root_data_folder/'PeakTables'

    if method == 'GCMS':
        Calib_file = cal_alloc.GCMS_allocations[exp_name]
    elif method == 'HPLC':
        Calib_file = cal_alloc.HPLC_allocations[exp_name]

    # read in modified boundaries
    mod_bd_file = root_data_folder/f'{exp_name}_local_assignments.csv'
    modified_bounds = file_import.read_local_assignments(mod_bd_file)

    calib_file_path = calib_path/e.data_type/Calib_file
    calib = Classes.Instrument_Calibration(file = calib_file_path)
    calib.modify_boundaries(modified_bounds)

    # Get conditions file
    info_path = root_data_folder/'{}_conditions.csv'.format(exp_name)
    conditions = Classes.Experiment_Conditions(information_file = info_path)
    # Get analysis file
    analysis_path = root_data_folder/'{}_analysis_details.csv'.format(exp_name)
    analysis = Classes.Analysis_Information(information_file = analysis_path)

    # Internal standard position to align to
    IS_pos = calib.internal_standard_position
    # Internal standard concentration
    IS_conc = analysis.internal_ref_concentration
    # dilution factor
    dil = analysis.dilution_factor

    # read in peak table
    peak_tables = []
    for file in os.listdir(peak_table_folder):
        peak_tables.append(Classes.PeakCollection(file = peak_table_folder/file))

    # Read in mass spectra
    # Point program to MS folder
    MS_folder = root_data_folder/'PeakMassSpectra'
    mass_spectra = []
    for file in os.listdir(MS_folder):
        MS_file = MS_folder/file
        mass_spectra.append(file_import.read_mass_spectra_report(MS_file))
    # Combine mass spectra with peak_tables
    for pt,ms in zip(peak_tables, mass_spectra):
        pt.mass_spectra = ms

    # Create series of peak collections
    series = Classes.PeakCollectionSeries(peak_tables, name = exp_name,
                                          conditions = conditions.conditions)

    # align the peaks using the internal standard
    series.align_peaks_to_IS(IS_pos)
    # reference integrals to internal standard
    series.reference_integrals_to_IS()
    # remove peaks below a certain threshold
    series.remove_peaks_below_threshold(peak_removal_limit)
    # cluster peaks
    series.get_peak_clusters(bound = peak_agglomeration_boundary)
    # assign peaks
    series.assign_peaks(calib.boundaries)
    # assign clusters
    series.assign_clusters(calib.boundaries)
    # Apply calibrations to assigned peaks
    series.apply_calibrations(analysis, calib)
    # add in error analysis
    series.calculate_conc_errors(calib, analysis)
    # Correct for sample dilution
    series.apply_peak_dilution_factors(analysis)

    # Create arrays for the series data
    series.make_integral_series()
    series.make_concentration_series()

    # Output data report
    os.makedirs(root_data_folder/'DataReports', exist_ok = True)
    dest_dir = root_data_folder/'DataReports'

    conc_data_report = series.create_conc_DataReport(calib)
    integral_data_report = series.create_integral_DataReport(calib)

    conc_data_report.write_to_file(path = dest_dir)
    integral_data_report.write_to_file(path = dest_dir)

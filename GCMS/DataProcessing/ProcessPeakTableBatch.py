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

start_token = 'FRN108'
end_token = 'FRN126'
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

    if method == 'GCMS':
        Calib_file = cal_alloc.GCMS_allocations[exp_name]
    elif method == 'HPLC':
        Calib_file = cal_alloc.HPLC_allocations[exp_name]

    # read in modified boundaries
    mod_bd_file = root_data_folder/'{}_local_assignments.csv'.format(exp_name)
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
    for file in os.listdir(data_folder):
        peak_tables.append(Classes.PeakCollection(file = data_folder/file))

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
    series.remove_peaks_below_threshold(0.0)
    # cluster peaks
    series.get_peak_clusters(bound = 0.025)
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

    import matplotlib.pyplot as plt
    from ChromProcess import info_params
    fig, ax = plt.subplots(ncols = 2)
    x_ax = integral_data_report.series_values
    for de in conc_data_report.data:
        nm = info_params.smiles_to_names[de.split('/')[0]]
        y_ax = conc_data_report.data[de]
        y_error = conc_data_report.errors[de]
        ax[1].plot(x_ax, y_ax, c = info_params.colour_assignments[de.split('/')[0]], label = nm)
        ax[1].fill_between(x_ax, y_ax+y_error, y2 = conc_data_report.data[de]-y_error,
                        color = info_params.colour_assignments[ de.split('/')[0]], alpha = 0.5)
        ax[1].errorbar(x_ax,y_ax, yerr = y_error, color = info_params.colour_assignments[ de.split('/')[0]], alpha = 0.5)

    for de in integral_data_report.data:
        name = de.split(' ')[0]
        if name in info_params.colour_assignments:
            ax[0].plot(series.series_values, integral_data_report.data[de],
                        c = info_params.colour_assignments[name])
        else:
            ax[0].plot(series.series_values, integral_data_report.data[de],
                        c = 'k', alpha = 0.5)
    plt.show()
    plt.close()

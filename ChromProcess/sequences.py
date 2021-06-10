from ChromProcess import Classes
from ChromProcess import file_import


def peak_table_process_1(experiment, root_data_folder,
                        Calib_file,
                        IS_align = True, integral_cutoff = 0.1,
                        peak_bin_threshold = 0.025):

    data_folder = root_data_folder/'PeakTables'

    exp_name = experiment.experiment_code
    method = experiment.data_type

    if method == 'GCMS':
        Calib_file = cal_alloc.GCMS_allocations[exp_name]
    elif method == 'HPLC':
        Calib_file = cal_alloc.HPLC_allocations[exp_name]

    # read in modified boundaries
    mod_bd_file = root_data_folder/'{}_local_assignments.csv'.format(exp_name)
    modified_bounds = file_import.read_local_assignments(mod_bd_file)


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
    if IS_align:
        series.align_peaks_to_IS(IS_pos)
    # reference integrals to internal standard
    series.reference_integrals_to_IS()
    # remove peaks below a certain threshold
    series.remove_peaks_below_threshold(integral_cutoff)
    # cluster peaks
    series.get_peak_clusters(bound = peak_bin_threshold)
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

    return series

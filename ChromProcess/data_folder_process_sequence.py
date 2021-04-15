def chrom_folder_process_sequence(source_folder, store_folder,
                                  conditions_file, analysis_file,
                                  copy_analysis = False,
                                  copy_conditions = False):
    '''
    Parameters
    ----------
    source_folder: str or pathlib Path object
        path to source of data
    store_folder: str or pathlib Path object
        path to extract data into
    conditions_file: str or pathlib Path object
        Path to conditions file
    analysis_file: str or pathlib Path object
        Path to analysis file

    Returns
    -------
    bool
        whether daat exctraction ran or not.
    '''
    import os
    from ChromProcess import Classes
    from ChromProcess import file_import
    from ChromProcess import file_output
    from ChromProcess import peak_operations as peak_ops
    from ChromProcess import processing_functions as p_f
    from ChromProcess import chromatogram_operations as chrom_ops

    # Create the store folder if it does not already exist
    os.makedirs(store_folder, exist_ok = True)

    # Get experiment conditions information
    if os.path.exists(conditions_file):
        conditions = Classes.Experiment_Conditions(information_file = conditions_file)
    else:
        print('Conditions file not found or parsing issues.')
        print('Passing data set {}.'.format(e.experiment_code))
        return False

    # Get experiment analysis details
    if os.path.exists(analysis_file):
        analysis = Classes.Analysis_Information(information_file = analysis_file)
    else:
        print('Analysis file not found or parsing issues.')
        print('Passing data set {}.'.format(e.experiment_code))
        return False

    # Read in the data files
    if analysis.analysis_type == 'GCMS':
        # load .cdf files from GCMS analysis
        chroms, _ = file_import.load_cdf_from_directory(source_folder,
                                                        ms = analysis.use_MS)
    elif analysis.analysis_type == 'HPLC':
        # load .txt files exported from Shimadzu LabSolutions for HPLC analysis
        chroms, _ = file_import.directoryLoadShimadzuASCII(source_folder)
    else:
        print('analysis_type provided is {}.')
        print('Please choose from GCMS or HPLC. Passing data set.')
        return False

    # Pre-process the chromatograms
    for c in chroms:
        if analysis.analysis_type == 'GCMS':
            # remove low intensity ion chromatograms and reconstitue
            # total ion chromatogram
            p_f.MS_intensity_threshold_chromatogram(c,
                                                 threshold = analysis.MS_cutoff)
    # Find integral information
    for c in chroms:
        if analysis.analysis_type == "GCMS":
            # Get interal reference integrals
            chrom_ops.internalRefIntegral(c, analysis.internal_ref_region)

        # Get peaks in regions of the chromatogram
        for r in analysis.regions:
            chrom_ops.pickPeaksRegion(c, r,
                                       threshold = analysis.peak_pick_threshold)

    # Output chromatograms to store folder
    os.makedirs(store_folder/'Chromatograms', exist_ok = True)
    dest_dir = store_folder/'Chromatograms'
    for c in chroms:
        file_output.chromatogram_to_csv_GCMS(c, filename = dest_dir/c.filename)

    # Output peak table
    os.makedirs(store_folder/'PeakTables', exist_ok = True)
    dest_dir = store_folder/'PeakTables'
    for c,v in zip(chroms, conditions.series_values):
        file_output.write_peak_table(c, filename = dest_dir/c.filename,
                                     value = v,
                                     series_unit = conditions.series_unit)
    # Output peak chromatograms
    if analysis.analysis_type == 'GCMS':
        # Output peak mass spectra
        os.makedirs(store_folder/'PeakMassSpectra', exist_ok = True)
        dest_dir = store_folder/'PeakMassSpectra'
        for c in chroms:
            for p in c.peaks:
                peak_ops.peakIonChromatogram(c.peaks[p],c)

            file_output.write_peak_mass_spectra(c, filename = dest_dir/c.filename)

    # Copy analysis and conditions information into target folder if
    # required
    if copy_analysis:
        analysis.write_to_file(directory = store_folder)

    if copy_conditions:
        conditions.experiment_code = e.experiment_code
        conditions.write_to_file(directory = store_folder)


    return True
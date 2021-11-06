def chrom_folder_process_sequence(source_folder, store_folder,
                                  conditions_file, analysis_file,
                                  copy_analysis = False,
                                  copy_conditions = False):
    '''
    Parameters
    ----------
    source_folder: str or pathlib Path object
        path to folder containing the source data

    store_folder: str or pathlib Path object
        path folder to extract data into.

    conditions_file: str or pathlib Path object
        Path to conditions file

    analysis_file: str or pathlib Path object
        Path to analysis file

    copy_analysis: bool
        Whether to copy the analysis file from the source folder to the
        targer folder or not.

    copy_conditions: bool
        Whether to copy the conditions file from the source folder to the
        targer folder or not.

    Returns
    -------
    False: if scrape fails
    chroms: list of ChromProcess Chromatogram objects
    '''
    import os
    from pathlib import Path
    from ChromProcess import Classes
    from ChromProcess import file_import
    from ChromProcess import file_output
    from ChromProcess import peak_operations as peak_ops
    from ChromProcess import processing_functions as p_f
    from ChromProcess import chromatogram_operations as chrom_ops

    if isinstance(source_folder, str):
        source_folder = Path(source_folder)
    if isinstance(store_folder,str):
        store_folder = Path(store_folder)
    if isinstance(conditions_file,str):
        conditions_file = Path(conditions_file)
    if isinstance(analysis_file,str):
        analysis_file = Path(analysis_file)

    # Create the store folder if it does not already exist
    os.makedirs(store_folder, exist_ok = True)

    # Get experiment conditions information
    if os.path.exists(conditions_file):
        conditions = Classes.Experiment_Conditions(information_file = conditions_file)
    else:
        print('Conditions file not found.')
        return False

    # Get experiment analysis details
    if os.path.exists(analysis_file):
        analysis = Classes.Analysis_Information(information_file = analysis_file)
    else:
        print('Analysis file not found.')
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
        try: 
            chroms = []
            for f in os.listdir(source_folder):
                chroms.append(Classes.Chromatogram(source_folder/f))
        except:
            print('Incorrect chromatogram type in analysis file.')
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
            chrom_ops.regionPeakPick(c, r, threshold = analysis.peak_pick_threshold)
            chrom_ops.integratePeaks(c)

    # Output chromatograms to store folder
    os.makedirs(store_folder/'Chromatograms', exist_ok = True)
    dest_dir = store_folder/'Chromatograms'
    for c in chroms:
        c.write_to_csv(filename = dest_dir/c.filename)

    # Output peak table
    os.makedirs(store_folder/'PeakTables', exist_ok = True)
    dest_dir = store_folder/'PeakTables'
    for c,v in zip(chroms, conditions.series_values):
        c.write_peak_table(filename = dest_dir/c.filename,
                            value = v,
                            series_unit = conditions.series_unit)
    # Output peak chromatograms
    if analysis.analysis_type == 'GCMS':
        # Output peak mass spectra
        os.makedirs(store_folder/'PeakMassSpectra', exist_ok = True)
        dest_dir = store_folder/'PeakMassSpectra'
        for c in chroms:
            for p in c.peaks:
                peak_ops.peakMassSpectrum(c.peaks[p],c)

            c.write_peak_mass_spectra(filename = dest_dir/c.filename)

    # Copy analysis and conditions information into target folder if
    # required
    if copy_analysis:
        analysis.write_to_file(directory = store_folder)

    if copy_conditions:
        conditions.experiment_code = analysis.experiment_code
        conditions.write_to_file(directory = store_folder)


    return chroms

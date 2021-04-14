import os
from ChromProcess import plotting
from ChromProcess import file_output
from ChromProcess import mass_spectra
from ChromProcess import series_operations


def series_from_chromatogram_set(series, information, threshold = 0.1,
                                 max_intensity = 1e100, min_intensity = -1e100,
                                 combine_assignments = False,
                                 plot_mass_spectra = False,
                                 label_assignments = True,
                                 cluster_dev = 0.01):

    '''
    Set of functions for integrating peaks form chromatograms and plotting results.
    Parameters
    ----------
    series: ChromProcess Series Object
        A series of chromatograms to be operated upon including experimental
        information.
    information: ChromProcess Calibration_File object
        Calibration factors and assignment boundaries for chromatograms.
    threshold: float
        Peaks below this fraction of the maximum of a signal in a given region
        will not be detected.
    max_intensity: float
        Peaks above this threshold will not be detected.
    min_intensity: float
        Peaks below this threshold will not be detected.
    combine_assignments: bool
        whether to call combine_assignments function or not.

    Returns
    -------
    None
    '''

    '''Input stage'''
    if series.chromatograms[0].c_type != information.type: # Checking the information provided versus chromatogram type.
        print("Warning: Data file type is considered {}, whilst the calibration file type is considered as {}.".format(series.chromatograms[0].c_type,information.type))
        print("Proceeding with analysis anyway.")

    bounds = information.boundaries # Get assignment boundaries
    calibrations = information.calibrations # Get calibration factors.

    '''Processing stage'''
    # GCMS data are usually recorded using an internal standard
    # so a function is called to get its parameters.
    if series.chromatograms[0].c_type == "GCMS":
        series_operations.get_internal_ref_integrals(series)

    # Find peaks in each chromatogram, region by region
    series_operations.pick_peaks(series, threshold = threshold, max_intensity = max_intensity, min_intensity = min_intensity)

    series_operations.get_integrals(series) # Adds intergral information into the picked peaks

    if series.chromatograms[0].c_type == "GCMS":
        for c in series.chromatograms: # Extracting ion chromatograms into Peak objects and integrating them.
            mass_spectra.peak_ion_chromatograms(c)
            mass_spectra.integrate_ion_chromatograms(c, threshold = 0.03)

    series_operations.generate_peak_series(series, stdev = cluster_dev) # Makes sequences of peaks

    if series.chromatograms[0].c_type == "GCMS":
        mass_spectra.get_peak_mass_spectra(series) # adds mass spectra into Peak objects

    series_operations.apply_calibration_to_peak_series(series, calibrations,bounds) # Converts TIC integrals to concentrations if calibrations exist.

    if combine_assignments: # Some compounds have two peaks in chromatograms. They may be combined numerically.
        series_operations.combine_common_assignments(series,bounds)

    if series.chromatograms[0].c_type == "GCMS":
        mass_spectra.ion_chromatogram_integral_series(series) # Makes sequences of ion chromatograph integrals.

    '''Output stage'''
    os.makedirs("Timecourses", exist_ok = True)
    os.chdir("Timecourses")

    plotting.plot_conc_series(series,bounds) # plots concentrations over series

    plotting.plot_integral_series(series,bounds) # plots integrals over series.

    os.chdir("..")

    os.makedirs("Data_reports", exist_ok = True)
    os.chdir("Data_reports")

    file_output.data_report_template_convert(series, information) # Creates a data report with concentration data

    file_output.report_all_peak_integrals(series, information) # Creates a data report with integral data
    os.chdir("..")

    '''Plot each chromatogram in each region of the series in a new folder.'''
    os.makedirs("Chromatograms", exist_ok = True)
    os.chdir("Chromatograms")
    for c in series.chromatograms:
        file_output.chromatogram_to_csv_GCMS(c) # Saves chromatograms as .csv files.

    plotting.plot_chromatograms(series, information) # Plots overlayed chromatograms with picked peaks + peak bounds, region by region.

    os.chdir("..")

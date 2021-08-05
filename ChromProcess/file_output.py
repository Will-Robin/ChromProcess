from ChromProcess import file_import, processing_functions
import numpy as np
from ChromProcess import info_params as i_p

'''
Functions for outputting files.
'''

def chromatogram_to_csv_HPLC(chromatogram):

    chrom_out = np.vstack((chromatogram.time,chromatogram.signal))
    chrom_out = chrom_out.T
    np.savetxt('{}_chromatogram.csv'.format(chromatogram.filename),chrom_out,delimiter = ',',header = 'time, signal')

def chromatogram_to_csv_GCMS(chromatogram, filename = 'chromatogram'):
    '''
    Converts a GCMS .cdf file to .csv chromatogram.
    Parameters
    ----------
    f: str
        file name of the GCMS .cdf file
    outname: str
        name for the output file

    '''
    time   = chromatogram.time
    signal = chromatogram.signal
    chrom_out = np.vstack((time,signal))
    chrom_out = chrom_out.T

    with open('{}.csv'.format(filename), 'w') as f:
        f.write('time/ min, signal/ total ion counts\n')
        for x in range(0,len(chrom_out)):
            f.write('{},{}\n'.format(chrom_out[x,0],chrom_out[x,1]))

def data_report_template_convert(series, information):
    '''
    Output data in report template.
    Parameters
    ----------

    Returns
    -------

    '''
    name = "{}_{}_{}_data_report.csv".format(series.set_name,"conc_series",series.chromatograms[0].c_type)
    with open(name, "w") as outfile:
        outfile.write("Dataset,{}\n".format(series.set_name))

        outfile.write("start_conditions\n")
        for c in series.conditions:
            outfile.write("{},".format(c))
            [outfile.write("{},".format(x)) for x in series.conditions[c]]
            outfile.write("\n")
        outfile.write("end_conditions\n")

        outfile.write("start_analysis_details\n")

        outfile.write("Calibrations_bounds_date,{}\n".format(information.date))
        outfile.write("Chromatography_method,{},{}\n".format(information.type, information.method))
        outfile.write("Derivatisation_method,{}\n".format(information.derivatisation))

        outfile.write("end_analysis_details\n")

        outfile.write("start_data\n")

        out = [series.x_series]
        p_header = [series.x_name]

        for p in series.conc_series:
            label_name = processing_functions.name_peak(float(p),information.boundaries)
            if processing_functions.int_test(label_name[-1]):
                label_name = label_name[:-2]

            species_smiles = i_p.canonical_SMILES[label_name]
            p_header.append(species_smiles+" M")
            out.append(series.conc_series[p])

        out = [list(i) for i in zip(*out)]

        [outfile.write("{},".format(x)) for x in p_header]
        outfile.write("\n")
        for o in out:
            [outfile.write("{},".format(x)) for x in o]
            outfile.write("\n")

        outfile.write("end_data\n")

def peak_ion_chromatograms_report(series, information):
    '''
    Output data in report template.
    Parameters
    ----------

    Returns
    -------

    '''

    for s in series.ion_series:

        name = "{}_{}_peak_report.csv".format(series.set_name,s)

        with open(name, "w") as outfile:
            outfile.write("Dataset,{}\n".format(series.set_name))
            outfile.write("Peak,{}\n".format(s))
            outfile.write("start_conditions\n")
            for c in series.conditions:
                outfile.write("{},".format(c))
                [outfile.write("{},".format(x)) for x in series.conditions[c]]
                outfile.write("\n")
            outfile.write("end_conditions\n")

            outfile.write("start_analysis_details\n")

            outfile.write("Calibrations_bounds_date,{}\n".format(information.date))
            outfile.write("Chromatography_method,{},{}\n".format(information.type, information.method))
            outfile.write("Derivatisation_method,{}\n".format(information.derivatisation))

            outfile.write("end_analysis_details\n")

            outfile.write("start_data\n")

            # create ion integral traces
            p_header = [series.x_name]
            out = [series.x_series]

            for m in series.ion_series[s]:
                p_header.append(m)
                out.append(series.ion_series[s][m])

            out = [list(i) for i in zip(*out)]

            [outfile.write("{},".format(x)) for x in p_header]
            outfile.write("\n")
            for o in out:
                [outfile.write("{},".format(x)) for x in o]
                outfile.write("\n")

            outfile.write("end_data\n")

def report_all_peak_integrals(series, information):

    '''
    Output data in report template.
    Parameters
    ----------

    Returns
    -------

    '''

    y_data = series.integral_series
    x_data = np.array([x.timepoint for x in series.chromatograms])
    out_unit = " integral"
    name = '{}'.format("{}_{}_all_integrals_report.csv".format(series.set_name,series.chromatograms[0].c_type))

    with open(name, "w") as outfile:
        outfile.write("Dataset,{}\n".format(series.set_name))

        outfile.write("start_conditions\n")
        outfile.write("dilution_factor, {}\n".format(series.dilution))
        outfile.write("internal_ref_concentration/ M, {}\n".format(series.internal_ref_concentration))
        for c in series.conditions:
            outfile.write("{},".format(c))
            [outfile.write("{},".format(x)) for x in series.conditions[c]]
            outfile.write("\n")
        outfile.write("end_conditions\n")

        outfile.write("start_analysis_details\n")

        outfile.write("Calibrations_bounds_date,{}\n".format(information.date))
        outfile.write("Chromatography_method,{},{}\n".format(information.type, information.method))
        outfile.write("Derivatisation_method,{}\n".format(information.derivatisation))

        outfile.write("end_analysis_details\n")

        outfile.write("start_data\n")

        out = x_data
        p_header = [series.x_name]

        for p in series.integral_series:
            label_name = processing_functions.name_peak(float(p),information.boundaries)
            out = np.vstack((out,series.integral_series[p]))
            p_header.append("{}({})  {}".format(label_name,p,out_unit))

        out = [list(i) for i in zip(*out)]

        [outfile.write("{},".format(x)) for x in p_header]
        outfile.write("\n")
        for o in out:
            [outfile.write("{},".format(x)) for x in o]
            outfile.write("\n")

        outfile.write("end_data\n")


def report_integral_log(series, information):
    '''
    Output integral log report.
    Parameters
    ----------

    Returns
    -------

    '''

    y_data = series.integral_series
    x_data = np.array([x.timepoint for x in series.chromatograms])
    out_unit = " integral"
    name = '{}'.format("{}_internal_reference_report.csv".format(series.set_name))

    with open(name, "w") as outfile:
        outfile.write("Dataset,{}\n".format(series.set_name))
        outfile.write("Start_Date,{}\n".format(i_p.dates[series.set_name]))
        out = x_data
        p_header = [series.x_name]

        for c,r in enumerate(y_data):
            for p in y_data[r]:
                label_name = processing_functions.name_peak(float(p),information.boundaries)
                out = np.vstack((out,y_data[r][p]))
                p_header.append("{}({})  {}".format(label_name,p,out_unit))
        out = [list(i) for i in zip(*out)]

        [outfile.write("{},".format(x)) for x in p_header]
        outfile.write("\n")
        for c,o in enumerate(out):
            [outfile.write("{},".format(x)) for x in o]
            outfile.write("\n")

    return 0

def write_deconvoluted_output(series, information):
    '''
    Output data in report template.
    Parameters
    ----------

    Returns
    -------

    '''
    name = "{}_deconvolution_report.csv".format(series.set_name)
    with open(name, "w") as outfile:
        outfile.write("Dataset,{}\n".format(series.set_name))

        outfile.write("start_conditions\n")
        for c in series.conditions:
            outfile.write("{},".format(c))
            [outfile.write("{},".format(x)) for x in series.conditions[c]]
            outfile.write("\n")
        outfile.write("end_conditions\n")

        outfile.write("start_analysis_details\n")

        outfile.write("Calibrations_bounds_date,{}\n".format(information.date))
        outfile.write("Chromatography_method,{},{}\n".format(information.type, information.method))
        outfile.write("Derivatisation_method,{}\n".format(information.derivatisation))

        outfile.write("end_analysis_details\n")

        outfile.write("start_data\n")

        out = [series.x_series]
        p_header = [series.x_name]

        for p in series.deconvoluted_series:
            print(p)
            p_header.append("{} integral".format(p))
            out.append(series.deconvoluted_series[p])

        out = [list(i) for i in zip(*out)]

        [outfile.write("{},".format(x)) for x in p_header]
        outfile.write("\n")
        for o in out:
            [outfile.write("{},".format(x)) for x in o]
            outfile.write("\n")

        outfile.write("end_data\n")

def store_chromatogram(chrom):
    '''
    Parameters
    ----------
    chrom: ChromProcess Chromatogram object`
        Chromatogram to be stored.
    Returns
    --------
        None
    '''
    import h5py

    with h5py.File("{}.hdf5".format(chrom.filename), 'w') as f:
        f.create_dataset('time_min', data=chrom.time)
        f.create_dataset('total_intensity_counts', data=chrom.signal)
        f.create_dataset("type", data=chrom.c_type)
        if chrom.mass_spectra:
            f.create_dataset("intensity_values_counts", data=chrom.mass_intensity)
            f.create_dataset("mass_values_mz", data=chrom.mass_values)
            f.create_dataset("scan_index", data=chrom.scan_indices)
            f.create_dataset("point_count", data=chrom.point_counts)
        else:
            f.create_dataset("intensity_values_counts", data=[])
            f.create_dataset("mass_values_mz", data=[])
            f.create_dataset("scan_index", data=[])
            f.create_dataset("point_count", data=[])

def write_modified_calibration_params(modified_bounds, series):
    with open("{}_local_assignments.csv".format(series.set_name), "w") as f:
        for b in modified_bounds:
            f.write("{},".format(b))
            f.write("{}\n".format(",".join([str(x) for x in modified_bounds[b]])))

def write_peak_table(chromatogram, filename = 'Peak_table', value = 0.0, series_unit = "None"):
    '''
    For writing peak integrals from a chromatogram to a .csv file.
    Parameters
    ----------
    chromatogram: ChromProcess Chromatogram object
        Parent chromatogram for peaks
    filename: str
        Name for the file
    value: float, str, bool
        Value assigned to the variable assigned to the chromatogram.
    series_unit: float, str, bool
        Units of the variable assigned to the chromatogram.
    Returns
    -------
    None
    '''
    with open("{}.csv".format(filename), "w") as f:
        f.write("{},{}\n".format(series_unit,value))
        f.write('IS_retention_time/ min,IS_integral,IS_peak start/ min,IS_peak end/ min\n')
        if chromatogram.internal_reference:
            st_ind = chromatogram.internal_reference.indices[0]
            end_ind = chromatogram.internal_reference.indices[-1]
            IS_lower_bound = chromatogram.time[st_ind]
            IS_upper_bound = chromatogram.time[end_ind]
            IS_RT = chromatogram.internal_reference.retention_time
            IS_integral = chromatogram.internal_reference.integral
        else:
            IS_RT, IS_integral, IS_lower_bound, IS_upper_bound = 'None', 'None', 'None', 'None'

        f.write("{},{},{},{}\n".format(IS_RT, IS_integral, IS_lower_bound, IS_upper_bound))
        f.write("Retention_time/ min,integral,peak start/ min,peak end/ min\n")

        for p in chromatogram.peaks:
            st_ind = chromatogram.peaks[p].indices[0]
            end_ind = chromatogram.peaks[p].indices[-1]
            f.write("{},{},{},{}\n".format(chromatogram.peaks[p].retention_time, chromatogram.peaks[p].integral, chromatogram.time[st_ind], chromatogram.time[end_ind]))

def write_peak_mass_spectra(chromatogram, filename = ''):
    '''
    For writing mass spectra of chromatogram peaks to a .csv file.

    Parameters
    ----------
    chromatogram: ChromProcess Chromatogram object
        Parent chromatogram for peaks
    filename: str
        Name for the file

    Returns
    -------
    None
    '''
    with open('{}.csv'.format(filename), 'w') as f:

        for p in chromatogram.peaks:
            if chromatogram.peaks[p].mass_spectrum:
                ms = chromatogram.peaks[p].mass_spectrum

                f.write('Peak retention time, {}\n'.format(chromatogram.peaks[p].retention_time))
                f.write('m/z,')
                [f.write('{},'.format(x)) for x in ms[0]]
                f.write('\n')
                f.write('relative abundance,')
                [f.write('{},'.format(x/max(ms[1]))) for x in ms[1]]
                f.write('\n')
            else:
                f.write('Peak retention time, {}\n'.format(chromatogram.peaks[p].retention_time))
                f.write('m/z,empty\n')
                f.write('relative abundance,empty\n')

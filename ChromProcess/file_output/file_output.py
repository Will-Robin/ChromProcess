from ChromProcess import file_import, processing_functions
import numpy as np
from ChromProcess import info_params as i_p

'''
Functions for outputting files.
'''


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

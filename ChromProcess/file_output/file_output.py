from ChromProcess import file_import, processing_functions
import numpy as np
from ChromProcess import info_params as i_p

'''
Functions for outputting files.
'''

def store_chromatogram(chrom):
    '''
    For storing chromatograms in HDF5 format.

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

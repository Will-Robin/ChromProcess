import os
from pathlib import Path
import numpy as np
from ChromProcess import Classes
from ChromProcess import simple_functions as s_f
from ChromProcess import file_import
from NorthNet import info_params

'''
Compiling the mass spectra and retention time assignments for all data sets
into one file.
'''

# database root folder
database_root = Path(r'C:\Users\willi\Documents\PrebioticDatabase')
calib_path = database_root/'Analysis_information'
storage_stem = database_root/'Data'

# Get paths to calibration data
exp_paths = Classes.DataPaths(database_root/'Data_Information'/'Data_information.csv')
alloc_file = calib_path/'Calibration_file_allocations.csv'
cal_alloc = file_import.importCalibrationFileAllocations(alloc_file)

# Choosing experiments to analyse.
exp_list = [*exp_paths.exp_code_path][25:94] # 25 to 94
'''
for e in exp_paths.exp_code_path:
    if '2020_calibrations' in str(exp_paths.exp_code_path[e].path):
        exp_list.append(e)
'''
# To keep track of assigned mass spectra
assigns = {}
# to keep track of all observed masses
all_masses = np.array([])
# counter for number of peaks
peak_no = 0
# Import calibration data
for el in exp_list:

    e = exp_paths.exp_code_path[el]
    exp_name = e.experiment_code
    print('Loading {}.'.format(exp_name))

    exp_folder = storage_stem/e.data_type/'FRN'/exp_name
    MS_folder = exp_folder/'PeakMassSpectra'
    PT_folder = exp_folder/'PeakTables'

    # read in peak table
    peak_tables = []
    for file in os.listdir(PT_folder):
        peak_tables.append(Classes.PeakCollection(file = PT_folder/file))

    # Read in mass spectra
    # Point program to MS folder
    MS_folder = exp_folder/'PeakMassSpectra'
    mass_spectra = []
    for file in os.listdir(MS_folder):
        MS_file = MS_folder/file
        mass_spectra.append(file_import.read_mass_spectra_report(MS_file))

    # Combine mass spectra with peak_tables
    for pt,ms in zip(peak_tables, mass_spectra):
        pt.add_mass_spectra(ms)

    # Get conditions file
    info_path = exp_folder/'{}_conditions.csv'.format(exp_name)
    conditions = Classes.Experiment_Conditions(information_file = info_path)

    # Create series of peak collections
    series = Classes.PeakCollectionSeries(peak_tables, name = exp_name,
                                          conditions = conditions.conditions)

    Calib_file = cal_alloc.GCMS_allocations[exp_name]
    calib_file_path = calib_path/e.data_type/Calib_file
    calib = Classes.Instrument_Calibration(file = calib_file_path)

    # read in modified boundaries
    mod_bd_file = exp_folder/'{}_local_assignments.csv'.format(exp_name)
    modified_bounds = file_import.read_local_assignments(mod_bd_file)
    calib.modify_boundaries(modified_bounds)

    IS_pos = calib.internal_standard_position

    # align the peaks using the internal standard
    series.align_peaks_to_IS(IS_pos)
    # remove peaks below a certain threshold
    series.remove_peaks_below_threshold(0.05)
    # assign peaks
    series.assign_peaks(calib.boundaries)

    for pc in series.peak_collections:
        for pk in pc.peaks:
            nom = pk.assignment.split(' ')[0]
            if not s_f.isfloat(nom) and nom in info_params.canonical_SMILES:
                if pk.mass_spectrum:
                    asgn = info_params.canonical_SMILES[nom]
                    if asgn in assigns:
                        assigns[asgn].append(pk)
                    else:
                        assigns[asgn] = [pk]

                    peak_no += 1

                    all_masses = np.hstack((all_masses, pk.mass_spectrum.mz))

print('Arranging data.')
# Arranging data
all_masses = np.unique(all_masses)
all_masses = np.sort(all_masses)
data_mat = np.zeros((peak_no, len(all_masses) + 1))
integral_mat = np.zeros(peak_no)

compound_labels = ['']*peak_no
#print(len(compound_labels), data_mat.shape, all_masses.shape, integral_mat.shape)
# keeping track of peak indices
pk_idx = 0
for c,a in enumerate(assigns):
    for c2,m in enumerate(assigns[a]):
        data_mat[pk_idx,0] = m.retention_time
        for mz,ra in zip(m.mass_spectrum.mz, m.mass_spectrum.relative_abundances):
            idx = np.where(all_masses == mz)[0]
            data_mat[pk_idx,idx+1] = ra

        integral_mat[pk_idx] = m.integral

        compound_labels[pk_idx] = a

        pk_idx += 1

print('Writing data to files.')
import pickle
collated_folder = database_root/'AggregatedData'
with open(collated_folder/'All_mass_spectra', 'wb') as f:
    pickle.dump(data_mat,f)

with open(collated_folder/'mass_spectra_labels', 'wb') as f:
    pickle.dump(compound_labels,f)

with open(collated_folder/'all_masses', 'wb') as f:
    pickle.dump(all_masses,f)

with open(collated_folder/'integral_mat', 'wb') as f:
    pickle.dump(integral_mat,f)

print('complete')

import os
from pathlib import Path
from ChromProcess import Classes
from ChromProcess import data_folder_process_sequence as dfps

'''
A program for loading sets of chromatograms from the source_folder and
extracting basic information from them (chromatogram traces, peak integrals,
mass spectra), which are deposited in the store_folder.

Minimal processing is performed on the chromatograms.

It is recommended to run a pre-scrape inspection on each source folder
before running this program. The inspection will indicate whether the
appropriate files are present in the source and target directories, as well as
aiding in inspecting chromatograms and selecting chromatogram regions for
integration.
'''
# Choose the experiment code (see Data_information.csv)
exp_name = 'allose_2_2021'
# whether local analysis files and condition files are to be used
local_analysis = True
local_conditions = True
analysis_cpy = False
conditions_cpy = False

# importing information
Path_file = r'C:\Users\willi\Documents\Data\Data_information.csv'
exp_paths = Classes.DataPaths(Path_file)
storage_stem = Path(r'C:\Users\willi\Documents\Data')

experiment = exp_paths.exp_code_path[exp_name]
data_type = experiment.data_type
exp_path = experiment.path
# Point the program to where the data are
source_folder = Path(exp_path)
# State directory in which to store results
store_folder = Path(storage_stem/data_type/'FRN'/exp_name)

# Create file names
if local_conditions:
    cond_file = store_folder/'{}_conditions.csv'.format(exp_name)
else:
    cond_file = source_folder/'{}_conditions.csv'.format(exp_name)

if local_analysis:
    analysis_details=store_folder/'{}_analysis_details.csv'.format(exp_name)
else:
    analysis_details=source_folder/'{}_analysis_details.csv'.format(exp_name)

# Run data scraping
dfps.chrom_folder_process_sequence(source_folder, store_folder,
                                   cond_file, analysis_details,
                                   experiment,
                                   copy_analysis = analysis_cpy,
                                   copy_conditions = conditions_cpy)

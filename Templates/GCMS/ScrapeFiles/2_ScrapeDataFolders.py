import os
from pathlib import Path
from ChromProcess import Classes
from ChromProcess import data_folder_process_sequence as dfps

'''
A program for loading sets of chromatograms from multiple sources and
extracting basic information from them (chromatogram traces, peak integrals,
mass spectra), which are deposited in storage folders.
Minimal processing is performed on the chromatograms.

It is recommended to run a pre-scrape inspection on each source folder
before running this program. The inspection will indicate whether the
appropriate files are present in the source and target directories, as well as
aiding in inspecting chromatograms and selecting chromatogram regions for
integration.
'''
#Choose where to start and where to end (see Data_information.csv)
start_token = ''
end_token = ''
# Choose whether to use conditions and analysis files from the source folder
# (False) or from the target storage folder (True)
local_analysis = True
local_conditions = True
# Choose whether to copy the loaded conditions and analysis files or not
analysis_cpy = False
conditions_cpy = False

# Get paths to source folders
Path_file = r'C:\Users\willi\Documents\Data\Data_information.csv'
exp_paths = Classes.DataPaths(Path_file)
storage_stem = Path(r'C:\Users\willi\Documents\Data')
path_list = [e for e in exp_paths.paths]

# Build a list of paths from start to end
idx_s = 0
idx_e = len(path_list)
for c,p in enumerate(exp_paths.experiment_codes):
    if p == start_token:
        idx_s = c
    if p == end_token:
        idx_e = c + 1

path_list = path_list[idx_s:idx_e]

pass_list = []
for e in path_list:
    exp_name = e.experiment_code
    # Point the program to where the data are
    source_folder = Path(e.path)
    # State directory in which to store results
    store_folder = Path(storage_stem/e.data_type/'FRN'/exp_name)

    # Create path to conditions file
    if local_conditions:
        cond_file = store_folder/'{}_conditions.csv'.format(exp_name)
    else:
        cond_file = source_folder/'{}_conditions.csv'.format(exp_name)

    # Create path to conditions file
    if local_analysis:
        analysis_details=store_folder/'{}_analysis_details.csv'.format(exp_name)
    else:
        analysis_details=source_folder/'{}_analysis_details.csv'.format(exp_name)

    if dfps.chrom_folder_process_sequence(source_folder, store_folder,
                                         cond_file, analysis_details,
                                         copy_analysis = analysis_cpy,
                                         copy_conditions = conditions_cpy):
        pass
    else:
        pass_list.append(exp_name)
        continue

print('Passed:')
if len(pass_list) == 0:
    print('none')
else:
    [print('{}'.format(x)) for x in pass_list]

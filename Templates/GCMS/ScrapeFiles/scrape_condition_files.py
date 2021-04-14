from ChromProcess import Classes
import shutil
from pathlib import Path
import os

exp_paths = Classes.DataPaths(r'C:\Users\willi\Documents\Data\Data_information.csv')
storage_stem = Path(r'C:\Users\willi\Documents\Data')
backup_location = Path(r'Z:\Huck\Prebiotic\FormoseReaction\Data\Conditions_backup')

path_list = [e for e in exp_paths.paths]

start_token = 'threose_1_2021'
end_token = ''
idx_s = 0
idx_e = len(path_list)
for c,p in enumerate(exp_paths.experiment_codes):
    if p == start_token:
        idx_s = c
    if p == end_token:
        idx_e = c + 1

path_list = path_list[idx_s:idx_e]

for e in path_list:

    store_folder = Path(storage_stem/e.data_type/'FRN'/e.experiment_code)
    os.makedirs(store_folder, exist_ok = True)

    # Load in experiment conditions file from the drive
    for f in os.listdir(e.path):
        if 'conditions' in f and not f.startswith("."):
            info = Classes.Information(e.path/f)
            cond_file = Classes.Experiment_Conditions()
            cond_file.read_from_file(e.path/f)
            cond_file.experiment_code = e.experiment_code

    # derive info from the file
    series_regions = info.regions
    internal_ref_region = info.internal_ref_region
    use_MS = True # Whether to extract mass spectra or not
    analysis_type = e.data_type
    i_thres_MS = 500 # Threshold for removing background mass spectra
    peak_pick_thres = 0.1 # Threshold for peak picking as a proportion of the max
    dil = info.dilution
    IS_c = info.internal_ref_concentration

    # create an analysis object
    analysis_file = Classes.Analysis_Information(regions= series_regions,
                        IS_region = internal_ref_region, use_MS = use_MS,
                        analysis_type = analysis_type, i_thres_MS = i_thres_MS,
                        peak_pick_thres = peak_pick_thres,
                        exp_name = e.experiment_code,
                        dilution_factor = dil, dil_err = 4.47E-02,
                        IS_conc = IS_c, IS_err = IS_c/100)

    analysis_file.write_to_file(directory = store_folder)
    cond_file.write_to_file(directory = store_folder)

    print(e.experiment_code)

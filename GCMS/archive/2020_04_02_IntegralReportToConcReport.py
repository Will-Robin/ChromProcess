import sys
sys.path.append(r'/Users/williamrobinson/Documents/Nijmegen/2020_03_13_Packages')
from ChromProcess import Classes, os, file_output, file_import, info_params, np, plt, processing_functions
from NetFit.file_in_out import get_data

def apply_calibration_to_dataset(dataset, calib_dict):
    '''
    Modifies dataset object, converting integrals into concenctrations
    '''
    concs = {}
    for d in dataset.dependents:
        if d in calib_dict:
            if calib_dict[d]["A"] == 0:
                conversion = lambda x : dataset.dilution*dataset.internal_ref_concentration*(x-calib_dict[d]["C"])/calib_dict[d]["B"]
            else:
                conversion = lambda x : dataset.dilution*dataset.internal_ref_concentration*(-calib_dict[d]["B"] + np.sqrt((calib_dict[d]["B"]**2) - (4*calib_dict[d]["A"]*(calib_dict[d]["C"]-x))))/(2*calib_dict[d]["A"])

            dataset.dependents[d] = np.nan_to_num(conversion(dataset.dependents[d]))

directory = r"/Users/williamrobinson/Documents/Nijmegen/Data/GCMS/FRN073"

calib_file = r"/Users/williamrobinson/Documents/Nijmegen/2020_03_13_Packages/info_files/2020_02_18_Combined_Species_Info.csv"

information = Classes.Calibration_File(calib_file, type = "GCMS") # get calibration and assignment information form file.

os.chdir(directory)
os.chdir("Data_reports")
print(os.getcwd())

f_list = os.listdir()
for f in f_list:
    if "all_integrals_assigned" in f:
        data = get_data(f, x_axis_key = "time/ s", flow_data = True)

apply_calibration_to_dataset(data, information.calibrations)


with open("{}_manual_assigned_conc_series.csv".format(data.name), "w") as outfile:
    outfile.write("Dataset,{}\n".format(data.name))

    outfile.write("start_conditions\n")

    for c in data.conditions:
        outfile.write("{},".format(c))
        [outfile.write("{},".format(x)) for x in data.conditions[c]]
        outfile.write("\n")
    outfile.write("end_conditions\n")

    outfile.write("start_analysis_details\n")

    outfile.write("Calibrations_bounds_date,{}\n".format(information.date))
    outfile.write("Chromatography_method,{},{}\n".format(information.type, information.method))
    outfile.write("Derivatisation_method,{}\n".format(information.derivatisation))

    outfile.write("end_analysis_details\n")

    outfile.write("start_data\n")

    out = [data.time]
    p_header = ["time/ s"]

    for p in data.dependents:
        label_name = p
        if processing_functions.int_test(label_name[-1]):
            label_name = label_name[:-2]

        species_smiles = info_params.canonical_SMILES[label_name]
        p_header.append(species_smiles+" M")
        out.append(data.dependents[p])

    out = [list(i) for i in zip(*out)]

    [outfile.write("{},".format(x)) for x in p_header]
    outfile.write("\n")
    for o in out:
        [outfile.write("{},".format(x)) for x in o]
        outfile.write("\n")

    outfile.write("end_data\n")

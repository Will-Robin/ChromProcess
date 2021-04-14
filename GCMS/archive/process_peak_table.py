import sys
sys.path.append(r'C:\Users\willi\Documents\Packages')
from ChromProcess import Classes, os, plt, series_builder, plotting, file_output, file_import
from ChromProcess import mass_spectra, calibration_functions, info_params
from ChromProcess import processing_functions
from ChromProcess import series_operations
import pickle
import numpy as np

def write_report(series_name, series_info, information, output_dict):

    with open("{}.csv".format(series_name), "w") as outfile:
        outfile.write("Dataset,{}\n".format(series_name))

        outfile.write("start_conditions\n")
        for c in series_info["conditions"]:
            outfile.write("{},".format(c))
            [outfile.write("{},".format(x)) for x in series_info["conditions"][c]]
            outfile.write("\n")
        outfile.write("end_conditions\n")

        outfile.write("start_analysis_details\n")

        outfile.write("Calibrations_bounds_date,{}\n".format(information.date))
        outfile.write("Chromatography_method,{},{}\n".format(information.type, information.method))
        outfile.write("Derivatisation_method,{}\n".format(information.derivatisation))

        outfile.write("end_analysis_details\n")

        outfile.write("start_data\n")
        [outfile.write("{},".format(o)) for o in output_dict]
        outfile.write("\n")

        out = []
        for o in output_dict:
            out.append(output_dict[o])

        out = [list(i) for i in zip(*out)]

        for o in out:
            [outfile.write("{},".format(x)) for x in o]
            outfile.write("\n")

        outfile.write("end_data\n")

def get_information_file(information_file):

    info_dict = {}

    with open(information_file, 'r') as f:
        for line in f:
            if "Dataset" in line:
                ins = line.strip("\n")
                info_dict["set_name"] = ins.split(",")[1]
            if "dilution_factor" in line:
                ins = line.strip("\n")
                info_dict["dilution"] = float(ins.split(",")[1])
            if "series_values" in line:
                ins = line.strip("\n")
                spl = ins.split(",")
                info_dict["x_series"] =  [float(x) for x in spl[1:] if x != ""]
            if "series_unit" in line:
                ins = line.strip("\n")
                info_dict["x_name"] = ins.split(",")[1]
            if "series_regions" in line:
                ins = line.strip("\n")
                reg = [float(x) for x in ins.split(",")[1:] if x != ""]
                info_dict["regions"] = [reg[x:x+2] for x in range(0, len(reg), 2)]
            if "internal_ref_region" in line:
                ins = line.strip("\n")
                reg = [float(x) for x in ins.split(",")[1:] if x != ""]
                info_dict["internal_ref_region"] = [reg[x:x+2] for x in range(0, len(reg), 2) if x != ""][0]
            if "internal_ref_concentration" in line:
                ins = line.strip("\n")
                info_dict["internal_ref_concentration"] = float(ins.split(",")[1])

        condset = []
        readstate = False
        with open(information_file, "r") as f:
            for c,line in enumerate(f):
                if "start_conditions" in line:
                    readstate = True
                    line = next(f)
                if "end_conditions" in line:
                    readstate = False
                if readstate:
                    newline = line.strip("\n")
                    condset.append([x for x in newline.split(",") if x != ""])
        c_out = {}
        for c in condset:
            c_out[c[0]] = [float(x) for x in c[1:]]

        info_dict["conditions"] = c_out

    return info_dict

def read_peak_table(file):
    chrom_table = {}
    with open(file, "r") as f:
        for c,line in enumerate(f):
            if c == 0:
                value = float(line.split(",")[1])
            elif c == 2:
                int_ref = [float(x) for x in line.strip("\n").split(",")]
            elif c < 3:
                pass
            else:
                rd = [float(x) for x in line.strip("\n").split(",")]
                chrom_table[round(rd[0],3)-int_ref[0]] = {"integral":rd[1], "start":rd[2], "end":rd[3]}

    return value, chrom_table

def apply_calibration_to_array(array, calibraton_factors, IS_conc, dilution_factor):
    if calibraton_factors["A"] == 0:
        conversion = lambda x : dilution_factor*IS_conc*(x-calibraton_factors["C"])/calibraton_factors["B"]
    elif calibraton_factors["A"] > 0:
        conversion = lambda x : dilution_factor*IS_conc*(-calibraton_factors["B"] + np.sqrt((calibraton_factors["B"]**2) - (4*calibraton_factors["A"]*(calibraton_factors["C"]-x))))/(2*calibraton_factors["A"])
    elif calibraton_factors["A"] < 0:
        conversion = lambda x : dilution_factor*IS_conc*(-calibraton_factors["B"] + np.sqrt((calibraton_factors["B"]**2) - (4*calibraton_factors["A"]*(calibraton_factors["C"]-x))))/(2*calibraton_factors["A"])

    concs = conversion(np.array(array))

    return concs

series_name = "FRN073"
IS_pos = 6.91
directory = r"C:\Users\willi\Documents\Data\GCMS\FRN\{}".format(series_name)

calib_file = r"C:\Users\willi\Documents\Packages\info_files\2020_02_18_Combined_Species_Info.csv"

information = Classes.Calibration_File(calib_file, type = "GCMS") # get calibration and assignment information form file.

for i in information.boundaries:
    information.boundaries[i] = [x - IS_pos for x in information.boundaries[i] ]

os.chdir(directory)

for f in os.listdir():
    if "conditions" in f:
        series_info = get_information_file(f)

os.chdir(r"Peak_tables")

chroms = {}
time = []
for f in os.listdir():
    if f.endswith(".csv"):
        t, peaks = read_peak_table(f)
        chroms[f.split(".")[0]] = peaks
        time.append(t)

integral_series = {}

for count,c in enumerate(chroms):
    for peak in chroms[c]:
        if peak in integral_series:
            integral_series[peak][count] = chroms[c][peak]["integral"]
        else:
            integral_series[peak] = np.zeros(len(chroms))
            integral_series[peak][count] = chroms[c][peak]["integral"]

integral_series = series_operations.bin_integral_series(integral_series, bound = 0.01)

output = {}
integrals = {}
integrals["time/ s"] = series_info["x_series"]

for i in integral_series:
    name = processing_functions.name_peak(i,information.boundaries)
    integrals[name + " ({}) integral".format(i+IS_pos)] = integral_series[i]
    if name in information.calibrations:

        y = apply_calibration_to_array(integral_series[i], information.calibrations[name], series_info["internal_ref_concentration"], series_info["dilution"])

        if info_params.canonical_SMILES[name.split(" ")[0]]+ " M" in output:
            output[info_params.canonical_SMILES[name.split(" ")[0]] + " {} M".format(1+[*output].count(info_params.canonical_SMILES[name.split(" ")[0]] + " M"))] = y
        else:
            output[info_params.canonical_SMILES[name.split(" ")[0]] + " M"] = y

        plt.plot(series_info["x_series"], y, c = info_params.colour_assignments[name], label = name)

    else:
        pass

plt.legend()
plt.title("Assigned")
plt.show()
plt.close()

for i in integral_series:
    name = processing_functions.name_peak(i,information.boundaries)
    if name in info_params.colour_assignments:
        plt.plot(series_info["x_series"], integral_series[i], c = info_params.colour_assignments[name], label = name)
    else:
        plt.plot(series_info["x_series"], integral_series[i], c = "r", linewidth = 5, label = i + IS_pos)

plt.legend()
plt.title("Unassigned")
plt.show()
plt.close()

sort_comps = sorted([*output])
output_2 = {}
output_2["time/ s"] = series_info["x_series"]
for s in sort_comps:
    output_2[s] = output[s]

os.makedirs(r"{}\Data_reports".format(directory), exist_ok = True)
os.chdir(r"{}\Data_reports".format(directory))
write_report("{}_concentrations".format(series_name), series_info, information, output_2)

write_report("{}_integrals".format(series_name), series_info, information, integrals)

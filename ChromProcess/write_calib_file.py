import Chromatography_processing as CP
from Chromatography_processing import GCMS_bounds
from Chromatography_processing import GCMS_calibrations

from Chromatography_processing import HPLC_bounds
from Chromatography_processing import HPLC_calibrations

with open("2019_12_02_GCMS_Calibrations.csv", "w") as f:
    f.write("Date,2019_12_02\n")
    f.write("Chromatography_type,GCMS\n")
    f.write("Method_file,Method 5\n")
    f.write("Derivatisation_details, 75 uL EtONH2.HCl; 30 mins @ 70 oC; 25 uL N,O-Bis(trimethylsilyl)trifluoroacetamide; 30 mins @ 70 oC; dilution to 200 uL with internal standard solution\n")
    f.write("internal_standard_concentration,0.0008\n")
    f.write("start_calibration\n")
    f.write("compound,gradient/ integral_unit/ M,intecept/ integral_unit \n")

    for x in GCMS_calibrations:
        f.write("{},{},{}\n".format(x, GCMS_calibrations[x][0],GCMS_calibrations[x][1]))

    f.write("end_calibration\n")
    f.write("start_boundaries\n")

    for x in GCMS_bounds:
        f.write("{},{},{}\n".format(x, GCMS_bounds[x][0], GCMS_bounds[x][1]))
    f.write("end_boundaries\n")

with open("2019_12_02_HPLC_Calibrations.csv", "w") as f:
    f.write("Date,2019_12_02\n")
    f.write("Chromatography_type,HPLC\n")
    f.write("Method_file, 2019_07_31_GIST_2\n")
    f.write("Derivatisation_details,60% MeCN saturated DNPH; 0.5% 2 M HCl; 19.5%MeCN; 20% sample; 30 mins\n")
    f.write("start_calibration\n")
    f.write("compound,gradient/ integral_unit/ M,intecept/ integral_unit \n")

    for x in HPLC_calibrations:
        f.write("{},{},{}\n".format(x, HPLC_calibrations[x][0],HPLC_calibrations[x][1]))

    f.write("end_calibration\n")
    f.write("start_boundaries\n")

    for x in HPLC_bounds:
        f.write("{},{},{}\n".format(x, HPLC_bounds[x][0], HPLC_bounds[x][1]))
    f.write("end_boundaries\n")

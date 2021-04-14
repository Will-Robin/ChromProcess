import sys
sys.path.append(r'C:\Users\willi\Documents\Packages')
from ChromProcess import Classes, os, plt, series_builder, plotting, file_output, file_import
from ChromProcess import mass_spectra, calibration_functions, info_params
from ChromProcess.series_operations import cluster

path_list = r"C:\Users\willi\Documents\Data\GCMS\2020_calibrations\Calibration_paths.csv"

path_dict = {}
wd= os.getcwd()
with open(path_list, "r") as f:
    for c,line in enumerate(f):
        if c > 0:
            x = line.strip("\n").split(",")
            path_dict[x[0]]= [z for z in x[1:] if z != ""]

calib_file = r"C:\Users\willi\Documents\Packages\info_files\2020_03_16_Combined_Species_Info.csv"
information = Classes.Calibration_File(calib_file, type = "GCMS") # get calibration and assignment information form file.
information.calibrations = {"none":[0,1]}

thold = 0.1 # threshold for peak detection as a fraction of the signal maximum
min_inten = 5000 # to cut off noise if needed.
use_mass_spec = True
#path_dict = {"X": "/Users/williamrobinson/Documents/Nijmegen/2020_calibrations/2020_03_12/threose/C"}
for p in path_dict:

    for directory in path_dict[p]:

        chroms, cond_file = file_import.load_cdf_from_directory(directory, ms = use_mass_spec) # get experimental data and conditions

        series = Classes.Chromatogram_Series(chroms, cond_file) # Convert chromatogram list and conditions into a Series object

        calibration_functions.calibration_series(series) # Convert the series into a calibration series (since experiment info is slightly different)
        '''
        for c in chroms: # plotting chromatograms for inspection.
            plt.plot(c.time,c.signal)
        plt.show()
        '''

        # Call a series of functions on the Series object.
        series_builder.series_from_chromatogram_set(series, information,
                                                    threshold = thold,
                                                    min_intensity = min_inten,
                                                    combine_assignments = False,
                                                    plot_mass_spectra = use_mass_spec)

        plotting.plot_peak_ion_chromatogram_from_series(series)
        plotting.write_peak_ion_chromatogram_from_series(series)
        #plotting.plot_ion_channels(series)
        '''
        combine_TMS = calibration_functions.combine_mz_73(series)
        with open("{}_combined_TMS.csv".format(series.set_name), "w") as f:
            f.write("{},I(m/z == 73)/I(stand)\n".format(series.x_name))
            for c,v in enumerate(series.x_series):
                f.write("{},{}\n".format(v,combine_TMS[c]))
        '''
        for c in series.chromatograms:
            print(c.internal_reference.retention_time)

        os.makedirs("Calibration_output", exist_ok = True)
        os.chdir("Calibration_output")
        cabs_1 = calibration_functions.fit_calibration_curve_TIC(series) # Get calibration factors from the total ion chromatogram.
        cabs_2 = calibration_functions.fit_calibration_curve_IC(series) # Get calibration factors from ion chromatograms.

        with open("calibration_params.csv", "w") as f:

            for p in cabs_1:
                f.write("Average_retention_time/ min,{}\n".format(p))
                f.write("class,A,B,C\n")
                f.write("{},".format("TIC"))
                [f.write("{},".format(c)) for c in cabs_1[p]]
                f.write("\n")
                for m in cabs_2[p]:
                    if m in info_params.frag_assignments:
                        f.write("{},".format(m))
                    else:
                        f.write("{},".format(m))
                    [f.write("{},".format(c)) for c in cabs_2[p][m]]
                    f.write("\n")
        calibration_functions.plot_calibration_curves_TIC(series,cabs_1)
        calibration_functions.plot_calibration_curves_IC(series,cabs_2)
        mass_spectra.ion_chromatograms_relative_to_peak(series)
        print("Finished.")

import sys
sys.path.append(r'C:\Users\willi\Documents\Packages')
from ChromProcess import Classes, os, plt, series_builder, plotting, file_output
from ChromProcess import np
from ChromProcess import mass_spectra

calib_file = r"C:\Users\willi\Documents\Packages\ChromProcess\Calibrations\2019_12_02_GCMS_Calibrations.csv"

information = Classes.Calibration_File(calib_file)

os.chdir(r"C:\Users\willi\Documents\Data\GCMS\2020_calibrations\All_sugar_overlay")# changes the directory to where files are stored

filelist = os.listdir() # get list of files in directory
filelist.sort() # sort the files

'''Writes the chromatograms to csv files and creates Chromatogram objects from
.cdf files'''

chroms = []
fname = "C7_glycoC7_5.cdf"
chroms.append(Classes.Chromatogram(fname, mass_spec = True))

for c in chroms:
    plt.plot(c.time,c.signal)
plt.show()

for f in filelist:
    if "conditions" in f:
        cond_file = f

'''Create a Chromatogram_Series object from the list of chromatograms.
Further information about the series is input via a .csv file
(see provided example)'''
series = Classes.Chromatogram_Series(chroms, cond_file)

os.makedirs("{}".format(fname.split(".")[0]), exist_ok = True)
os.chdir("{}".format(fname.split(".")[0]))

thold = 0.1 # threshold for peak detection as a fraction of the signal maximum
min_inten = 2000 # to cut off noise if needed.
series_builder.series_from_chromatogram_set(series, information, threshold = thold, min_intensity = min_inten, combine_assignments = False)

mass_spectra.get_peak_mass_spectra(series)

os.makedirs("Peak_mass_spectra", exist_ok = True)
os.chdir("Peak_mass_spectra")

mass_spectra.write_mass_spectra(series)
mass_spectra.mass_sepctrum_peak_picker(series, mass_split = 7)

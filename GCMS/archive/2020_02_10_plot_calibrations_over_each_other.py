import sys
sys.path.append(r'/Users/williamrobinson/Documents/Nijmegen/2020_03_13_Packages')
from ChromProcess import Classes, os, plt, series_builder, plotting, file_output, info_params

os.chdir(r"/Users/williamrobinson/Documents/Nijmegen/2020_calibrations/sugar_overlay_after_feb") # changes the directory to where files are stored

calib_file = r"/Users/williamrobinson/Documents/Nijmegen/2020_03_13_Packages/info_files/2020_03_16_Combined_Species_Info.csv"
information = Classes.Calibration_File(calib_file, type = "GCMS")

filelist = os.listdir() # get list of files in directory
filelist.sort() # sort the files

'''Writes the chromatograms to csv files and creates Chromatogram objects from
.cdf files'''
chroms = []
names = []
specs = []
for f in filelist:
    if f.endswith(".cdf"):
        print("loading", f)
        chroms.append(Classes.Chromatogram(f, mass_spec = True))
        names.append(f.split("_")[0])
        specs.append(f[3:-6])
    if "conditions" in f:
        cond_file = f

'''Create a Chromatogram_Series object from the list of chromatograms.
Further information about the series is input via a .csv file
(see provided example)'''

series = Classes.Chromatogram_Series(chroms, cond_file)
series.x_series = names

for c in chroms:
    plt.plot(c.time,c.signal, label  = c.filename)
plt.legend()
plt.show()

series.x_series = [s for s in range(len(specs))]
print(specs)
plotting.region_heatmap(series, information)

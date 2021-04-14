import sys
sys.path.append(r'C:\Users\willi\Documents\Packages')
from ChromProcess import Classes, os, plt, series_builder, plotting, file_output

directory = r"C:\Users\willi\Documents\Data\GCMS\FRN\FRN087"

calib_file = r"C:\Users\willi\Documents\Packages\ChromProcess\Calibrations\2019_12_02_GCMS_Calibrations.csv"
information = Classes.Calibration_File(calib_file)
os.chdir(directory) # changes the directory to where files are stored
filelist = os.listdir() # get list of files in directory
filelist.sort() # sort the files

chroms = []
for f in filelist:
    if f.endswith(".cdf"):
        print("loading", f)
        chroms.append(Classes.Chromatogram(f, mass_spec = False))
    if "conditions" in f:
        cond_file = f

'''Create a Chromatogram_Series object from the list of chromatograms.
Further information about the series is input via a .csv file
(see provided example)'''

series = Classes.Chromatogram_Series(chroms, cond_file)

for c in chroms:
    plt.plot(c.time,c.signal)
plt.show()

thold = 0.1 # threshold for peak detection as a fraction of the signal maximum
min_inten = 4000 # to cut off noise if needed.

'''Get plots and spreadsheets of the peaks picked from each region of each
chromatogram'''
series_builder.series_from_chromatogram_set(series, information, threshold = thold, min_intensity = min_inten, combine_assignments = False)

os.makedirs("assigned_chromatograms", exist_ok = True)
os.chdir("assigned_chromatograms")
for c in chroms:
    plotting.plot_chromatogram_with_assignments(c,series,information.boundaries, labelling = True)
os.chdir("..")

os.makedirs("Time_progresses", exist_ok = True)
os.chdir("Time_progresses")
plotting.region_heatmap(series, information)
os.chdir("..")

'''Plot each chromatogram in each region of the series in a new folder.'''
os.makedirs("Chromatograms", exist_ok = True)
os.chdir("Chromatograms")
for c in chroms:
    print("writing chromatogram {} to csv file".format(c.filename))
    file_output.chromatogram_to_csv_GCMS(c)
plotting.plot_chromatograms(series)
print("Plotted chromatograms.")
os.chdir("..")

os.makedirs("Mass_Spectra", exist_ok = True)
os.chdir("Mass_Spectra")

print("Finished.")

import sys
sys.path.append(r'C:\Users\willi\Documents\Packages')
from ChromProcess import Classes, os, plt, series_builder, plotting, file_output

def calibration_series(series):

    '''
    Replace sample numbers with concentrations.
    '''
    stock_conc = (series.conditions["Mass/ g"][0]/ series.conditions["Mr/ g/mol"][0])/series.conditions["Stock_vol/ L"][0]

    mol_series = [stock_conc*x for x in series.conditions["Titration_series/ L"]]

    series_conc = [x/series.conditions["Sample_volume/ L"][0] for x in mol_series]

    series.x_series = series_conc

    series.x_name = "concentration/ M"

    series.conditions["internal_ref_concentration/ M"] = [series.internal_ref_concentration]

    for c,v in enumerate(series.chromatograms):
        v.timepoint = series.x_series[c]

calib_file = r"C:\Users\willi\Documents\Packages\ChromProcess\Calibrations\2020_02_18_GCMS_Calibrations.csv"

information = Classes.Calibration_File(calib_file)

os.chdir(r"C:\Users\willi\Documents\Data\GCMS\2020_calibrations\2020_05_27\arabitol")# changes the directory to where files are stored
sugar_name = "arabitol"

keys = [*information.boundaries]
for i in keys:
    if sugar_name in i:
        pass
    else:
        del information.boundaries[i]


filelist = os.listdir() # get list of files in directory
filelist.sort() # sort the files

'''Writes the chromatograms to csv files and creates Chromatogram objects from
.cdf files'''
chroms = []

for f in filelist:
    if f.endswith(".cdf"):
        print("loading", f)
        chroms.append(Classes.Chromatogram(f, mass_spec = True))
    if "conditions" in f:
        cond_file = f

'''Create a Chromatogram_Series object from the list of chromatograms.
Further information about the series is input via a .csv file
(see provided example)'''

series = Classes.Chromatogram_Series(chroms, cond_file)
calibration_series(series)

for c in chroms:
    plt.plot(c.time,c.signal, label = c.filename)
plt.legend()
plt.show()

thold = 0.1 # threshold for peak detection as a fraction of the signal maximum
min_inten = 2000 # to cut off noise if needed.

'''Get plots and spreadsheets of the peaks picked from each region of each
chromatogram'''
series_builder.series_from_chromatogram_set(series, information, threshold = thold, min_intensity = min_inten, combine_assignments = False)

plotting.region_heatmap(series, information)

'''Plot each chromatogram in each region of the series in a new folder.'''
if not os.path.exists("Chromatograms"):
    os.mkdir("Chromatograms")

os.chdir("Chromatograms")

for c in chroms:
    print("writing chromatogram {} to csv file".format(c.filename))
    file_output.chromatogram_to_csv_GCMS(c)

plotting.plot_chromatograms(series)
print("Plotted chromatograms.")

os.chdir("..")
if not os.path.exists("Mass_Spectra"):
    os.mkdir("Mass_Spectra")
os.chdir("Mass_Spectra")


print("Finished.")

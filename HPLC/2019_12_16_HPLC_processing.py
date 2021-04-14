import sys
sys.path.append(r'/Users/williamrobinson/documents/nijmegen/packages')
from ChromProcess import Classes, os, file_import, plt, series_builder, plotting, file_output

folder = r"/Users/williamrobinson/Documents/Nijmegen/Data/HPLC/FRN108" # You  can paste the folder path here, but it must be preceded by an 'r' (in Windows)

calib_file = r"/Users/williamrobinson/Documents/Nijmegen/Packages/info_files/2020_03_16_Combined_Species_Info.csv"

information = Classes.Calibration_File(calib_file, type = "HPLC")
information.boundaries["formaldehyde"] = [10,11]

'''
# X3 calibration 01/2021
'''

information.calibrations = {'formaldehyde'    : {'A': 0.0, 'B': 10965324.05, 'C': 777.1823564},
                            'dihydroxyacetone' : {'A':0.0, 'B':16980777.48, 'C':	1432.809754}}
information.boundaries = {"formaldehyde":[10.4,11], 'dihydroxyacetone':[4.8,5.8]}


os.chdir(folder) # changes the working directory
drlist = os.listdir() # get list of files in directory
drlist.sort() # sorts the files

'''Please make sure that your files are named so that they are read in the right
order to be aligned with the time or concentration series array. Otherwise the
automatic plotting/ series generation will not work properly!'''

chromatograms = []
for f in drlist:
    if f.endswith(".txt"): # this is how the programs finds the data files: .txt is assumed to be a HPLC file, .cdf is assumed to be GCMS
        print(f)
        chromatograms.append(Classes.Chromatogram(f,channel_select = '215nm'))
    if "conditions" in f:
        cond_file = f

n = 1
fig,ax= plt.subplots(nrows = 1, sharex=True)
for count,c in enumerate(chromatograms):
     ax.plot(c.time, c.signal, label = n)
     n+=1
plt.legend()
plt.show()

'''Create a series object from the list of chromatograms and an information file.'''
Series = Classes.Chromatogram_Series(chromatograms, cond_file)

'''Thresholds for peak detection'''
thold  = 0.01 # threshold for peak detection as a fraction of the signal maximum
maxten = 1e100 # peaks above this level will be ignored. Set so that the DNPH peak is ignored
minten = -1e100
'''
Get plots and spreadsheets of the peaks picked from each region of each chromatogram
'''
series_builder.series_from_chromatogram_set(Series, information, threshold = thold, max_intensity = maxten, min_intensity = minten,cluster_dev = 0.05)
plotting.region_heatmap(Series,information)

'''
Plot each chromatogram in each region of the series.
'''
#os.makedirs("assigned_chromatograms", exist_ok = True)
#os.chdir("assigned_chromatograms")
#for c in chromatograms:
#    plotting.plot_chromatogram_with_assignments(c,Series,information.boundaries)
#os.chdir("..")

if not os.path.exists("Chromatograms"):
    os.mkdir("Chromatograms")

os.chdir("Chromatograms")

for c in chromatograms:
    print("writing chromatogram {} to csv file".format(c.filename))
    file_output.chromatogram_to_csv_HPLC(c)

plotting.plot_chromatograms(Series,information)

print("Finished.")

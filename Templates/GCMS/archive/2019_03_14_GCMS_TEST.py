import sys
sys.path.append(r'/Users/williamrobinson/Documents/Nijmegen/2020_03_13_Packages/')
from ChromProcess import Classes, os, plt, series_builder, plotting, file_output, processing_functions
from ChromProcess import mass_spectra, np
from scipy.stats import linregress

directory = r"/Users/williamrobinson/Documents/Nijmegen/2020_calibrations/ribose"

os.chdir(directory) # changes the directory to where files are stored
chroms = []
f_list = os.listdir()
for f in f_list:
    if f.endswith(".cdf"):
        chroms.append(Classes.Chromatogram(f, mass_spec = True))
    else:
        pass

region = [10.9,11.2]

x_axis = [0.1,0.25,0.5,1,1.5,2.0]
integrals = {}

for c,chromatogram in enumerate(chroms):
    '''
    1. Get peak parameters from TIC
    '''
    region_inds = np.where((chromatogram.time > region[0]) & (chromatogram.time < region[1]))[0]
    peak = processing_functions.Peak_finder(chromatogram.signal[region_inds], thres= 0.1, min_dist=0.1)

    print("peak indices",peak)
    time = chromatogram.time[region_inds]
    signal  = chromatogram.signal[region_inds]
    rt    = peak['Peak_indices'][0]
    start = peak['Peak_start_indices'][0]
    end   = peak['Peak_end_indices'][0]

    '''
    2. Extract ion chromatograms
    '''
    mass_spectra.get_ion_chromatograms(series)


    integral_track = {}
    for m in ion_chromatgrams:
        if np.trapz(ion_chromatgrams[m]) > 30000:
            if m in integrals:
                integrals[m][c] = np.trapz(ion_chromatgrams[m])
            else:
                integrals[m] = np.zeros(len(x_axis))
                integrals[m][c] = np.trapz(ion_chromatgrams[m])

line_vals = {}
for i in integrals:
    coef = linregress(x_axis[:-1], integrals[i][:-1])
    line = [x*coef[0] + coef[1] for x in x_axis[:-1]]
    line_vals[i] = [coef[0],coef[1]]
    # poly1d_fn is now a function which takes in x and returns an estimate for y
    plt.scatter(x_axis[:-1], integrals[i][:-1])
    plt.plot(x_axis[:-1], line, label = i)

plt.legend()
plt.show()

with open("cal_by_mass.csv", "w") as f:
    f.write("mass,grad,intercept\n")
    for i in line_vals:
        f.write("{},{},{}\n".format(i,line_vals[i][0],line_vals[i][1]))

print(integrals)

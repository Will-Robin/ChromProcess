import sys
sys.path.append(r'C:\Users\willi\Documents\Packages')
from ChromProcess import Classes, os, plt, series_builder, plotting, file_output, file_import
from ChromProcess import mass_spectra, calibration_functions, info_params
from ChromProcess import processing_functions
import pickle

directory = r"C:\Users\willi\Documents\Data\GCMS\FRN\FRN100\A"

calib_file = r"C:\Users\willi\Documents\Packages\info_files\2020_03_16_Combined_Species_Info.csv"

thold = 0.05 # threshold for peak detection as a fraction of the signal maximum
min_inten = 10000 # to cut off noise if needed.
use_mass_spectra = True

number_of_peaks = 2

information = Classes.Calibration_File(calib_file, type = "GCMS") # get calibration and assignment information form file.

chroms, cond_file = file_import.load_cdf_from_directory(directory, ms = use_mass_spectra) # get experimental data and conditions

for c in chroms: # plotting chromatograms for inspection.
    #c.signal = processing_functions.savitzky_golay(c.signal, 7, 3, deriv=0, rate=1)
    processing_functions.MS_intensity_threshold_chromatogram(c, threshold = 500)
    plt.plot(c.time,c.signal)
plt.show()


os.chdir("deconvolution")

for f in os.listdir():
    if "conditions" in f:
        cond_file = f


from ChromProcess import deconvolution
from ChromProcess import np
from ChromProcess import series_operations

series = Classes.Chromatogram_Series(chroms, cond_file) # Convert chromatogram list and conditions into a Series object

series_operations.pick_peaks(series, threshold = 0.1, max_intensity = 1e100, min_intensity = min_inten)
os.makedirs("plots",exist_ok = True)
os.chdir("plots")

series_operations.get_internal_ref_integrals(series)

series_integrals = []
for count,c in enumerate(series.chromatograms):

    for r in series.regions:
        inds = np.where((c.time > r[0])&(c.time < r[1]))[0]

        popt,pcov = deconvolution.deconvolute_region(c, r, num_peaks = number_of_peaks)
        indivs = [popt[x:x+3] for x in range(0, len(popt), 3)]

        indivs = sorted(indivs, key = lambda x:x[1])
        sum = np.zeros(len(c.time[inds]))

        for i in indivs:

            component = deconvolution._1gaussian(c.time[inds], *i)
            plt.plot(c.time[inds],component, c = 'b', label = "component")
            sum += component
            integral = np.trapz(component, x = c.time[inds])/c.internal_reference.integral

            if i[1] in series.deconvoluted_series:
                series.deconvoluted_series[i[1]][count] = integral
            else:
                series.deconvoluted_series[i[1]] = [0.0 for x in series.chromatograms]
                series.deconvoluted_series[i[1]][count] = integral

        plt.plot(c.time[inds], sum, c = 'r', label = "sum components")
        plt.plot(c.time[inds], c.signal[inds], c = "k", label = "data")
        plt.legend()
        plt.savefig("{}".format(c.filename))
        plt.close()

series.deconvoluted_series = series_operations.bin_integral_series(series.deconvoluted_series,bound = 0.02)

series_operations.apply_calibration_to_deconv_series(series,information.calibrations,information.boundaries)
file_output.write_deconvoluted_output(series, information)

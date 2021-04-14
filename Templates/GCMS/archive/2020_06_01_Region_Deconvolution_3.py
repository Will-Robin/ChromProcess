from sklearn.decomposition import PCA

import sys
sys.path.append(r'C:\Users\willi\Documents\Packages')
from ChromProcess import Classes, os, plt, series_builder, plotting, file_output, file_import
from ChromProcess import mass_spectra, calibration_functions, info_params
from ChromProcess import processing_functions
import pickle

from ChromProcess import deconvolution
from ChromProcess import np
from ChromProcess import series_operations


directory = r"C:\Users\willi\Documents\Data\GCMS\FRN\FRN089"

calib_file = r"C:\Users\willi\Documents\Packages\info_files\2020_03_16_Combined_Species_Info.csv"

thold = 0.05 # threshold for peak detection as a fraction of the signal maximum
min_inten = 5e4 # to cut off noise if needed.
use_mass_spectra = False
information = Classes.Calibration_File(calib_file, type = "GCMS") # get calibration and assignment information form file.

chroms, cond_file = file_import.load_cdf_from_directory(directory, ms = use_mass_spectra) # get experimental data and conditions
'''
for c in chroms: # plotting chromatograms for inspection.
    plt.plot(c.time,c.signal)
plt.show()
'''

os.chdir("deconvolution")

for f in os.listdir():
    if "conditions" in f:
        cond_file = f

series = Classes.Chromatogram_Series(chroms, cond_file) # Convert chromatogram list and conditions into a Series object

os.makedirs("plots",exist_ok = True)
os.chdir("plots")

stacking = np.array([])

min_length = 1e100
for c in series.chromatograms:
    if len(c.time) < min_length:
        min_length = len(c.time)
        time = c.time

ser = np.array(series.x_series)
time_region = np.where((ser > 0)&(ser < 12000))[0]

heat_stack = np.empty(( len(series.chromatograms), len(series.chromatograms[0].signal[:min_length] )))

for c in range(0,len(series.chromatograms)):
    heat_stack[c] = series.chromatograms[c].signal[:min_length]

reg = series.regions[0]
inds = np.where((time > reg[0])&(time < reg[1]))[0]
time = time[inds]

heat_stack =  heat_stack[time_region]
stacking = heat_stack[:,inds]

import matplotlib.pyplot as plt
from tensorly.base import tensor_to_vec, partial_tensor_to_vec
from tensorly.datasets.synthetic import gen_image
from tensorly.random import check_random_state
from tensorly.regression.kruskal_regression import KruskalRegressor
import tensorly as tl

# shape of the images
patterns = [stacking]
# ranks to test
ranks = [1, 2, 3, 4, 5]

# Generate random samples
rng = check_random_state(1)
X = tl.tensor(rng.normal(size=(2400,stacking.shape[0],stacking.shape[1]), loc=0, scale=1))

# Parameters of the plot, deduced from the data
n_rows = 1
n_columns = len(ranks) + 1
# Plot the three images
fig = plt.figure()

for i, pattern in enumerate(patterns):

    # Generate the original image
    weight_img = tl.tensor(pattern)


    # Generate the labels
    print(tensor_to_vec(weight_img).shape)
    y = tl.dot(partial_tensor_to_vec(X, skip_begin=1), tensor_to_vec(weight_img))

    # Plot the original weights
    ax = fig.add_subplot(n_rows, n_columns, i*n_columns + 1)
    ax.imshow(tl.to_numpy(weight_img), cmap=plt.cm.OrRd, interpolation='nearest')
    ax.set_axis_off()
    if i == 0:
        ax.set_title('Original\nweights')

    for j, rank in enumerate(ranks):

        # Create a tensor Regressor estimator
        estimator = KruskalRegressor(weight_rank=rank, tol=10e-7, n_iter_max=100, reg_W=1, verbose=0)

        # Fit the estimator to the data
        estimator.fit(X, y)

        ax = fig.add_subplot(n_rows, n_columns, i*n_columns + j + 2)
        ax.imshow(tl.to_numpy(estimator.weight_tensor_), cmap=plt.cm.OrRd, interpolation='nearest')
        ax.set_axis_off()

        if i == 0:
            ax.set_title('Learned\nrank = {}'.format(rank))

plt.suptitle("Kruskal tensor regression")
plt.show()

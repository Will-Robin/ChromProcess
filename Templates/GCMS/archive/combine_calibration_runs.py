import sys
sys.path.append('/Users/williamrobinson/documents/nijmegen/packages')
import os
from matplotlib.pyplot as plt
import numpy as np
from ChromProcess import curve_fit
from shutil import copyfile
from scipy.optimize import curve_fit # for the peak fitting process
from pathlib import Path

def quadratic(x, A,B):
    return A*(x**2) + B*x

def fit_upper_lower_bounds(conc_arr, peak_arr):

    popt,pcov = curve_fit(quadratic, conc_arr, peak_arr, p0 = [0.5,0.5])
    perr = np.sqrt(np.diag(pcov))
    p_high = popt + perr
    p_low = popt - perr
    upper = quadratic(conc_arr,*p_high)# upper bound on fit
    lower = quadratic(conc_arr,*p_low) # lower bound on fit
    fit = quadratic(conc_arr,*popt) #  fit

    return lower, fit, upper

calib_folder = Path('/Users/williamrobinson/Documents/Nijmegen/Data/2020_calibrations')
os.chdir(calib_folder)

paths = {}
with open("Calibration_paths.csv", "r") as f:
    next(f)
    for line in f:
        l = line.strip("\n")
        l = l.split(",")

        paths[l[0]] = [x for x in l[1:] if x != ""]

calibration_parameters = {"Sugar":["A","B","C"]}

for p in paths:
    print(p)
    os.chdir(calib_folder/"results/{}".format(p))
    concentrations = []
    peak_1 = []
    peak_2 = []
    readstate = False
    for file in os.listdir():
        if 'Combined_Results' in file: continue
        if file.endswith('.png'): continue
        dataset = []
        with open(file, 'r',encoding='latin1') as f:
            for n,line in enumerate(f,0):
                if "start_data" in line:
                    readstate = True
                    line = next(f)
                if "end_data" in line:
                    readstate = False
                if readstate:
                    newline = line.strip("\n")
                    dataset.append([x for x in newline.split(",") if x != ""])

            e = [list(i) for i in zip(*dataset)]
            d_out = {}
            for s in e:
                d_out[s[0]] = np.array([0 if x == 'nan' else float(x) for x in s[1:]])

            i = np.argsort(d_out['[sample]/[internal standard]'])
            for d in d_out:
                d_out[d] = d_out[d][i]

            concentrations.append(d_out['[sample]/[internal standard]'])

            del d_out['[sample]/[internal standard]']

            highest = 0
            for d in d_out:
                pos = float(d.split('(')[0])
                if pos > highest:
                    highest = pos

            for d in d_out:
                pos = float(d.split('(')[0])
                if pos == highest:
                    peak_2_name = d.split('(')[0]
                    peak_2.append(d_out[d])
                else:
                    peak_1_name = d.split('(')[0]
                    peak_1.append(d_out[d])

    max_length = max([len(x) for x in concentrations])

    conc_arr = np.zeros((len(concentrations),max_length))
    for cn,c in enumerate(concentrations):
        for x in range(0,len(c)):
            conc_arr[cn,x] = c[x]

    peak1_arr = np.zeros((len(concentrations),max_length))
    for cn,c in enumerate(peak_1):
        for x in range(0,len(c)):
            peak1_arr[cn,x] = c[x]

    peak2_arr = np.zeros((len(concentrations),max_length))
    for cn,c in enumerate(peak_2):
        for x in range(0,len(c)):
            peak2_arr[cn,x] = c[x]

    conc_av = np.mean(conc_arr, axis = 0)
    conc_stdev = np.std(conc_arr, axis = 0, ddof = 1)

    peak1_av = np.mean(peak1_arr, axis = 0)
    peak1_stdev = np.std(peak1_arr, axis = 0, ddof = 1)

    peak2_av = np.mean(peak2_arr, axis = 0)
    peak2_stdev = np.std(peak2_arr, axis = 0, ddof = 1)

    conc_stdev = np.nan_to_num(conc_stdev)
    peak1_stdev = np.nan_to_num(peak1_stdev)
    peak2_stdev = np.nan_to_num(peak2_stdev)

    with open('{}_Combined_Results.csv'.format(p), 'w') as f:
        f.write('Average concentration/ M, standard deviation/ M, peak1 integral, peak1 standard deviation,  peak2 integral, peak2 standard deviation')
        f.write('\n')
        for x in range(0,len(conc_av)):
            f.write('{},{},{},{},{},{}\n'.format(conc_av[x], conc_stdev[x], peak1_av[x], peak1_stdev[x], peak2_av[x], peak2_stdev[x]))

    fig, ax = plt.subplots(figsize = (10/2.54,10/2.54))

    '''ax.errorbar(conc_av, peak1_av, xerr = conc_stdev, yerr = peak1_stdev)
    ax.errorbar(conc_av, peak2_av, xerr = conc_stdev, yerr = peak2_stdev)
    '''

    conc_ax = np.ravel(conc_arr)
    i = np.argsort(conc_ax)
    conc_ax = conc_ax[i]
    p1 = np.ravel(peak1_arr)
    p1  = p1[i]
    p2 = np.ravel(peak2_arr)
    p2 = p2[i]

    if np.sum(peak1_arr) != 0:
        lower, fit, upper = fit_upper_lower_bounds(conc_ax, p1)
        ax.fill_between(conc_ax,upper, lower, facecolor = "k", alpha = 0.3)
        ax.plot(conc_ax, fit, c= 'k')
    if np.sum(peak2_arr) != 0:
        lower, fit, upper = fit_upper_lower_bounds(conc_ax, p2)
        ax.fill_between(conc_ax,upper, lower, facecolor = "r", alpha = 0.3)
        ax.plot(conc_ax, fit, c= 'r')

    for x  in range(0,len(conc_arr)):
        if np.sum(peak1_arr) != 0:
            ax.scatter(conc_arr[x],peak1_arr[x], c = 'k', label = peak_1_name)
        if np.sum(peak2_arr) != 0:
            ax.scatter(conc_arr[x],peak2_arr[x], c = 'r', label = peak_2_name)

    plt.legend()
    plt.savefig('{}Calibration_plot.png'.format(p))
    plt.close()

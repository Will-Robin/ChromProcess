import numpy as np
from scipy.optimize import curve_fit
import matplotlib.pyplot as plt
from ChromProcess import info_params
from ChromProcess import calibration_functions as calib_f

def calibration_series(series):
    return CreateCalibrationSeries(series)

def CreateCalibrationSeries(series):

    '''
    Replace sample numbers with concentrations.

    Parameters
    -----------

    Returns
    -------

    '''
    stock_conc = (series.conditions["Mass/ g"][0]/ series.conditions["Mr/ g/mol"][0])/series.conditions["Stock_vol/ L"][0]

    mol_series = [stock_conc*x for x in series.conditions["Titration_series/ L"]]

    series_conc = [(x/series.conditions["Sample_volume/ L"][0])/series.internal_ref_concentration for x in mol_series]

    series.x_series = series_conc

    series.x_name = "[sample]/[internal standard]"

    series.conditions["internal_ref_concentration/ M"] = [series.internal_ref_concentration]

    for c,v in enumerate(series.chromatograms):
        v.timepoint = series.x_series[c]

def fit_calibration_curve_TIC(series):

    '''
    Parameters
    -----------


    Returns
    -------

    '''
    x_vals = np.array(series.x_series)
    calibration_params = {}
    for p in series.integral_series:

        popt,pcov = curve_fit(quadratic, x_vals, series.integral_series[p], p0 = [0.5,0.5,0])
        perr = np.sqrt(np.diag(pcov))
        p_high = popt + perr
        p_low = popt - perr
        upper = calib_f.quadratic(x_vals,*p_high)# upper bound on fit
        lower = calib_f.quadratic(x_vals,*p_low) # lower bound on fit

        calibration_params[p] = popt

    return calibration_params

def fit_calibration_curve_IC(series):

    '''
    Parameters
    -----------


    Returns
    -------

    '''
    x_vals = np.array(series.x_series)
    calibration_params = {}
    for p in series.ion_series:
        calibration_params[p] = {}
        for m in series.ion_series[p]:

            popt,pcov = curve_fit(quadratic, x_vals, series.ion_series[p][m], p0 = [0.5,0.5,0])
            perr = np.sqrt(np.diag(pcov))
            p_high = popt + perr
            p_low = popt - perr
            upper = calib_f.quadratic(x_vals,*p_high)# upper bound on fit
            lower = calib_f.quadratic(x_vals,*p_low) # lower bound on fit

            calibration_params[p][m] = popt

    return calibration_params

def plot_calibration_curves_TIC(series,calib_params):
    fig = plt.figure(figsize = (info_params.across + 4,info_params.up+5))
    ax = plt.subplot(111)
    x_vals = np.array(series.x_series)
    for s in series.integral_series:
        ax.scatter(x_vals, series.integral_series[s], label = round(s,3))
        line = calib_f.quadratic(x_vals,*calib_params[s])

        ax.plot(x_vals, line, c = "k")

    ax.tick_params(axis='both', which='major', labelsize=info_params.labels)
    ax.set_xlabel("{}".format(series.x_name), fontsize = info_params.font)
    ax.set_ylabel("{}".format("I$_S$/I$_I$$_S$"), fontsize = info_params.font)
    # Shrink current axis's height by 10% on the bottom
    box = ax.get_position()
    ax.set_position([box.x0, box.y0 + box.height * 0.3, box.width, box.height * 0.7])
    # Put a legend below current axis
    ax.legend(loc='upper center', bbox_to_anchor=(0.5, -0.15),
              fancybox=True, shadow=True, ncol=3, fontsize = info_params.labels-4)

    plt.savefig("{}_TIC_calibration_lines.png".format(series.set_name))
    plt.clf()
    plt.close()

def plot_calibration_curves_IC(series,calib_params):
    x_vals = np.array(series.x_series)
    for s in series.ion_series:
        fig = plt.figure(figsize = (info_params.across + 4,info_params.up+5))
        ax = plt.subplot(111)

        for m in series.ion_series[s]:
            ax.scatter(x_vals, series.ion_series[s][m], label = m)
            line = calib_f.quadratic(x_vals,*calib_params[s][m])
            ax.plot(x_vals, line, c = "k")

        ax.tick_params(axis='both', which='major', labelsize=info_params.labels)
        ax.set_xlabel("{}".format(series.x_name), fontsize = info_params.font)
        ax.set_ylabel("{}".format("I$_I$$_C$/I$_I$$_S$"), fontsize = info_params.font)
        # Shrink current axis's height by 10% on the bottom
        box = ax.get_position()
        ax.set_position([box.x0, box.y0 + box.height * 0.3, box.width, box.height * 0.7])
        # Put a legend below current axis
        ax.legend(loc='upper center', bbox_to_anchor=(0.5, -0.15),
                  fancybox=True, shadow=True, ncol=3, fontsize = info_params.labels-4)

        plt.savefig("{}_{}_IC_calibration_line.png".format(series.set_name, round(s,3)))
        plt.clf()
        plt.close()

def combine_mz_73(series):

    combined = []
    for p in series.ion_series:
        for m in series.ion_series[p]:
            if int(m) == 73:
                combined.append(series.ion_series[p][m])

    comb = np.array([list(i) for i in zip(*combined)])

    return np.sum(comb, axis = 1)

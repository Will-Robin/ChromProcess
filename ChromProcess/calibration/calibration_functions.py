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

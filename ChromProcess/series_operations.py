import math
import numpy as np
from scipy import interpolate
import matplotlib.pyplot as plt

from ChromProcess import Classes
from ChromProcess import deconvolution
from ChromProcess import processing_functions
from ChromProcess import series_operations as s_o

def cluster(values, bound = 0.1):

    import numpy as np

    cluster = []
    for m in range(0,len(values)):
        '''
        # removed since a cluster need not contain more than one item
        if len(cluster) < 2:
            cluster.append(values[m])
            continue
        '''
        if  abs(values[m]-np.average(cluster)) > bound:
            yield cluster
            cluster = []

        cluster.append(values[m])

    yield cluster


def pick_peaks(series, threshold = 0.01, max_intensity = 1e100, min_intensity = -1e100):
    '''
    For picking peaks in a Chromatogram_Series. Operates on the objects attributes.

    Parameters
    ----------
    series: Chromatogram_Series object
        Object containing chromatograms and associated series data which is
        modified by the function.
    thres : float between [0., 1.]
        Normalized threshold. Only the peaks with amplitude higher than the
        threshold will be detected.
    max_intensity: float
        peaks will not be detected above this threshold.
    min_intensity: float
        peaks will not be detected below this threshold.

    '''
    if len(series.chromatograms) == 0:
        print('No chromatograms in series.')
        exit()
    for c in series.chromatograms:
        for v,b in enumerate(series.regions):

            region_inds = np.where((c.time > b[0]) & (c.time < b[1]))[0]

            peak_inds = processing_functions.Peak_finder(c.signal[region_inds], thres= threshold, min_dist=0.1, max_inten = max_intensity, min_inten = min_intensity)

            for p in range(0,len(peak_inds["Peak_indices"])):

                ret_ind = peak_inds["Peak_indices"][p]
                st_ind  = peak_inds["Peak_start_indices"][p]
                end_ind = peak_inds["Peak_end_indices"][p]

                indices = np.where( (c.time >= c.time[region_inds][st_ind])&(c.time <= c.time[region_inds][end_ind]) )[0]

                if len(indices) > 0:
                    c.peaks[c.time[region_inds][ret_ind]]  = Classes.Peak(c.time[region_inds][ret_ind], indices)

def get_integrals(series):
    '''
    Get integrals of peaks for a Chromatogram_Series
    Parameters
    ----------
    series: Chromatogram_Series object
        Object containing chromatograms and associated series data which is
        modified by the function.

    '''

    for c_count, c in enumerate(series.chromatograms):
        for v,b in enumerate(series.regions):

            region_inds = np.where((c.time > b[0]) & (c.time < b[1]))[0]
            time = c.time[region_inds]
            signal = c.signal[region_inds]

            for p in c.peaks:

                time = c.time[c.peaks[p].indices]
                signal = c.signal[c.peaks[p].indices]

                if c.internal_reference:
                    #linterp = np.interp(time,[time[0],time[-1]],[signal[0],signal[-1]])

                    c.peaks[p].integral = ( np.trapz(signal, x = time) )/c.internal_reference.integral

                else:
                    #linterp = np.interp(time,[time[0],time[-1]],[signal[0],signal[-1]])

                    c.peaks[p].integral = ( np.trapz(signal, x = time) )

def get_internal_ref_integrals(series):
    '''
    Get integrals of internal references for a Chromatogram_Series
    Parameters
    ----------
    series: Chromatogram_Series object
        Object containing chromatograms and associated series data which is
        modified by the function.
    '''
    if series.internal_ref_region:

        for c in series.chromatograms:
            region = series.internal_ref_region

            region_inds = np.where((c.time > region[0]) & (c.time < region[1]))[0]

            time = c.time[region_inds]
            signal = c.signal[region_inds]

            peaks = processing_functions.Peak_finder(signal, thres=0.1, min_dist=0.1, max_inten = 1e100, min_inten = -1e100)
            if len(peaks['Peak_indices']) == 0:
                pass
            else:
                rt    = time[peaks['Peak_indices'][0]]
                start = peaks['Peak_start_indices'][0]
                end   = peaks['Peak_end_indices'][0]

                peak_inds = np.where( (c.time >= c.time[region_inds][start])&(c.time <= c.time[region_inds][end]) )[0]
                c.internal_reference = Classes.Peak(rt, peak_inds)
                #linterp = np.interp(c.time[c.internal_reference.indices],[c.time[c.internal_reference.indices][0],c.time[c.internal_reference.indices][-1]],[c.signal[c.internal_reference.indices][0],c.signal[c.internal_reference.indices][-1]])

                c.internal_reference.integral = ( np.trapz(c.signal[c.internal_reference.indices], x = c.time[c.internal_reference.indices])   )

                series.internal_ref_integrals.append(c.internal_reference.integral)
                series.internal_ref_heights.append(signal[peaks['Peak_indices'][0]])

    else:
        pass

def bin_integral_series(series_log,bound = 0.005):

    '''
    Combine peaks from chromatograms which are close in retention time.

    Parameters
    ----------
    series_log: dict
        dictionary of peak integrals
    bound_dict: dict
        dictionary of peak assignments {"name": [upper,lower ]}
    bound: float
        If two peaks are within this value, they are combined into the same
        integral series.(beware: not perfect so there is a possibility that this
        could combine many unrelated peaks into one series)

    Returns
    -------
    region_log: dict
        A modified dictionary of time series chromatogram integrals divided by
        chromatogram region.
    '''

    ordered_peaks = sorted([*series_log])

    clusters = []
    for c in s_o.cluster(ordered_peaks, bound = bound):
        clusters.append(c)

    out_log = {}
    for c in clusters:
        position = round(np.average(c),3)
        newline = [0 for x in range(0,len(series_log[ordered_peaks[0]]))]

        for o in range(0,len(ordered_peaks)):
            if ordered_peaks[o] in c:
                newline = [x+y for x,y in zip(series_log[ ordered_peaks[o] ], newline) ]

        out_log[position] = newline

    return out_log

def bin_peak_series(series_log,standard_dev = 0.005):

    '''
    Combine peaks from chromatograms which are close in retention time.

    Parameters
    ----------
    series_log: dict
        dictionary of peak integrals
    bound_dict: dict
        dictionary of peak assignments {"name": [upper,lower ]}
    bound: float
        If two peaks are within this value, they are combined into the same
        integral series.(beware: not perfect so there is a possibility that this
        could combine many unrelated peaks into one series)

    Returns
    -------
    region_log: dict
        A modified dictionary of time series chromatogram integrals divided by
        chromatogram region.
    '''

    ordered_peaks = sorted([*series_log])

    clusters = []
    for c in s_o.cluster(ordered_peaks, bound = standard_dev):
        clusters.append(c)

    out_log = {}
    for c in clusters:
        position = round(np.average(c),3)
        newline = [Classes.Peak(0,[0]) for x in range(0,len(series_log[ordered_peaks[0]]))]

        out_log[position] = []

        for o in range(0,len(ordered_peaks)):
            if ordered_peaks[o] in c:
                if len(out_log[position]) == 0:
                    out_log[position] = series_log[ ordered_peaks[o] ]
                else:
                    for count,p in enumerate(out_log[position]):
                        if p.integral == False:
                            out_log[position][count] = series_log[ ordered_peaks[o] ][count]
                        else:
                            pass

    return out_log

def generate_integral_series(series):
    '''
    Generate a time series of integrals from a Chromatogram_Series

    Parameters
    ----------
    series: ChromProcess Chromatogram_Series object
        Object containing chromatograms and associated series data which is
        modified by the function.

    Returns
    -------

    '''

    for count,c in enumerate(series.chromatograms):
        for peak in c.peaks:
            if peak in series.integral_series:
                series.integral_series[peak][count] = c.peaks[peak].integral
            else:
                series.integral_series[peak] = [0 for x in series.x_series]
                series.integral_series[peak][count] = c.peaks[peak].integral

    series.integral_series = s_o.bin_integral_series(series.integral_series, bound = 0.1)

    # get rid of nan values (convert to 0)
    for s in series.integral_series:
        series.integral_series[s] = np.nan_to_num(series.integral_series[s])

def generate_peak_series(series, stdev = 0.05):
    '''
    Generate a time series of integrals from a Chromatogram_Series

    Parameters
    ----------
    series: ChromProcess Chromatogram_Series object
        Object containing chromatograms and associated series data which is
        modified by the function.

    Returns
    -------

    '''

    for count,c in enumerate(series.chromatograms):
        for peak in c.peaks:

            if peak in series.peak_series:
                series.peak_series[peak][count] = c.peaks[peak]
            else:
                series.peak_series[peak] = [Classes.Peak(0,[0]) for x in series.x_series]
                series.peak_series[peak][count] = c.peaks[peak]

    series.peak_series = s_o.bin_peak_series(series.peak_series, standard_dev = stdev)

    for p in series.peak_series:
        series.integral_series[p] = [x.integral for x in series.peak_series[p]]

    # get rid of nan values (convert to 0)
    for s in series.integral_series:
        series.integral_series[s] = np.nan_to_num(series.integral_series[s])

def apply_calibration_to_peak_series(series,calib_dict,bound_dict):
    '''
    Convert integral values into concentrations using supplied calibrations.

    Calibrations come from the calib_dict, which has parameters for a quadratic function:

    Ax^2 + Bx + C

    Parameters
    ----------
    series: ChromProcess Chromatogram_Series object
        Object containing chromatograms and associated series data which is
        modified by the function.

    calib_dict:

    bound_dict: dict
        dictionary of peak assignments {"name": [upper,lower ]}


    Returns
    -------

    '''
    for p in series.peak_series:
        assign = processing_functions.name_peak(float(p),bound_dict)
        if assign in calib_dict:

            if calib_dict[assign]["A"] == 0:
                conversion = lambda x : series.dilution*series.internal_ref_concentration*(x-calib_dict[assign]["C"])/calib_dict[assign]["B"]
            elif calib_dict[assign]["A"] > 0:
                conversion = lambda x : series.dilution*series.internal_ref_concentration*(-calib_dict[assign]["B"] + np.sqrt((calib_dict[assign]["B"]**2) - (4*calib_dict[assign]["A"]*(calib_dict[assign]["C"]-x))))/(2*calib_dict[assign]["A"])
            elif calib_dict[assign]["A"] < 0:
                conversion = lambda x : series.dilution*series.internal_ref_concentration*(-calib_dict[assign]["B"] + np.sqrt((calib_dict[assign]["B"]**2) - (4*calib_dict[assign]["A"]*(calib_dict[assign]["C"]-x))))/(2*calib_dict[assign]["A"])

            concs = np.zeros(len(series.peak_series[p]))
            for x in range(0,len(series.peak_series[p])):
                if series.peak_series[p][x].integral == 0.0:
                    concs[x] = 0.0
                else:
                    concs[x] = conversion(series.peak_series[p][x].integral)

            series.conc_series[p] = concs

        else:
            pass

def bin_dictionary(in_dict, bound = 0.1, rounding = 2):

    sorted_ = sorted([*in_dict])

    clusters = []
    for c in s_o.cluster(sorted_, bound = bound):
        clusters.append(c)

    out_log = {}
    for c in clusters:

        position = round(np.average(c),rounding)

        out_log[position] = []

        for o in range(0,len(sorted_)):
            if sorted_[o] in c:
                if len(out_log[position]) == 0:
                    out_log[position] = in_dict[ sorted_[o] ]
                else:
                    for count,p in enumerate(out_log[position]):
                        if p == 0:
                            out_log[position][count] = in_dict[ sorted_[o] ][count]
                        else:
                            pass

    return out_log

def apply_calibration_to_series(series,calib_dict,bound_dict):
    '''
    Convert integral values into concentrations using supplied calibrations.

    Calibrations come from the calib_dict, which has parameters for a quadratic function:

    Ax^2 + Bx + C

    Parameters
    ----------
    series: ChromProcess Chromatogram_Series object
        Object containing chromatograms and associated series data which is
        modified by the function.

    calib_dict:

    bound_dict: dict
        dictionary of peak assignments {"name": [upper,lower ]}


    Returns
    -------

    '''
    for p in series.integral_series:
        assign = processing_functions.name_peak(float(p),bound_dict)
        if assign in calib_dict:

            if calib_dict[assign]["A"] == 0:
                conversion = lambda x : series.dilution*series.internal_ref_concentration*(x-calib_dict[assign]["C"])/calib_dict[assign]["B"]
            else:
                conversion = lambda x : series.dilution*series.internal_ref_concentration*(-calib_dict[assign]["B"] + np.sqrt((calib_dict[assign]["B"]**2) - (4*calib_dict[assign]["A"]*(calib_dict[assign]["C"]-x))))/(2*calib_dict[assign]["A"])

            concs = np.zeros(len(series.integral_series[p]))
            for x in range(0,len(series.integral_series[p])):
                concs[x] = conversion(series.integral_series[p][x])

            series.conc_series[p] = concs

        else:
            pass

def apply_calibration_to_deconv_series(series,calib_dict,bound_dict):
    '''
    Convert integral values into concentrations using supplied calibrations.

    Calibrations come from the calib_dict, which has parameters for a quadratic function:

    Ax^2 + Bx + C

    Parameters
    ----------
    series: ChromProcess Chromatogram_Series object
        Object containing chromatograms and associated series data which is
        modified by the function.

    calib_dict:

    bound_dict: dict
        dictionary of peak assignments {"name": [upper,lower ]}


    Returns
    -------

    '''
    for p in series.deconvoluted_series:
        assign = processing_functions.name_peak(float(p),bound_dict)
        if assign in calib_dict:

            if calib_dict[assign]["A"] == 0:
                conversion = lambda x : series.dilution*series.internal_ref_concentration*(x-calib_dict[assign]["C"])/calib_dict[assign]["B"]
            else:
                conversion = lambda x : series.dilution*series.internal_ref_concentration*(-calib_dict[assign]["B"] + np.sqrt((calib_dict[assign]["B"]**2) - (4*calib_dict[assign]["A"]*(calib_dict[assign]["C"]-x))))/(2*calib_dict[assign]["A"])

            concs = np.zeros(len(series.deconvoluted_series[p]))
            for x in range(0,len(series.deconvoluted_series[p])):
                concs[x] = conversion(series.deconvoluted_series[p][x])

            series.deconvoluted_series[p] = concs

        else:
            pass

def combine_common_assignments(series, bound_dict):
    '''
    Combines the integral series of signals which have been assigned to the
    same species.

    Parameters
    ----------
    series:

    bound_dict: dict
        dictionary of peak assignments {"name": [upper,lower]}


    Returns
    -------

    '''
    peaks = [*series.conc_series]
    del_list = []
    for p1 in range(0,len(peaks)):
        assign1 = str(processing_functions.name_peak(float(peaks[p1]), bound_dict))
        if p1 in del_list:
            pass
        else:
            for p2 in range(p1,len(peaks)):
                assign2 = str(processing_functions.name_peak(float(peaks[p2]), bound_dict))
                if assign1[:-2] == assign2[:-2] and assign1 != assign2:
                    series.conc_series[peaks[p1]] = [x+y for x,y in zip(series.conc_series[peaks[p1]],series.conc_series[peaks[p2]])]
                    del_list.append(peaks[p2])

    del_list = list(set(del_list))
    for d in del_list:
        del series.conc_series[d]

def base_subtract_series(series):
    for c in series.chromatograms:
        c.signal = processing_functions.baseline_subtraction(c)
    return

def set_integral_series_to_TMS_integral(series, mass_channel = 73):
    '''
    Must be called after ion chromatograms have been integrated
    '''

    # Converting integral series from ion chromatogram to series
    # corresponding to the intergal of m/z = 73 ion chromatograms
    for peak_pos in series.ion_series:
        if mass_channel not in series.ion_series[peak_pos]:
            series.integral_series[peak_pos] = [0 for x in series.integral_series[peak_pos]]
        for i in series.ion_series[peak_pos]:
            if int(i) == mass_channel:
                series.integral_series[peak_pos] = series.ion_series[peak_pos][i]
            else:
                pass

def create_timecourse_matrix(series):
    '''
    Parameters
    -----------
    series: ChromProcess series object

    Returns
    -------
    grid: 2D numpy array
        chromatograms stacked in a grid.
    '''


    # Create a time axis in a region common to all chromatograms
    max_time = 1e100
    min_time = 0
    max_len = 0
    for c in series.chromatograms:
        if c.time[-1] < max_time:
            max_time = c.time[-1]
        if c.time[0] > min_time:
            min_time = c.time[0]
        if len(c.time) > max_len:
            max_len = len(c.time)

    t_axis = np.linspace(min_time,max_time, num = max_len)

    grid = np.empty((len(series.chromatograms),len(t_axis)))

    for c in range(0,len(series.chromatograms)):
        # linear inerpolation of chromatogram trace
        f = interpolate.interp1d(series.chromatograms[c].time,series.chromatograms[c].signal, kind = "linear")
        # use t_axis to create the interpolated signal
        ynew = f(t_axis)
        # add into grid
        grid[c] = ynew

    return t_axis, grid

def create_region_timecourse_matrix(series, region):
    '''
    Parameters
    -----------
    series: ChromProcess series object

    region: list
        [lower bound, upper bound]

    Returns
    -------
    grid: 2D numpy array
        chromatograms stacked in a grid.
    '''

    # Create a time axis in a region common to all chromatograms
    max_time = region[1]
    min_time = region[0]
    max_len = 0
    for c in series.chromatograms:
        if c.time[-1] < max_time:
            max_time = c.time[-1]
        if c.time[0] > min_time:
            min_time = c.time[0]
        if len(c.time) > max_len:
            max_len = len(c.time)

    t_axis = np.linspace(min_time,max_time, num = 500)

    grid = np.empty((len(series.chromatograms),len(t_axis)))

    for c in range(0,len(series.chromatograms)):
        # linear inerpolation of chromatogram trace
        f = interpolate.interp1d(series.chromatograms[c].time,series.chromatograms[c].signal, kind = "linear")
        # use t_axis to create the interpolated signal
        ynew = f(t_axis)
        # add into grid
        grid[c] = ynew

    return t_axis, grid

def correlate_peak(t_pos, grid, time):
    '''
    Parameters
    -----------


    Returns
    -------

    '''

    pos = np.where(np.round(time, 3) == round(t_pos,3))[0][0]
    x_i = grid[:,pos] - grid[:,pos].mean()

    corr = np.zeros(len(grid[0]))
    s_corr = np.zeros(len(grid[0]))

    for n in range(0,len(grid[0])):
        x_j = grid[:,n] - grid[:,n].mean()
        con = np.correlate(x_i,x_j)
        ccor = con/(len(x_i)*x_i.std()*x_j.std())

        corr[n] = con # cross correlation
        s_corr[n] = ccor # scaled cross correlation

    return s_corr

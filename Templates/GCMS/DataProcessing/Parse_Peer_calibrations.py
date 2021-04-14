import os
from pathlib import Path
import numpy as np
from ChromProcess import Classes
from ChromProcess import simple_functions as s_f
from ChromProcess import calibration_operations as cal_ops
from ChromProcess import file_import
from NorthNet import Global_formatting
from NorthNet import info_params
import matplotlib.pyplot as plt

# Get paths to calibration data
database_root = Path(r'C:\Users\willi\Documents\PrebioticDatabase')
exp_paths = Classes.DataPaths(database_root/r'Data_Information\2021_calibration_paths.csv')
storage_stem = database_root/'Data'

calib_output_folder = database_root/r'Analysis_information\GCMS'
calib_dict = {}

# Import calibration data
for e in exp_paths.paths:

    exp_name = e.experiment_code
    print(exp_name)
    exp_folder = storage_stem/e.data_type/'FRN'/exp_name
    fn = [f for f in os.listdir(exp_folder) if 'conditions' in f and not f.startswith('.')][0]
    info_file = exp_folder/fn
    conditions = Classes.Experiment_Conditions(information_file = info_file)

    analysis_file = exp_folder/'{}_analysis_details.csv'.format(exp_name)
    analysis = Classes.Analysis_Information(information_file = analysis_file)

    report_folder = exp_folder/'DataReports'
    report_name = '{}_GCMS_integral_report.csv'.format(exp_name)
    report = Classes.DataReport(file = report_folder/report_name)

    # remove low-concentration components below a bound
    report.remove_entries_below_threshold(0.1)
    '''
    # Get the correct x_axis
    cal_ax = cal_ops.calibration_series_values(conditions, analysis)

    # order everything by increasing series_values
    report.series_values = cal_ax
    report.series_unit = '[sample]/[internal standard]'

    i = np.argsort(report.series_values)
    report.series_values = report.series_values[i]
    for d in report.data:
        report.data[d] = report.data[d][i]
    '''
    # remove internal standard
    for k in [*report.data]:
        if 'CCCCCCCCCCCCCC' in k:
            del report.data[k]
        if k == '7.53':
            del report.data[k]

    calib_dict[exp_name] = report

all_compounds = [x.split('_')[0] for x in calib_dict]
all_compounds = list(set(all_compounds))
groups = {}
for a in all_compounds:
    groups[a] = []
    for c in calib_dict:
        if a in c and len(a) == len(c.split('_')[0]):
            groups[a].append(c)

calib_results = {}
for a,g in zip(all_compounds,groups):

    concentrations = np.array([])
    peak_1 = np.array([])
    peak_2 = np.array([])

    mx_dset_len = 0
    for dr in groups[g]:
        d_set = calib_dict[dr]
        if len(d_set.series_values) > mx_dset_len:
            mx_dset_len = len(d_set.series_values)

    for dr in groups[g]:
        d_set = calib_dict[dr]
        input = np.zeros(mx_dset_len)
        for x in range(0,len(d_set.series_values)):
            input[x] = d_set.series_values[x]
        concentrations = np.hstack((concentrations,input))

    # Choose peak 1 and peak 2 based on size of integral
    peak_1_rts = []
    peak_2_rts = []
    for dr in groups[g]:
        d_set = calib_dict[dr]

        highest = 0
        highest_val = 0
        for d in d_set.data:
            if np.sum(d_set.data[d]) > highest_val:
                highest = d
                highest_val = np.sum(d_set.data[d])

        peak_1_ins = np.zeros(mx_dset_len)
        peak_2_ins = np.zeros(mx_dset_len)
        for d in d_set.data:
            if d == highest:
                peak_1_rts.append(s_f.get_rt_from_header(d))
                for x in range(0,len(d_set.series_values)):
                    peak_1_ins[x] = d_set.data[d][x]
            else:
                peak_2_rts.append(s_f.get_rt_from_header(d))
                for x in range(0,len(d_set.series_values)):
                    peak_2_ins[x] = d_set.data[d][x]

        peak_1_rts = list(set(peak_1_rts))
        peak_2_rts = list(set(peak_2_rts))
        peak_1 = np.hstack((peak_1,peak_1_ins))
        peak_2 = np.hstack((peak_2,peak_2_ins))

    idx = np.argsort(concentrations)
    concentrations = concentrations[idx]
    peak_1 = peak_1[idx]
    peak_2 = peak_2[idx]

    factors1, cov1 = cal_ops.AnalyseCalibrationCurve(concentrations,peak_1,
                                    function = s_f.QuadraticFunctionNoIntercept,
                                    func_p0 = [0.5,0.5],
                                    func_bounds = ((0,0),(5,5)))

    factors2, cov2 = cal_ops.AnalyseCalibrationCurve(concentrations,peak_2,
                                    function = s_f.QuadraticFunctionNoIntercept,
                                    func_p0 = [0.5,0.5],
                                    func_bounds = ((0,0),(5,5)) )

    A1,B1 = factors1[0], factors1[1]
    A1_var, B1_var = cov1[0,0],cov1[1,1]
    A1_B1_var = cov1[0,1]

    A2,B2 = factors2[0], factors2[1]
    A2_var, B2_var = cov2[0,0],cov2[1,1]
    A2_B2_var = cov2[1,0]

    calib_results['{}_1'.format(a)] = {'A':A1, 'B':B1, 'Avar':A1_var, 'Bvar':B1_var, 'ABvar':A1_B1_var, 'RT':peak_1_rts}
    calib_results['{}_2'.format(a)] = {'A':A2, 'B':B2, 'Avar':A2_var, 'Bvar':B2_var, 'ABvar':A2_B2_var, 'RT':peak_2_rts}

    fig, ax = plt.subplots(ncols = 2)
    plt.title(a)
    ax[0].scatter(concentrations,peak_1,c='k')
    ax[1].scatter(concentrations,peak_2,c='r')
    plt.show()
    plt.close()

# write calibration file output
with open(calib_output_folder/'2021_calibrations_analysis.csv', 'w') as f:
    f.write('compound, A, A_variance, B, B_variance, AB_convariance, average_RT/ min\n')
    for c in calib_results:
        f.write('{},'.format(c.split('_')[0]))
        f.write('{},'.format(calib_results[c]['A']))
        f.write('{},'.format(calib_results[c]['Avar']))
        f.write('{},'.format(calib_results[c]['B']))
        f.write('{},'.format(calib_results[c]['Bvar']))
        f.write('{},'.format(calib_results[c]['ABvar']))
        rts  = calib_results[c]['RT']
        if len(rts) > 0:
            print(c, rts)
            f.write('{},'.format(sum(rts)/len(rts)))
        else:
            f.write(',')
        f.write('\n')

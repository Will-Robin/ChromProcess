import sys
sys.path.append('/Users/williamrobinson/documents/nijmegen/packages')
import numpy as np
import matplotlib.pyplot as plt
from ChromProcess.simple_functions import QuadraticFunction
from ChromProcess.calibration_operations import AnalyseCalibrationCurve, CalculateRSE, PredictionSE

def calculate_calibration_error(calc_conc, integral_error,
                                calib_factors,x):

    unew = PredictionSE(calc_conc, integral_error**2,
                    calib_factors['A'], calib_factors['B'], calib_factors['C'],
                    (params['A'][1])**2, (params['B'][1])**2, (params['C'][1])**2,
                    covariances['sab'], covariances['sac'], covariances['sbc'])

    return unew

# Read data
m= np.array([[0.126757475, 0.069048553,	0.031374743],
             [0.315328761, 0.341076431,	0.103120546],
             [0.630657523, 1.02225455,	0.15901761],
             [1.261315045, 2.354055858,	0.412895594],
             [1.891972567, 3.901107043,	0.259454872],
             [2.52263009, 5.103139139,	0.679966042]])
m = m.T
x = m[0]
y = m[1]
y_err = m[2]

# Analyse the calibration curve
ans = AnalyseCalibrationCurve(x,y,y_err)
params = ans['Parameters']
covariances = ans['Covariances']

# Generate additional points on the x-axis for plotting.
xnew = np.linspace(0, 3, num = 20)
ynew = QuadraticFunction(xnew, params['A'][0], params['B'][0], params['C'][0])

# Calculate the error for the predicted value.
unew = PredictionSE(ynew,
                    [np.amax(y_err)**2 for x in ynew],
                    params['A'][0], params['B'][0], params['C'][0],
                    (params['A'][1])**2, (params['B'][1])**2, (params['C'][1])**2,
                    covariances['sab'], covariances['sac'], covariances['sbc'])

fig, ax = plt.subplots()
ax.fill_between(xnew, ynew + unew, y2 = ynew - unew,
                color = 'r', alpha = 0.5,
                label = 'error boundary')
ax.plot(xnew, ynew,
        c = 'r',
        label = 'calibration line')
ax.errorbar(x, y, yerr = y_err,
            marker = 'o', ls = 'none', c = 'k', capsize = 2,
            label = 'data (n = 3)')
ax.set_xlabel('[analyte]\[internal standard]')
ax.set_ylabel('analyte integral/ internal standard integral')
ax.legend()
plt.show()

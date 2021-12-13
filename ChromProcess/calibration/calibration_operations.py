'''
Tools for analysing calibration data.
'''
from ChromProcess import simple_functions as s_f 

def AnalyseCalibrationCurve(x,y,error = [.0,.0],
                            function = s_f.QuadraticFunction,
                            func_p0 = [0.5,0.5,0],
                            func_bounds = ((0,0,0),(2,2,1e-12))):
    '''
    Fits a function to x,y data.

    Parameters
    ----------
    x: array_like
        independent variable for calibration
    y: array_like
        dependent variable for calibration

    Returns
    -------
    results: dict
        dictionary of the results of the analysis.
    '''
    import numpy as np
    from scipy.optimize import curve_fit

    abs_sigma = True
    if np.sum(error) == 0.0:
        error = None
        abs_sigma = False

    # Generate quadratic regression curve
    popt,pcov = curve_fit(function, x, y,
                          p0 = func_p0,
                          bounds = func_bounds,
                          sigma = error, absolute_sigma = abs_sigma)

    return popt, pcov

def CalculateRSE(data, calc):
    '''
    Calculation the residual squared error between two arrays.

    Parameters
    ----------
    data: numpy array
        Data
    calc: numpy array
        Calculated values

    Return
    ------
    rse: float
        residual squared error
    '''
    import numpy as np

    RSS = np.sum(np.square(data - calc))
    rse = np.sqrt(RSS/ (len(data) - 2))
    return rse

def QuadraticPredictionSE(yhat, sy2,
                          a, b, c,
                          sa2, sb2, sc2,
                          sab, sac, sbc):
    '''
    The equation for the quadratic calibration curve is:
    f = sqrt(-b + (b**2 - 4*c*(a-y)))/(2*c)

    Calculation of the standard error of the estimation from a quadratic
    calibration curve.

    Parameters
    ----------
    yhat: float or numpy array
        average of measurement.
    sy2: float or numpy array
        Variance of measurement.
    a, b, c: float
        (second order, first order 0th order terms) average values of calibration
        parameters.
    sa2, sb2, sc2: float
        (second order, first order 0th order terms) variance values of calibration
        parameters.
    sab, sac, sbc: float
        Covariances of calibration parameters.

    Returns
    -------

    unew: float or numpy array
        Calculated standard error of the estimation.

    '''
    import numpy as np

    # The partial derivatives of f with respect to Y is:
    dfdy = 1/np.sqrt(b**2 - 4*a*(c-yhat))

    # The other partial derivatives are:
    b4ac = np.sqrt(b**2 - 4*a*(c-yhat))
    dfda = -1/b4ac
    dfdb = (-1 + b/b4ac)/(2*a)
    dfdc = ( (-4*c + 4*yhat)/(b4ac*(4*a)) - (-b + b4ac)/(2*(a**2)) )

    # The standard deviation of without covariances
    u = np.sqrt(dfdy**2*sy2 + dfda**2*sa2 + dfdb**2*sb2 + dfdc**2*sc2)
    # Adding uncertainty in covariances to standard deviation
    unew = np.sqrt(u**2 + 2*dfda*dfdb*sab + 2*dfda*dfdc*sac + 2*dfdb*dfdc*sbc)

    return unew

def UpperLowerBoundsQuadratic(X, A, B, C):
    '''
    Get the upper and lower bounds for a fit for quadratic calibration.

    Parameters
    ----------
    X: array
        variable for calculation
    A, B, C: tuple
        paramter values as: (average, standard deviation)
    Returns
    -------
    lower, upper: arrays
        Calculation of the lower and upper confidence intervals of the fit.
    '''
    from ChromProcess.simple_functions import QuadraticFunction

    lower = QuadraticFunction(X, A[0] - A[1], B[0] - B[1], C[0] - C[1])
    upper = QuadraticFunction(X, sum(A), sum(B), sum(C))

    return lower, upper

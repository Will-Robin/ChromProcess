import numpy as np

def LinearFunction(X, A, B):
    '''
    Linear function

    Parameters
    ----------
    X: numpy array, float, int
        independent variable
    A, B: float
        Quadratic function parameters.

    Returns
    -------
    Operation on X
    '''
    return A*X + B

def LinearCalibrationFunction(integ, B, C, dil, IS):
    '''
    Solution to linear function
    Parameters
    ----------
    integ: float
        integral value
    B, C: float
        (Gradient, intercept) Calibration parameters
    dil: float
        Dilution factor.
    IS: float
        Internal standard concentration.

    Returns
    -------
    Operation on integ
    '''
    return dil*IS*(integ-C)/B

def QuadraticFunction(X, A, B, C):
    '''
    Quadratic function
    Parameters
    ----------
    X: numpy array, float, int
        independent variable
    A, B: float
        Quadratic function parameters.

    Returns
    -------
    Operation on X
    '''

    return np.nan_to_num(A*(X**2) + B*X + C)

def QuadraticFunctionNoIntercept(X, A, B):
    '''
    Quadratic function
    Parameters
    ----------
    X: numpy array, float, int
        independent variable
    A, B: float
        Quadratic function parameters.

    Returns
    -------
    Operation on X
    '''
    return np.nan_to_num(A*(X**2) + B*X)

def QuadraticCalibrationFunction(integ, A, B, C, dil, IS):
    '''
    Solution to quadratic function
    Parameters
    ----------
    integ: float
        integral value
    A, B, C: float
        (second order, first order 0th order terms) Calibration parameters
    dil: float
        Dilution factor.
    IS: float
        Internal standard concentration.

    Returns
    -------
    Operation on integ
    '''
    return dil*IS*(-B + np.sqrt((B**2) - (4*A*(C-integ))))/(2*A)

def _1gaussian(x, amp1, cen1, sigma1):
    '''
    A single gaussian function
    Parameters
    ----------
    x: array
        x axis data
    amp1: float
        amplitude of the function
    cen1: float
        centre of the function (mean)
    sigma1: float
        width of the function (standard deviation)
    Returns
    -------
    function: numpy array
        y values for the function
    '''
    return amp1*(1/(sigma1*(np.sqrt(2*np.pi))))*(np.exp(-((x-cen1)**2)/((2*sigma1)**2)))

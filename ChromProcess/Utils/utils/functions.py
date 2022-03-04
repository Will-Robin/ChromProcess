import numpy as np


def linear(X, A, B):
    """
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
    """
    return A * X + B


def inverse_linear(variable, gradient, intercept, factor=1.0):
    """
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
    """
    return factor * (variable - intercept) / gradient


def quadratic_function(X, A, B, C):
    """
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
    """
    return np.nan_to_num(A * (X**2) + B * X + C)


def inverse_quadratic(integ, A, B, C, factor=1.0):
    """
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
    """
    return factor * (-B + np.sqrt((B**2) - (4 * A * (C - integ)))) / (2 * A)


def inverse_quadratic_standard_error(yhat, sy2, a, b, c, sa2, sb2, sc2, sab, sac, sbc):
    """
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
    """

    # The partial derivatives of f with respect to Y is:
    dfdy = 1 / np.sqrt(b**2 - 4 * a * (c - yhat))

    # The other partial derivatives are:
    b4ac = np.sqrt(b**2 - 4 * a * (c - yhat))
    dfda = -1 / b4ac
    dfdb = (-1 + b / b4ac) / (2 * a)
    dfdc = (-4 * c + 4 * yhat) / (b4ac * (4 * a)) - (-b + b4ac) / (2 * (a**2))

    # The standard deviation of without covariances
    u = np.sqrt(dfdy**2 * sy2 + dfda**2 * sa2 + dfdb**2 * sb2 + dfdc**2 * sc2)
    # Adding uncertainty in covariances to standard deviation
    unew = np.sqrt(
        u**2 + 2 * dfda * dfdb * sab + 2 * dfda * dfdc * sac + 2 * dfdb * dfdc * sbc
    )

    return unew


def residual_squared_error(data_1, data_2):
    """
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
    """

    RSS = np.sum(np.square(data_1 - data_2))
    rse = np.sqrt(RSS / (len(data_1) - 2))
    return rse


def _1gaussian(x, amp1, cen1, sigma1):
    """
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
    """
    return (
        amp1
        * (1 / (sigma1 * (np.sqrt(2 * np.pi))))
        * (np.exp(-((x - cen1) ** 2) / ((2 * sigma1) ** 2)))
    )

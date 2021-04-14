#import numpy as np
from scipy.optimize import curve_fit # for the peak fitting process
from ChromProcess import deconvolution as d_c
import numpy as np

'''
A bunch of gaussian sums for fitting and the fitting function.
'''

def _1gaussian_grad(x, amp1, cen1, sigma1):
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
    return amp1*(1/(sigma1*(np.sqrt(2*np.pi))))*(np.exp(-((x-cen1)**2)/((2*sigma1)**2))), (x-cen1)*np.exp((-(x-cen1)**2)/2*(sigma1**2))/(np.sqrt(2*np.pi)*sigma1**3)


def autogauss(bounds,N):
    '''
    Any number of exponentials...?
    '''
    '''
    arrays = [[a1 , a2]]*N + [[min(y), max(y)]]*N
    Sort this out into a vectorised function.
    np.dot(amp1,(1/(sigma1*(np.sqrt(2*np.pi))))*(np.exp(-((-cen1+x)**2)/((2*sigma1)**2))))
    '''
    pass

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

def _2gaussian(x, amp1, cen1, sigma1, amp2, cen2, sigma2):
    '''
    A double gaussian function
    Parameters
    ----------
    x: array
        x axis data
    ampn: float
        amplitude of a component gaussian function
    cenn: float
        centre of a component gaussian function (mean)
    sigman: float
        width of a component gaussian function (standard deviation)
    Returns
    -------
    function: numpy array
        y values for the function
    '''
    return d_c._1gaussian(x, amp1, cen1, sigma1) + d_c._1gaussian(x, amp2, cen2, sigma2)

def _3gaussian(x, amp1, cen1, sigma1,amp2, cen2, sigma2, amp3, cen3, sigma3):
    return d_c._1gaussian(x, amp1, cen1, sigma1) + d_c._2gaussian(x, amp2, cen2, sigma2, amp3, cen3, sigma3)
    #return amp1*(1/(sigma1*(np.sqrt(2*np.pi))))*(np.exp(-((x-cen1)**2)/((2*sigma1)**2))) + amp2*(1/(sigma2*(np.sqrt(2*np.pi))))*(np.exp(-((x-cen2)**2)/((2*sigma2)**2))) + amp3*(1/(sigma3*(np.sqrt(2*np.pi))))*(np.exp(-((x-cen3)**2)/((2*sigma3)**2)))

def _4gaussian(x, amp1, cen1, sigma1,amp2, cen2, sigma2, amp3, cen3, sigma3, amp4, cen4, sigma4):
    return d_c._3gaussian(x, amp1, cen1, sigma1,amp2, cen2, sigma2, amp3, cen3, sigma3) + d_c._1gaussian(x, amp4, cen4, sigma4)
    #return amp1*(1/(sigma1*(np.sqrt(2*np.pi))))*(np.exp(-((x-cen1)**2)/((2*sigma1)**2))) + amp2*(1/(sigma2*(np.sqrt(2*np.pi))))*(np.exp(-((x-cen2)**2)/((2*sigma2)**2))) + amp3*(1/(sigma3*(np.sqrt(2*np.pi))))*(np.exp(-((x-cen3)**2)/((2*sigma3)**2))) + amp4*(1/(sigma4*(np.sqrt(2*np.pi))))*(np.exp(-((x-cen4)**2)/((2*sigma4)**2)))

def _5gaussian(x, amp1, cen1, sigma1,amp2, cen2, sigma2, amp3, cen3, sigma3, amp4, cen4, sigma4, amp5, cen5, sigma5):
    return d_c._4gaussian(x, amp1, cen1, sigma1,amp2, cen2, sigma2, amp3, cen3, sigma3, amp4, cen4, sigma4) + d_c._1gaussian(x, amp5, cen5, sigma5)
    #return amp1*(1/(sigma1*(np.sqrt(2*np.pi))))*(np.exp(-((x-cen1)**2)/((2*sigma1)**2))) + amp2*(1/(sigma2*(np.sqrt(2*np.pi))))*(np.exp(-((x-cen2)**2)/((2*sigma2)**2))) + amp3*(1/(sigma3*(np.sqrt(2*np.pi))))*(np.exp(-((x-cen3)**2)/((2*sigma3)**2))) + amp4*(1/(sigma4*(np.sqrt(2*np.pi))))*(np.exp(-((x-cen4)**2)/((2*sigma4)**2))) + amp5*(1/(sigma5*(np.sqrt(2*np.pi))))*(np.exp(-((x-cen5)**2)/((2*sigma5)**2)))

def _6gaussian(x, amp1, cen1, sigma1,amp2, cen2, sigma2, amp3, cen3, sigma3, amp4, cen4, sigma4, amp5, cen5, sigma5, amp6, cen6, sigma6):
    return d_c._5gaussian(x, amp1, cen1, sigma1,amp2, cen2, sigma2, amp3, cen3, sigma3, amp4, cen4, sigma4, amp5, cen5, sigma5) + d_c._1gaussian(x, amp6, cen6, sigma6)
    #return amp1*(1/(sigma1*(np.sqrt(2*np.pi))))*(np.exp(-((x-cen1)**2)/((2*sigma1)**2))) + amp2*(1/(sigma2*(np.sqrt(2*np.pi))))*(np.exp(-((x-cen2)**2)/((2*sigma2)**2))) + amp3*(1/(sigma3*(np.sqrt(2*np.pi))))*(np.exp(-((x-cen3)**2)/((2*sigma3)**2))) + amp4*(1/(sigma4*(np.sqrt(2*np.pi))))*(np.exp(-((x-cen4)**2)/((2*sigma4)**2))) + amp5*(1/(sigma5*(np.sqrt(2*np.pi))))*(np.exp(-((x-cen5)**2)/((2*sigma5)**2))) + amp6*(1/(sigma6*(np.sqrt(2*np.pi))))*(np.exp(-((x-cen6)**2)/((2*sigma6)**2)))

def _7gaussian(x, amp1, cen1, sigma1,amp2, cen2, sigma2, amp3, cen3, sigma3, amp4, cen4, sigma4, amp5, cen5, sigma5, amp6, cen6, sigma6, amp7, cen7, sigma7):
    return d_c._6gaussian(x, amp1, cen1, sigma1,amp2, cen2, sigma2, amp3, cen3, sigma3, amp4, cen4, sigma4, amp5, cen5, sigma5, amp6, cen6, sigma6) + d_c._1gaussian(x, amp7, cen7, sigma7)
    #return amp1*(1/(sigma1*(np.sqrt(2*np.pi))))*(np.exp(-((x-cen1)**2)/((2*sigma1)**2))) + amp2*(1/(sigma2*(np.sqrt(2*np.pi))))*(np.exp(-((x-cen2)**2)/((2*sigma2)**2))) + amp3*(1/(sigma3*(np.sqrt(2*np.pi))))*(np.exp(-((x-cen3)**2)/((2*sigma3)**2))) + amp4*(1/(sigma4*(np.sqrt(2*np.pi))))*(np.exp(-((x-cen4)**2)/((2*sigma4)**2))) + amp5*(1/(sigma5*(np.sqrt(2*np.pi))))*(np.exp(-((x-cen5)**2)/((2*sigma5)**2))) + amp6*(1/(sigma6*(np.sqrt(2*np.pi))))*(np.exp(-((x-cen6)**2)/((2*sigma6)**2))) + amp7*(1/(sigma7*(np.sqrt(2*np.pi))))*(np.exp(-((x-cen7)**2)/((2*sigma7)**2)))

def _8gaussian(x, amp1, cen1, sigma1, amp2, cen2, sigma2, amp3, cen3, sigma3, amp4, cen4, sigma4, amp5, cen5, sigma5, amp6, cen6, sigma6, amp7, cen7, sigma7, amp8, cen8, sigma8):
    return d_c._7gaussian(x, amp1, cen1, sigma1,amp2, cen2, sigma2, amp3, cen3, sigma3, amp4, cen4, sigma4, amp5, cen5, sigma5, amp6, cen6, sigma6, amp7, cen7, sigma7) + d_c._1gaussian(x, amp8, cen8, sigma8)
    #return amp1*(1/(sigma1*(np.sqrt(2*np.pi))))*(np.exp(-((x-cen1)**2)/((2*sigma1)**2))) + amp2*(1/(sigma2*(np.sqrt(2*np.pi))))*(np.exp(-((x-cen2)**2)/((2*sigma2)**2))) + amp3*(1/(sigma3*(np.sqrt(2*np.pi))))*(np.exp(-((x-cen3)**2)/((2*sigma3)**2))) + amp4*(1/(sigma4*(np.sqrt(2*np.pi))))*(np.exp(-((x-cen4)**2)/((2*sigma4)**2))) + amp5*(1/(sigma5*(np.sqrt(2*np.pi))))*(np.exp(-((x-cen5)**2)/((2*sigma5)**2))) + amp6*(1/(sigma6*(np.sqrt(2*np.pi))))*(np.exp(-((x-cen6)**2)/((2*sigma6)**2))) + amp7*(1/(sigma7*(np.sqrt(2*np.pi))))*(np.exp(-((x-cen7)**2)/((2*sigma7)**2))) + amp8*(1/(sigma8*(np.sqrt(2*np.pi))))*(np.exp(-((x-cen8)**2)/((2*sigma8)**2)))

def _9gaussian(x, amp1, cen1, sigma1, amp2, cen2, sigma2, amp3, cen3, sigma3, amp4, cen4, sigma4, amp5, cen5, sigma5, amp6, cen6, sigma6, amp7, cen7, sigma7, amp8, cen8, sigma8,amp9,cen9,sigma9):
    return d_c._8gaussian(x, amp1, cen1, sigma1, amp2, cen2, sigma2, amp3, cen3, sigma3, amp4, cen4, sigma4, amp5, cen5, sigma5, amp6, cen6, sigma6, amp7, cen7, sigma7, amp8, cen8, sigma8) + d_c._1gaussian(x,amp9,cen9,sigma9)
    #return amp1*(1/(sigma1*(np.sqrt(2*np.pi))))*(np.exp(-((x-cen1)**2)/((2*sigma1)**2))) + amp2*(1/(sigma2*(np.sqrt(2*np.pi))))*(np.exp(-((x-cen2)**2)/((2*sigma2)**2))) + amp3*(1/(sigma3*(np.sqrt(2*np.pi))))*(np.exp(-((x-cen3)**2)/((2*sigma3)**2))) + amp4*(1/(sigma4*(np.sqrt(2*np.pi))))*(np.exp(-((x-cen4)**2)/((2*sigma4)**2))) + amp5*(1/(sigma5*(np.sqrt(2*np.pi))))*(np.exp(-((x-cen5)**2)/((2*sigma5)**2))) + amp6*(1/(sigma6*(np.sqrt(2*np.pi))))*(np.exp(-((x-cen6)**2)/((2*sigma6)**2))) + amp7*(1/(sigma7*(np.sqrt(2*np.pi))))*(np.exp(-((x-cen7)**2)/((2*sigma7)**2))) + amp8*(1/(sigma8*(np.sqrt(2*np.pi))))*(np.exp(-((x-cen8)**2)/((2*sigma8)**2))) + amp9*(1/(sigma9*(np.sqrt(2*np.pi))))*(np.exp(-((x-cen9)**2)/((2*sigma9)**2)))

def _10gaussian(x, amp1, cen1, sigma1, amp2, cen2, sigma2, amp3, cen3, sigma3, amp4, cen4, sigma4, amp5, cen5, sigma5, amp6, cen6, sigma6, amp7, cen7, sigma7, amp8, cen8, sigma8,amp9,cen9,sigma9,amp10,cen10,sigma10):
    return d_c._9gaussian(x, amp1, cen1, sigma1, amp2, cen2, sigma2, amp3, cen3, sigma3, amp4, cen4, sigma4, amp5, cen5, sigma5, amp6, cen6, sigma6, amp7, cen7, sigma7, amp8, cen8, sigma8,amp9,cen9,sigma9) + d_c._1gaussian(x, amp10,cen10,sigma10)
    #return amp1*(1/(sigma1*(np.sqrt(2*np.pi))))*(np.exp(-((x-cen1)**2)/((2*sigma1)**2))) + amp2*(1/(sigma2*(np.sqrt(2*np.pi))))*(np.exp(-((x-cen2)**2)/((2*sigma2)**2))) + amp3*(1/(sigma3*(np.sqrt(2*np.pi))))*(np.exp(-((x-cen3)**2)/((2*sigma3)**2))) + amp4*(1/(sigma4*(np.sqrt(2*np.pi))))*(np.exp(-((x-cen4)**2)/((2*sigma4)**2))) + amp5*(1/(sigma5*(np.sqrt(2*np.pi))))*(np.exp(-((x-cen5)**2)/((2*sigma5)**2))) + amp6*(1/(sigma6*(np.sqrt(2*np.pi))))*(np.exp(-((x-cen6)**2)/((2*sigma6)**2))) + amp7*(1/(sigma7*(np.sqrt(2*np.pi))))*(np.exp(-((x-cen7)**2)/((2*sigma7)**2))) + amp8*(1/(sigma8*(np.sqrt(2*np.pi))))*(np.exp(-((x-cen8)**2)/((2*sigma8)**2))) + amp9*(1/(sigma9*(np.sqrt(2*np.pi))))*(np.exp(-((x-cen9)**2)/((2*sigma9)**2))) + amp10*(1/(sigma10*(np.sqrt(2*np.pi))))*(np.exp(-((x-cen10)**2)/((2*sigma10)**2)))

def _11gaussian(x, amp1, cen1, sigma1, amp2, cen2, sigma2, amp3, cen3, sigma3, amp4, cen4, sigma4, amp5, cen5, sigma5, amp6, cen6, sigma6, amp7, cen7, sigma7, amp8, cen8, sigma8,amp9,cen9,sigma9,amp10,cen10,sigma10,amp11,cen11,sigma11):
    return d_c._10gaussian(x, amp1, cen1, sigma1, amp2, cen2, sigma2, amp3, cen3, sigma3, amp4, cen4, sigma4, amp5, cen5, sigma5, amp6, cen6, sigma6, amp7, cen7, sigma7, amp8, cen8, sigma8,amp9,cen9,sigma9,amp10,cen10,sigma10) + d_c._1gaussian(x, amp11,cen11,sigma11)
    #return amp1*(1/(sigma1*(np.sqrt(2*np.pi))))*(np.exp(-((x-cen1)**2)/((2*sigma1)**2))) + amp2*(1/(sigma2*(np.sqrt(2*np.pi))))*(np.exp(-((x-cen2)**2)/((2*sigma2)**2))) + amp3*(1/(sigma3*(np.sqrt(2*np.pi))))*(np.exp(-((x-cen3)**2)/((2*sigma3)**2))) + amp4*(1/(sigma4*(np.sqrt(2*np.pi))))*(np.exp(-((x-cen4)**2)/((2*sigma4)**2))) + amp5*(1/(sigma5*(np.sqrt(2*np.pi))))*(np.exp(-((x-cen5)**2)/((2*sigma5)**2))) + amp6*(1/(sigma6*(np.sqrt(2*np.pi))))*(np.exp(-((x-cen6)**2)/((2*sigma6)**2))) + amp7*(1/(sigma7*(np.sqrt(2*np.pi))))*(np.exp(-((x-cen7)**2)/((2*sigma7)**2))) + amp8*(1/(sigma8*(np.sqrt(2*np.pi))))*(np.exp(-((x-cen8)**2)/((2*sigma8)**2))) + amp9*(1/(sigma9*(np.sqrt(2*np.pi))))*(np.exp(-((x-cen9)**2)/((2*sigma9)**2))) + amp10*(1/(sigma10*(np.sqrt(2*np.pi))))*(np.exp(-((x-cen10)**2)/((2*sigma10)**2))) + amp11*(1/(sigma11*(np.sqrt(2*np.pi))))*(np.exp(-((x-cen11)**2)/((2*sigma11)**2)))

def _12gaussian(x, amp1, cen1, sigma1, amp2, cen2, sigma2, amp3, cen3, sigma3, amp4, cen4, sigma4, amp5, cen5, sigma5, amp6, cen6, sigma6, amp7, cen7, sigma7, amp8, cen8, sigma8,amp9,cen9,sigma9,amp10,cen10,sigma10,amp11,cen11,sigma11,amp12,cen12,sigma12):
    return d_c._11gaussian(x, amp1, cen1, sigma1, amp2, cen2, sigma2, amp3, cen3, sigma3, amp4, cen4, sigma4, amp5, cen5, sigma5, amp6, cen6, sigma6, amp7, cen7, sigma7, amp8, cen8, sigma8,amp9,cen9,sigma9,amp10,cen10,sigma10,amp11,cen11,sigma11) + d_c._1gaussian(x, amp12,cen12,sigma12)

def _13gaussian(x, amp1, cen1, sigma1, amp2, cen2, sigma2, amp3, cen3, sigma3, amp4, cen4, sigma4, amp5, cen5, sigma5, amp6, cen6, sigma6, amp7, cen7, sigma7, amp8, cen8, sigma8,amp9,cen9,sigma9,amp10,cen10,sigma10,amp11,cen11,sigma11,amp12,cen12,sigma12,amp13,cen13,sigma13):
    return d_c._12gaussian(x, amp1, cen1, sigma1, amp2, cen2, sigma2, amp3, cen3, sigma3, amp4, cen4, sigma4, amp5, cen5, sigma5, amp6, cen6, sigma6, amp7, cen7, sigma7, amp8, cen8, sigma8,amp9,cen9,sigma9,amp10,cen10,sigma10,amp11,cen11,sigma11,amp12,cen12,sigma12) + d_c._1gaussian(x, amp13,cen13,sigma13)

def _14gaussian(x, amp1, cen1, sigma1, amp2, cen2, sigma2, amp3, cen3, sigma3, amp4, cen4, sigma4, amp5, cen5, sigma5, amp6, cen6, sigma6, amp7, cen7, sigma7, amp8, cen8, sigma8,amp9,cen9,sigma9,amp10,cen10,sigma10,amp11,cen11,sigma11,amp12,cen12,sigma12,amp13,cen13,sigma13,amp14,cen14,sigma14):
    return d_c._13gaussian(x, amp1, cen1, sigma1, amp2, cen2, sigma2, amp3, cen3, sigma3, amp4, cen4, sigma4, amp5, cen5, sigma5, amp6, cen6, sigma6, amp7, cen7, sigma7, amp8, cen8, sigma8,amp9,cen9,sigma9,amp10,cen10,sigma10,amp11,cen11,sigma11,amp12,cen12,sigma12,amp13,cen13,sigma13) + d_c._1gaussian(x, amp14,cen14,sigma14)

def _15gaussian(x, amp1, cen1, sigma1, amp2, cen2, sigma2, amp3, cen3, sigma3, amp4, cen4, sigma4, amp5, cen5, sigma5, amp6, cen6, sigma6, amp7, cen7, sigma7, amp8, cen8, sigma8,amp9,cen9,sigma9,amp10,cen10,sigma10,amp11,cen11,sigma11,amp12,cen12,sigma12,amp13,cen13,sigma13,amp14,cen14,sigma14,amp15,cen15,sigma15):
    return d_c._14gaussian(x, amp1, cen1, sigma1, amp2, cen2, sigma2, amp3, cen3, sigma3, amp4, cen4, sigma4, amp5, cen5, sigma5, amp6, cen6, sigma6, amp7, cen7, sigma7, amp8, cen8, sigma8,amp9,cen9,sigma9,amp10,cen10,sigma10,amp11,cen11,sigma11,amp12,cen12,sigma12,amp13,cen13,sigma13,amp14,cen14,sigma14) + d_c._1gaussian(x, amp15,cen15,sigma15)

def _16gaussian(x, amp1, cen1, sigma1, amp2, cen2, sigma2, amp3, cen3, sigma3, amp4, cen4, sigma4, amp5, cen5, sigma5, amp6, cen6, sigma6, amp7, cen7, sigma7, amp8, cen8, sigma8,amp9,cen9,sigma9,amp10,cen10,sigma10,amp11,cen11,sigma11,amp12,cen12,sigma12,amp13,cen13,sigma13,amp14,cen14,sigma14,amp15,cen15,sigma15,amp16,cen16,sigma16):
    return d_c._15gaussian(x, amp1, cen1, sigma1, amp2, cen2, sigma2, amp3, cen3, sigma3, amp4, cen4, sigma4, amp5, cen5, sigma5, amp6, cen6, sigma6, amp7, cen7, sigma7, amp8, cen8, sigma8,amp9,cen9,sigma9,amp10,cen10,sigma10,amp11,cen11,sigma11,amp12,cen12,sigma12,amp13,cen13,sigma13,amp14,cen14,sigma14,amp15,cen15,sigma15) + d_c._1gaussian(x, amp16,cen16,sigma16)

def _17gaussian(x, amp1, cen1, sigma1, amp2, cen2, sigma2, amp3, cen3, sigma3, amp4, cen4, sigma4, amp5, cen5, sigma5, amp6, cen6, sigma6, amp7, cen7, sigma7, amp8, cen8, sigma8,amp9,cen9,sigma9,amp10,cen10,sigma10,amp11,cen11,sigma11,amp12,cen12,sigma12,amp13,cen13,sigma13,amp14,cen14,sigma14,amp15,cen15,sigma15,amp16,cen16,sigma16,amp17,cen17,sigma17):
    return d_c._16gaussian(x, amp1, cen1, sigma1, amp2, cen2, sigma2, amp3, cen3, sigma3, amp4, cen4, sigma4, amp5, cen5, sigma5, amp6, cen6, sigma6, amp7, cen7, sigma7, amp8, cen8, sigma8,amp9,cen9,sigma9,amp10,cen10,sigma10,amp11,cen11,sigma11,amp12,cen12,sigma12,amp13,cen13,sigma13,amp14,cen14,sigma14,amp15,cen15,sigma15,amp16,cen16,sigma16) + d_c._1gaussian(x, amp17,cen17,sigma17)

def _18gaussian(x, amp1, cen1, sigma1, amp2, cen2, sigma2, amp3, cen3, sigma3, amp4, cen4, sigma4, amp5, cen5, sigma5, amp6, cen6, sigma6, amp7, cen7, sigma7, amp8, cen8, sigma8,amp9,cen9,sigma9,amp10,cen10,sigma10,amp11,cen11,sigma11,amp12,cen12,sigma12,amp13,cen13,sigma13,amp14,cen14,sigma14,amp15,cen15,sigma15,amp16,cen16,sigma16,amp17,cen17,sigma17,amp18,cen18,sigma18):
    return d_c._17gaussian(x, amp1, cen1, sigma1, amp2, cen2, sigma2, amp3, cen3, sigma3, amp4, cen4, sigma4, amp5, cen5, sigma5, amp6, cen6, sigma6, amp7, cen7, sigma7, amp8, cen8, sigma8,amp9,cen9,sigma9,amp10,cen10,sigma10,amp11,cen11,sigma11,amp12,cen12,sigma12,amp13,cen13,sigma13,amp14,cen14,sigma14,amp15,cen15,sigma15,amp16,cen16,sigma16,amp17,cen17,sigma17) + d_c._1gaussian(x, amp18,cen18,sigma18)

def _19gaussian(x, amp1, cen1, sigma1, amp2, cen2, sigma2, amp3, cen3, sigma3, amp4, cen4, sigma4, amp5, cen5, sigma5, amp6, cen6, sigma6, amp7, cen7, sigma7, amp8, cen8, sigma8,amp9,cen9,sigma9,amp10,cen10,sigma10,amp11,cen11,sigma11,amp12,cen12,sigma12,amp13,cen13,sigma13,amp14,cen14,sigma14,amp15,cen15,sigma15,amp16,cen16,sigma16,amp17,cen17,sigma17,amp18,cen18,sigma18,amp19,cen19,sigma19):
    return d_c._18gaussian(x, amp1, cen1, sigma1, amp2, cen2, sigma2, amp3, cen3, sigma3, amp4, cen4, sigma4, amp5, cen5, sigma5, amp6, cen6, sigma6, amp7, cen7, sigma7, amp8, cen8, sigma8,amp9,cen9,sigma9,amp10,cen10,sigma10,amp11,cen11,sigma11,amp12,cen12,sigma12,amp13,cen13,sigma13,amp14,cen14,sigma14,amp15,cen15,sigma15,amp16,cen16,sigma16,amp17,cen17,sigma17,amp18,cen18,sigma18) + d_c._1gaussian(x, amp19,cen19,sigma19)

def fit_gaussian_peaks(time, sig, peaks, d_type = 'null'):
    '''
    Fitting sums of gaussian peaks to data using supplied peak indices.``
    Parameters
    ----------
    time: array
        time values
    sig: array
        signal to be fit to
    peaks: list of peak positions
        list peak positions in the time
    d_type: str
        HPLC or GCMS. Used to determine which sets of parameters to use.
    Returns
    -------
    popt: ndarray
        list of fitted values [[amplitude, centre, width],]
    pcov: array
        correlation matrix
    '''

    guess = [] # amp, cen, sig
    lbds = []
    ubds = []

    if d_type == 'HPLC':
        initial_guess = [1000                , 1, 0.01]
        lowerbounds   = [0                   , 0,   0.0]
        upperbounds   = [np.amax(sig)*1.1    , 0,  0.15]

    elif d_type == 'GCMS':
        initial_guess = [10000              , 1, 0.005]
        lowerbounds   = [0                  , 0,  0.0]
        upperbounds   = [np.amax(sig)*1.1, 1,  0.025]

    else:
        print('Gaussian function fit: No (or unrecognised) chromatogram type provided. Using default initialisations.')
        initial_guess = [1000                , 1, 0.001]
        lowerbounds   = [0                   , 0,   0.0]
        upperbounds   = [np.amax(sig)*1.1    , 0,  0.025]

    for p in range(0,len(peaks)): # extend the boundary
        initial_guess[0] = np.amax(sig)
        initial_guess[1] = peaks[p]
        lowerbounds[1] = time[0]
        upperbounds[1] = time[-1]
        guess.extend(initial_guess)
        lbds.extend(lowerbounds)
        ubds.extend(upperbounds)

    boundarr  = [lbds,ubds]
    try:
        if len(peaks) == 1:
            popt, pcov = curve_fit(d_c._1gaussian, time, sig, p0=guess, bounds = boundarr)
        elif len(peaks) == 2:
            popt, pcov = curve_fit(d_c._2gaussian, time, sig, p0=guess, bounds = boundarr)
        elif len(peaks) == 3:
            popt, pcov = curve_fit(d_c._3gaussian, time, sig, p0=guess, bounds = boundarr)
        elif len(peaks) == 4:
            popt, pcov = curve_fit(d_c._4gaussian, time, sig, p0=guess, bounds = boundarr)
        elif len(peaks) == 5:
            popt, pcov = curve_fit(d_c._5gaussian, time, sig, p0=guess, bounds = boundarr)
        elif len(peaks) == 6:
            popt, pcov = curve_fit(d_c._6gaussian, time, sig, p0=guess, bounds = boundarr)
        elif len(peaks) == 7:
            popt, pcov = curve_fit(d_c._7gaussian, time, sig, p0=guess, bounds = boundarr)
        elif len(peaks) == 8:
            popt, pcov = curve_fit(d_c._8gaussian, time, sig, p0=guess, bounds = boundarr)
        elif len(peaks) == 9:
            popt, pcov = curve_fit(d_c._9gaussian, time, sig, p0=guess, bounds = boundarr)
        elif len(peaks) == 10:
            popt, pcov = curve_fit(d_c._10gaussian, time, sig, p0=guess, bounds = boundarr)
        elif len(peaks) == 11:
            popt, pcov = curve_fit(d_c._11gaussian, time, sig, p0=guess, bounds = boundarr)
        elif len(peaks) == 12:
            popt, pcov = curve_fit(d_c._12gaussian, time, sig, p0=guess, bounds = boundarr)
        elif len(peaks) == 13:
            popt, pcov = curve_fit(d_c._13gaussian, time, sig, p0=guess, bounds = boundarr)
        elif len(peaks) == 14:
            popt, pcov = curve_fit(d_c._14gaussian, time, sig, p0=guess, bounds = boundarr)
        elif len(peaks) == 15:
            popt, pcov = curve_fit(d_c._15gaussian, time, sig, p0=guess, bounds = boundarr)
        elif len(peaks) == 16:
            popt, pcov = curve_fit(d_c._16gaussian, time, sig, p0=guess, bounds = boundarr)
        elif len(peaks) == 17:
            popt, pcov = curve_fit(d_c._17gaussian, time, sig, p0=guess, bounds = boundarr)
        elif len(peaks) == 18:
            popt, pcov = curve_fit(d_c._18gaussian, time, sig, p0=guess, bounds = boundarr)
        elif len(peaks) == 19:
            popt, pcov = curve_fit(d_c._19gaussian, time, sig, p0=guess, bounds = boundarr)
        else:
            print('error', len(peak_indices))
            popt, pcov = [0,0,0],0
    except:
        popt, pcov = [0,0,0],0

    return popt, pcov

def deconvolute_region(chromatogram, region, num_peaks = 1):

    '''
    chromatogram: ChromProcess Chromatogram object
        Chromatogram
    region: list
        region of chromatogram under operation [lower bound, upper bound]
    '''

    upper = region[1]
    lower = region[0]

    inds = np.where((chromatogram.time > lower)&(chromatogram.time < upper))[0]

    time = chromatogram.time[inds]
    signal = chromatogram.signal[inds]

    signal = signal-np.average(signal[-5:-1])

    peak_list = np.array([*chromatogram.peaks])
    peak_inds = np.where((peak_list > lower)&(peak_list < upper))[0]

    peaks = peak_list[peak_inds]

    while len(peaks) < num_peaks:
        peaks = np.append(peaks, np.average(peaks))

    if len(peaks) > num_peaks:
        peaks = peaks[:num_peaks]


    return d_c.fit_gaussian_peaks(time, signal, peaks, d_type = 'GCMS')

def deconvolute_peak(peak, chromatogram, num_peaks = 2):

    '''
    chromatogram: ChromProcess Chromatogram object
        Chromatogram
    region: list
        region of chromatogram under operation [lower bound, upper bound]
    '''

    peaks = [peak.retention_time for p in range(num_peaks)]
    time = chromatogram.time[peak.indices]
    signal = chromatogram.signal[peak.indices]
    baseline = np.interp(time, [time[0],time[-1]], [signal[0],signal[-1]])

    signal = signal-baseline

    return d_c.fit_gaussian_peaks(time, signal, peaks, d_type = 'GCMS')

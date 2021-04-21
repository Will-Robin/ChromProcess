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
    import numpy as np
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
    import numpy as np
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
    import numpy as np

    return dil*IS*(-B + np.sqrt((B**2) - (4*A*(C-x))))/(2*A)

def _1gaussian(x, amp1, cen1, sigma1):
    import numpy as np
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

def runSVD(matrix):
    ''''
    Parameters
    ----------
    matrix: array

    Returns
    -------
    U: numpy array
        Unitary matrix ith left singular vectors

    S: numpy array
        array of singular values

    Vh: numpy array
        Unitary matrix with right singular vectors
    '''

    import numpy as np
    from ChromProcess import simple_functions as simp_func

    # perform SVD using numpy
    U, S, Vh = np.linalg.svd(matrix, full_matrices = False,
                             compute_uv = True,
                             hermitian = False)

    # Flip signs so that largest vals are positive based on U
    U, Vh = simp_func.SVDflipSigns(U,Vh)

    return U, S, Vh

def SVDflipSigns(U,V):
    '''
    Adapted from sklearn code
    https://github.com/scikit-learn/

    use columns of U, rows of V

    Flip signs so that largest vals are positive based on U

    '''
    import numpy as np

    max_abs_cols = np.argmax(np.abs(U), axis=0)
    signs = np.sign(U[max_abs_cols, range(U.shape[1])])

    U = U*signs
    V = V*signs[:, np.newaxis]

    return U, V

def runPCA(matrix, n_components = 10):

    from ChromProcess import simple_functions as simp_func
    samples = matrix.shape[0]

    U, S, Vh = simp_func.runSVD(matrix)


    exp_variance = (S**2)/(samples - 1)
    exp_var_ratio = exp_variance/exp_variance.sum()

    results = {'LH_singular_vectors':U[:,:n_components],
               'Singular_values':S[:n_components],
               'RH_singular_vectors':Vh[:n_components],
               'Explained_variance':exp_variance[:n_components],
               'Explained_variance_ratio':exp_var_ratio[:n_components]}

    return results

def reconstitute_full_SVD(U,S,Vh):
    import numpy as np
    return np.dot(U*S,Vh)

def SVD_transform(U,S,Vh,components = 1):
    import numpy as np
    Vh = Vh[:components]
    U = U[:, :components]
    U *= S[:components]
    return np.dot(U,Vh)

def cluster(values, bound = 0.1):

    import numpy as np

    values = np.sort(values)

    cluster = []
    for m in range(0,len(values)):
        if len(cluster) == 0:
            pass
        else:
            clust_av = np.average(cluster)

            if  abs(values[m]-clust_av) > bound:
                yield cluster
                cluster = []

        cluster.append(values[m])

    yield cluster

def cluster_indices(values, bound = 0.1):

    import numpy as np

    values = np.sort(values)

    cluster = []
    for m in range(0,len(values)):
        if len(cluster) == 0:
            pass
        else:
            clust_av = np.average(values[cluster])

            if  abs(values[m]-clust_av) > bound:
                yield cluster
                cluster = []

        cluster.append(m)

    yield cluster

def isfloat(thing):
    '''
    Test if an thing (e.g. str) can be converted to a float.
    '''
    try:
        float(thing)
        return True
    except ValueError:
        return False

def get_rt_from_header(element):
    from ChromProcess import simple_functions as s_f

    if s_f.isfloat(element):
        position = float(element)
    else:
        spl = element.split('(')[-1].strip(')')
        position = float(spl[0:-1])
    return position

def upper_tri_no_diag(arr):
    m = arr.shape[0]
    r,c = np.triu_indices(m,1)
    return arr[r,c]

def sum_error_prop(errors):
    import numpy as np
    sq_err = [x**2 for x in errors]

    error = np.sqrt(sum(sq_err))

    return error

def mult_div_error_prop(averages, errors):
    import numpy as np

    sq_err = [(b/a)**2 for a,b in zip(averages, errors)]
    error = np.sqrt(sum(sq_err))

    return error

def weighted_stdev(samples, weights):
    '''
    Parameters
    ----------
    samples, weights: numpy arrays (1D)

    Returns
    -------
    numpy array (1D)
    '''
    import numpy as np

    averages = np.mean(samples)
    x_x_bar = (samples - averages)**2
    wt_av_2 = weights*x_x_bar
    sum = np.sum(wt_av_2)

    denom = ((len(weights)-1)/ len(weights))*np.sum(weights)

    return sum/denom

def weighted_cov(samples, weights, ddof = 1):

    import numpy as np

    v1 = np.sum(weights)
    w_samples = np.zeros(samples.shape)
    for x in range(0,len(samples)):
        w_samples[x] = samples[x]*weights[x]

    cov = np.dot(w_samples, w_samples.T) * v1 / (v1**2 - ddof * v1)

    return cov

def weighted_corrcoef(samples,weights):
    from ChromProcess import simple_functions as s_f
    import numpy as np

    cov = s_f.weighted_cov(samples, weights, ddof = 1)
    v = np.sqrt(np.diag(cov))
    outer_v = np.outer(v, v)
    correlation = cov / outer_v
    correlation[cov == 0] = 0

    return correlation


def sparse_corrcoef(A):

    A = A.astype(np.float64)
    n = A.shape[1]

    # Compute the covariance matrix
    rowsum = A.sum(1)
    centering = rowsum.dot(rowsum.T.conjugate()) / n
    C = (A.dot(A.T.conjugate()) - centering) / (n - 1)

    # The correlation coefficients are given by
    # C_{i,j} / sqrt(C_{i} * C_{j})
    d = np.diag(C)
    coeffs = C / np.sqrt(np.outer(d, d))

    return coeffs

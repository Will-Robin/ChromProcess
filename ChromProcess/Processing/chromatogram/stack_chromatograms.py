import numpy as np

def stack_chromatograms(chromatogram_list):
    '''
    A bit crude. Could also stack chromatograms using interpolations.

    Parameters
    ----------
    chromatogram_list: list of ChromProcess Chromatogram objects.
        List of chromatograms to stack.

    Returns
    -------
    chrom_stack: numpy array (2D)
        Stack of chromatograms
        shape = (len(chromatogram_list), minimum chromatogram time axis length)
    '''
    import numpy as np

    min_length = 1e100
    for c in chromatogram_list:
        if len(c.time) < min_length:
            min_length = len(c.time)

    chrom_stack = np.empty((len(chromatogram_list),min_length))

    for c in range(0,len(chromatogram_list)):
        chrom_stack[c] = chromatogram_list[c].signal[:min_length]

    return chrom_stack

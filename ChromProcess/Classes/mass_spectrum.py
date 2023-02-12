class MassSpectrum:
    """
    Container for mass spectral data.
    """

    def __init__(self, mz, inten, pos=None):
        """
        Parameters
        ----------
        mz: numpy.ndarray[np.float64]
        inten: numpy.ndarray[np.float64]
        pos: None or float

        Returns
        -------
        None
        """

        self.mz = mz
        self.relative_abundances = inten
        self.retention_time = pos

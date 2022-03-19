from ChromProcess.Writers import mass_spectrum_to_string_rows


class MassSpectrum:
    def __init__(self, mz, inten, pos=None):
        """
        Parameters
        ----------
        fname: str or pathlib Path
        mz: array
        inten: array
        pos: None or float

        Returns
        -------
        None
        """

        self.mz = mz
        self.relative_abundances = inten
        self.retention_time = pos

    def write_to_string(self):
        """
        Write a row-wise csv for the mass spectrum.

        Parameters
        ----------

        Returns
        -------
        ms_string: str
        """

        ms_string = mass_spectrum_to_string_rows(self)

        return ms_string

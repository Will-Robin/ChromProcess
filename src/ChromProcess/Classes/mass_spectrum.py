import numpy as np


class MassSpectrum:
    """
    Container for mass spectral data.
    """

    def __init__(
        self,
        mz: np.ndarray,
        inten: np.ndarray,
        pos: float | None = None,
    ):
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

    def __str__(self):
        str_repr = f"""
        retention time: {self.retention_time}
        m/z values: {len(self.mz)}
        relative abundance values: {len(self.relative_abundances)}
        """
        return str_repr

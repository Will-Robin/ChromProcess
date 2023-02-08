import numpy as np


class Peak:
    """
    Container for peak information.
    """

    def __init__(
        self,
        retention_time: float,
        start: float,
        end: float,
        indices: list[int] = [],
        integral: float = None,
        height: float = None,
        parent: str = "",
        mass_spectrum: list[np.ndarray] = None,
    ):
        """
        Creates a Peak object using a retention time
        and the indices of the places in the data's
        parent chromatogram from which the time and
        signal of the peak can be obtained.

        Parameters
        ----------
        retention_time: float
        start: float
        end: float
        indices: list
        integral: float
        parent: str
        mass_spectrum: bool or MassSpectrum object

        Attributes
        ----------
        self.retention_time: float
        self.start: float
        self.end: float
        self.indices: array or list
        self.integral: float
        self.height: float
        self.ion_chromatograms: dict
        self.ion_integrals: dict
        self.mass_spectrum: bool or MassSpectrum object.
        self.deconvolution: list
        self.assignment: str
        self.concentration: bool or float
        self.conc_error: bool or float
        self.parent_peak_collection: str
        """

        self.retention_time = retention_time

        self.start = start

        self.end = end

        self.indices = indices

        if integral == None:
            self.integral = 0.0
        else:
            self.integral = integral

        if height == None:
            self.height = 0.0
        else:
            self.height = height

        self.ion_chromatograms = dict()

        self.ion_integrals = dict()

        self.deconvolution = []

        self.assignment: str = None

        self.concentration: float = None

        self.conc_error: float = None

        self.mass_spectrum: list[np.ndarray] = mass_spectrum

        self.parent_peak_collection: str = parent

    def set_retention_time(self, retention_time: float):
        """
        Set the peak retention time.

        Parameters
        ----------
        retention_time: float
        """
        self.retention_time = retention_time

    def set_integral(self, integral: float):
        """
        Set Integral.

        Parameters
        ----------
        integral: float
        """

        self.integral = integral

    def set_height(self, height: float):
        """
        Set the height of the peak.

        Parameters
        ----------
        height: float

        Returns
        -------

        """
        self.height = height

    def set_concentration(self, concentration: float):
        """
        Set the peak concentration.

        Parameters
        ----------
        concentration: float

        Returns
        -------

        """
        self.concentration = concentration

    def set_conc_error(self, error: float):
        """
        Set the concentration error of the peak.

        Parameters
        ----------
        error: float

        Returns
        -------

        """
        self.conc_error = error

    def set_assignment(self, assignment: str):
        """
        Get the assignment of a peak.

        Parameters
        ----------
        assignment: str

        Returns
        -------
        None
        """
        self.assignment = assignment

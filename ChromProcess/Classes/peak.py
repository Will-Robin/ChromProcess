import numpy as np
from typing import Optional
from .mass_spectrum import MassSpectrum
from .deconvolution_parameters import Deconvolution


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
        integral: float = 0.0,
        height: float = 0.0,
        mass_spectrum: MassSpectrum | None = None,
        deconvolution_params: Deconvolution | None = None,
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
        """

        self.retention_time = retention_time

        self.start: float = start

        self.end: float = end

        self.indices: list[int] = indices

        self.integral: float = 0.0
        if integral is not None:
            self.integral = integral

        self.height: float = 0.0
        if height is not None:
            self.height = height

        self.ion_chromatograms: dict[float, np.ndarray] = dict()

        self.ion_integrals: dict[float, float] = dict()

        self.deconvolution: list = []

        self.assignment: str = ""

        self.concentration: float = 0.0

        self.conc_error: float = 0.0

        self.mass_spectrum: Optional[MassSpectrum] = mass_spectrum

        self.deconvolution_params = deconvolution_params

    def __str__(self):
        str_repr = f"""
        retention time: {self.retention_time}
        start: {self.start}
        end: {self.end}
        height: {self.height}
        integral: {self.integral}
        ion chromatograms: {len(self.ion_chromatograms)}
        assignment: {self.assignment}
        """
        return str_repr

    def set_retention_time(self, retention_time: float) -> None:
        """
        Set the peak retention time.

        Parameters
        ----------
        retention_time: float
        """
        self.retention_time = retention_time

    def set_integral(self, integral: float) -> None:
        """
        Set Integral.

        Parameters
        ----------
        integral: float
        """

        self.integral = integral

    def set_height(self, height: float) -> None:
        """
        Set the height of the peak.

        Parameters
        ----------
        height: float

        Returns
        -------

        """
        self.height = height

    def set_concentration(self, concentration: float) -> None:
        """
        Set the peak concentration.

        Parameters
        ----------
        concentration: float

        Returns
        -------

        """
        self.concentration = concentration

    def set_conc_error(self, error: float) -> None:
        """
        Set the concentration error of the peak.

        Parameters
        ----------
        error: float

        Returns
        -------

        """
        self.conc_error = error

    def set_assignment(self, assignment: str) -> None:
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

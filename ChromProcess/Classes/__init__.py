"""
These Classes are Python objects designed to encapsulate data. They provide an
intuitive wrapper around pure Python and numpy code for interaction with the
various kinds of information required for data analysis.

The objects also have some methods for writing to files, conversions and
getting data. Often, they call code from another submodule of ChromProcess, and
can be considered as shorthand for common operations.
"""

from .peak import Peak
from .chromatogram import Chromatogram
from .mass_spectrum import MassSpectrum

from .instrument_calibration import InstrumentCalibration
from .analysis_information import AnalysisInformation
from .experiment_conditions import ExperimentConditions

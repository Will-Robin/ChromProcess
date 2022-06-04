"""
These Classes are Python objects designed to encapsulate data. They provide an
intuitive wrapper around pure Python and numpy code for interaction with the
various kinds of information required for data analysis.

The objects also have some methods for writing to files, conversions and
getting data. Often, they call code from another submodule of ChromProcess, and
can be considered as shorthand for common operations.
"""

from ChromProcess.Classes.Peak import Peak
from ChromProcess.Classes.Chromatogram import Chromatogram
from ChromProcess.Classes.MassSpectrum import MassSpectrum

from ChromProcess.Classes.PeakCollection import PeakCollection
from ChromProcess.Classes.PeakCollectionSeries import PeakCollectionSeries

from ChromProcess.Classes.InstrumentCalibration import InstrumentCalibration
from ChromProcess.Classes.AnalysisInformation import AnalysisInformation
from ChromProcess.Classes.ExperimentConditions import ExperimentConditions

from ChromProcess.Classes.DataReport import DataReport

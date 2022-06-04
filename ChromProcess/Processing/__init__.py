"""
This submodule contains functions which operate on data contained in
ChromProcess.Classes objects.

Most functions do not modify the objects, but some for example,
.chromatogram.modify_chromatogram.py, do perform in place modifications.
"""

from .chromatogram.find_peaks import find_peaks_in_region
from .chromatogram.background_subtraction import ic_background_subtraction

from .chromatogram.modify_chromatogram import integrate_chromatogram_peaks
from .chromatogram.modify_chromatogram import add_peaks_to_chromatogram
from .chromatogram.modify_chromatogram import internal_standard_integral

from .chromatogram.stack_chromatograms import stack_chromatograms

from .peak.assign_peak import assign_retention_time

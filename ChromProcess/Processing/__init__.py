"""
This submodule contains functions which operate on data contained in
ChromProcess.Classes objects.

Most functions do not modify the objects, but some for example,
.chromatogram.modify_chromatogram.py, do perform in place modifications.
"""

from .chromatogram.find_peaks import find_peaks_in_region
from .chromatogram.deconvolution import deconvolute_region
from .chromatogram.background_subtraction import ic_background_subtraction
from .chromatogram.background_subtraction import background_subtraction

from .chromatogram.stack_chromatograms import stack_chromatograms

from .peak.assign_peak import assign_retention_time

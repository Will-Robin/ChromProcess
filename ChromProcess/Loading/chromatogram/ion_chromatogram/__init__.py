"""
Functions for generating ion chromatograms from GC-MS data.
"""
from .ion_chromatogram_from_peak import ion_chromatogram_from_peak
from .ion_chromatogram_from_region import ion_chromatogram_from_region

__all__ = ["ion_chromatogram_from_peak", "ion_chromatogram_from_region"]

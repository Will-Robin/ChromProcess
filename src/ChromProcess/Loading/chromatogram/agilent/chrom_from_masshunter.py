try:
    import rainbow as rb
except ImportError as err:
    _rainbow_import_error = err
    rb = None
else:
    _rainbow_import_error = None

from typing import Union
from pathlib import Path
from ChromProcess.Classes import Chromatogram
from ChromProcess.Loading.chromatogram.adapters import chrom_from_rainbow


def chrom_from_masshunter(filename: Union[str, Path]) -> Chromatogram:
    """
    Load a chromatogram from an Agilent MassHunter .D directory.

    Parameters
    ----------
    filename: str or pathlib.Path

    Returns
    -------
    Chromatogram
    """

    if rb is None:
        raise ImportError(
            "Loading MassHunter files requires the optional dependency "
            "'rainbow-api'. Install it with:\n\n"
            "pip install ChromProcess[masshunter]"
        ) from _rainbow_import_error

    data = rb.read(filename)

    return chrom_from_rainbow(data)

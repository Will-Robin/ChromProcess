from ChromProcess.Classes import Peak
from ChromProcess.Classes import Chromatogram
from ChromProcess.Classes import MassSpectrum


def mass_spectrum_from_peak(peak: Peak, chromatogram: Chromatogram) -> MassSpectrum:
    """
    Get the mass spectrum at the apex of a peak.

    Parameters
    ----------
    peak: Peak
        Peak to find mass spectrum for.
    chromatogram: Chromatogram
        Parent chromatogram for the peak

    Returns
    -------
    mass_spectrum: MassSpectrum
    """

    mz, inten = chromatogram.get_mass_spectrum(peak.retention_time)

    valid = inten > 0

    spectrum = MassSpectrum(
        mz[valid],
        inten[valid],
        pos=peak.retention_time,
    )

    peak.mass_spectrum = spectrum

    return spectrum

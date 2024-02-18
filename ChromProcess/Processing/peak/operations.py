import numpy as np

from ChromProcess.Classes import Peak
from ChromProcess.Classes import Chromatogram

from ChromProcess.Utils.utils.functions import inverse_quadratic
from ChromProcess.Utils.utils.functions import inverse_quadratic_standard_error
from ChromProcess.Utils.utils.functions import inverse_linear

from ChromProcess.Processing.peak.assign_peak import assign_retention_time
from ChromProcess.Utils.utils.error_propagation import mult_div_error_prop


def get_peak_mass_spectrum(peak: Peak, chromatogram: Chromatogram):
    """
    Get the mass spectrum at the apex of the peak.

    Parameters
    ----------
    peak: Peak
    chromatogram: Chromatogram
        Parent chromatogram of the peak.

    Return
    ------
    mass_spectrum: list[np.ndarray]
        Mass spectrum.
    """

    mass_spectrum = []
    if len(chromatogram.scan_indices) > 0:
        time = chromatogram.time

        ind = np.where(time == peak.retention_time)[0]

        start = chromatogram.scan_indices[ind][0]
        end = start + chromatogram.point_counts[ind][0]

        mass_spectrum = [
            np.round(chromatogram.mz_values[start:end], 2),
            chromatogram.mz_intensity[start:end],
        ]

    return mass_spectrum


def get_peak_integral(
    peak: Peak, chromatogram: Chromatogram, baseline_subtract: bool = False
) -> float:
    """
    Get the integral of the peak using its parent chromatogram. Note
    that an arbitray chromatogram can be passed to this method, meaning
    it is not secure. The baseline substraction substracts a baseliner
    interpolated linearly between the start and the end of the peak.

    Parameters
    ----------
    peak: Peak
    chromatogram: Chromatogram
    baseline_subtract: bool

    Returns
    -------
    integral: float
        Integral of the peak.
    """

    time = chromatogram.time[peak.indices]
    signal = chromatogram.signal[peak.indices]

    if baseline_subtract:
        time_bound = [time[0], time[-1]]
        signal_bound = [signal[0], signal[-1]]
        linterp = np.interp(time, time_bound, signal_bound)
        integral = np.trapz(signal - linterp, x=time)
    else:
        integral = np.trapz(signal, x=time)

    return integral


def get_peak_height(
    peak: Peak, chromatogram: Chromatogram, baseline_subtract: bool = False
) -> float:
    """
    Get the height of the peak.

    Parameters
    ----------
    peak: Peak
    chromatogram: Chromatogram
    baseline_subtract: bool

    Returns
    -------
    height: float
        Height of the peak.
    """

    if baseline_subtract:
        time = chromatogram.time[peak.indices]
        signal = chromatogram.signal[peak.indices]

        idx = np.abs(time - peak.retention_time).argmin()
        time_bound = [time[0], time[-1]]
        signal_bound = [signal[0], signal[-1]]

        linterp = np.interp(time, time_bound, signal_bound)
        baseline_at_peak = linterp[idx]
        height = signal[idx] - baseline_at_peak
    else:
        idx = np.abs(chromatogram.time - peak.retention_time).argmin()
        height = chromatogram.signal[idx]

    return height


def apply_linear_calibration(
    peak: Peak, A: float, B: float, internal_standard: float = 1.0
):
    """
    Apply a linear calibration conversion to the peak integral to obtain a
    concentration value.

    $integral/IS integral = A*concentration + B$

    Parameters
    ----------
    peak: Peak
    A: float
    B: float
    internal_standard: float

    Returns
    -------
    concentration: float
    """

    c1 = inverse_linear(peak.integral, A, B)

    concentration = internal_standard * c1

    return concentration


def apply_quadratic_calibration(
    peak: Peak, A: float, B: float, C: float, internal_standard: float = 1.0
):
    """
    Apply a quadratic calibration conversion to the peak integral to obtain
    a concentration value.

    $integral/IS integral = A*concentration^2 + B*concentration + C$

    Parameters
    ----------
    peak: Peak
    A, B, C, internal_standard: float

    Returns
    -------
    concentration: float
    """

    c1 = inverse_quadratic(peak.integral, A, B, C)

    concentration = internal_standard * c1

    if np.isnan(peak.concentration):
        concentration = apply_linear_calibration(
            peak, B, C, internal_standard=internal_standard
        )

    return concentration


def calculate_concentration_error(
    peak: Peak,
    A: float,
    B: float,
    C: float,
    A_variance: float,
    B_variance: float,
    C_variance: float,
    AB_covariance: float,
    AC_covariance: float,
    BC_covariance: float,
    IS_conc: float,
    IS_conc_err: float,
) -> float:
    """
    Calculation of the standard error on a concentration estimation from
    the calibration.

    $integral/IS integral = A*concentration^2 + B*concentration + C$

    Modifies Peak object attributes.

    Parameters
    ----------
    peak: Peak
    A: float,
    B: float,
    C: float,
    A_variance: float,
    B_variance: float,
    C_variance: float,
    AB_covariance: float,
    AC_covariance: float,
    BC_covariance: float,
    IS_conc: float
        Concentration of the internal standard.
    IS_conc_err: float
        Concentration error for the internal standard.

    Returns
    -------
    error: float
    """

    yhat = peak.integral
    sy2 = 1e-10  # Assume small error on integral

    err = inverse_quadratic_standard_error(
        yhat,
        sy2,
        A,
        B,
        C,
        A_variance,
        B_variance,
        C_variance,
        AB_covariance,
        AC_covariance,
        BC_covariance,
    )
    err = np.nan_to_num(err)
    val = inverse_quadratic(yhat, A, B, C)
    err = IS_conc * val * mult_div_error_prop([val, IS_conc], [err, IS_conc_err])

    error = np.nan_to_num(err)

    return error


def dilution_correction(
    peak: Peak, factor: float, factor_error: float
) -> tuple[float, float]:
    """
    Apply a correction to obtain the sample concentration considering its
    dilution before analysis.

    Parameters
    ----------
    peak: Peak
    factor: float
        Factor by which the concentration value must be multiplied to
        obtain the sample concentration before dilution.
    factor_error: float
        Error for the dilution factor.

    Returns
    -------
    (corr_conc, err) : (float, float)
    """

    err_prop = mult_div_error_prop(
        [peak.concentration, factor], [peak.conc_error, factor_error]
    )

    corr_conc: float = peak.concentration * factor
    err: float = err_prop * corr_conc

    return (corr_conc, err)


def get_peak_assignment(peak: Peak, boundaries: dict[str, list[float]]) -> str:
    """
    Assign a name to the peak using boundaries.

    Parameters
    ----------
    peak: Peak
    boundaries: dict
        {'compound name': [lower bound, upper bound]}

    Returns
    -------
    str
    """

    return assign_retention_time(peak, boundaries)

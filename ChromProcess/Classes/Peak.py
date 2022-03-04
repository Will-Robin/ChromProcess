import numpy as np
from ChromProcess.Utils.utils import functions
from ChromProcess.Processing.peak import assign_peak
from ChromProcess.Utils.utils import error_propagation as error_prop


class Peak:
    def __init__(
        self,
        retention_time,
        start,
        end,
        indices=[],
        integral=None,
        parent="",
        mass_spectrum=False,
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

        self.height = 0.0

        self.ion_chromatograms = {}

        self.ion_integrals = {}

        self.deconvolution = []

        self.assignment = "Unknown"

        self.concentration = False

        self.conc_error = False

        self.mass_spectrum = mass_spectrum

        self.parent_peak_collection = parent

    def get_integral(self, chromatogram, baseline_subtract=False):
        """
        Get the integral of the peak using a chromatogram. Note that an
        arbitray chromatogram can be passed to this method, meaning it is not
        secure. The baseline substraction substracts a baseliner interpolated
        linearly between the start and the end of the peak.

        Parameters
        ----------
        chromatogram: ChromProcess Chromatogram object
        baseline_subtract: bool

        Returns
        -------
        self.integral: float
            Integral of the peak.
        """

        time = chromatogram.time[self.indices]
        signal = chromatogram.signal[self.indices]

        if baseline_subtract:
            time_bound = [time[0], time[-1]]
            signal_bound = [signal[0], signal[-1]]
            linterp = np.interp(time, time_bound, signal_bound)
            self.integral = np.trapz(signal - linterp, x=time)
        else:
            self.integral = np.trapz(signal, x=time)

        return self.integral

    def get_height(self, chromatogram):
        """
        Get the height of the peak.

        Parameters
        ----------
        chromatogram: ChromProcess Chromatogram object

        Returns
        -------
        self.height: float
            Height of the peak.
        """

        self.height = chromatogram.signal[self.indices]

        return self.height

    def get_mass_spectrum(self, chromatogram):
        """
        Get the mass spectrum at the apex of the peak.

        Parameters
        ----------
        chromatogram: ChromProcess Chromatogram object
            Parent chromatogram of the peak.

        Return
        ------
        self.mass_spectrum
            Mass spectrum.
        """

        if len(chromatogram.scan_indices) == 0:
            # no mass spectral information in the chromatogram
            pass
        else:
            time = chromatogram.time

            ind = np.where(time == self.retention_time)[0]

            start = chromatogram.scan_indices[ind][0]
            end = start + chromatogram.point_counts[ind][0]

            self.mass_spectrum = [
                np.round(chromatogram.mass_values[start:end], 2),
                chromatogram.mass_intensity[start:end],
            ]

        return self.mass_spectrum

    def reference_integral_to_IS(self, IS_integral):
        """
        Divide the peak integral by IS_integral.

        Parameters
        ----------
        IS_integral: float

        Returns
        -------
        None
        """

        if IS_integral > 0.0:
            self.integral = self.integral / IS_integral
        else:
            pass

    def assign_peak(self, boundaries):
        """
        Assign a name to the peak using boundaries.

        Parameters
        ----------
        boundaries: dict
            {'compound name': [lower bound, upper bound]}

        Returns
        -------
        None
        """

        self.assignment = assign_peak.assign_retention_time(
            self.retention_time, boundaries
        )

    def apply_linear_calibration(self, A, B, internal_standard=1.0):
        """
        Apply a linear calibration conversion to the peak integral to obtain a
        concentration value.

        Parameters
        ----------
        A, B, internal_standard: float

        Returns
        -------
        None
        """

        c1 = functions.inverse_linear(self.integral, A, B, factor=1.0)
        self.concentration = internal_standard * c1

    def apply_quadratic_calibration(self, A, B, C, internal_standard=1.0):
        """
        Apply a quadratic calibration conversion to the peak integral to obtain
        a concentration value.

        Parameters
        ----------
        A, B, C, internal_standard: float

        Returns
        -------
        None
        """

        c1 = functions.inverse_quadratic(self.integral, A, B, C, factor=1.0)

        self.concentration = internal_standard * c1

        if np.isnan(self.concentration):
            self.apply_linear_calibration(B, C, internal_standard=internal_standard)

    def calculate_error(self, calibrations, IS_conc, IS_conc_err):
        """
        Calculation of the standard error on a concentration estimation from
        the calibration.

        Modifies Peak object attributes.

        Parameters
        ----------
        calibrations: ChromProcess Instrument_Calibration object
            Contains calibration information.
        IS_conc: float
            Concentration of the internal standard.
        IS_conc_err: float
            Concentration error for the internal standard.

        Returns
        -------
        None
        """

        assign = self.assignment
        yhat = self.integral
        sy2 = 1e-10

        if assign in calibrations.calibration_factors:

            a = calibrations.calibration_factors[assign]["A"]
            b = calibrations.calibration_factors[assign]["B"]
            c = calibrations.calibration_factors[assign]["C"]

            sa2 = calibrations.calibration_factors[assign]["A_variance"]
            sb2 = calibrations.calibration_factors[assign]["B_variance"]
            sc2 = calibrations.calibration_factors[assign]["C_variance"]

            sab = calibrations.calibration_factors[assign]["AB_covariance"]
            sac = calibrations.calibration_factors[assign]["AC_covariance"]
            sbc = calibrations.calibration_factors[assign]["BC_covariance"]

            err = functions.inverse_quadratic_standard_error(
                yhat, sy2, a, b, c, sa2, sb2, sc2, sab, sac, sbc
            )
            err = np.nan_to_num(err)
            val = functions.inverse_quadratic(yhat, a, b, c, factor=1.0)
            err = (
                IS_conc
                * val
                * error_prop.mult_div_error_prop([val, IS_conc], [err, IS_conc_err])
            )

            self.conc_error = np.nan_to_num(err)

    def dilution_correction(self, factor, factor_error):
        """
        Apply a correction to obtain the sample concentration considering its
        dilution before analysis.

        Parameters
        ----------
        factor: float
            Factor by which the concentration value must be multiplied to
            obtain the sample concentration before dilution.
        factor_error: float
            Error for the dilution factor.

        Returns
        -------
        None
        """

        err = error_prop.mult_div_error_prop(
            [self.concentration, factor], [self.conc_error, factor_error]
        )

        corr_conc = self.concentration * factor
        err *= corr_conc
        self.concentration = corr_conc
        self.conc_error = err

import numpy as np
from ChromProcess.Utils.utils import functions
from ChromProcess.Processing.peak import assign_peak
from ChromProcess.Utils.utils import error_propagation as error_prop


class PeakCollectionElement:
    """
    Information on individual peaks
    """

    def __init__(
        self,
        position,
        integral,
        start,
        end,
        parent="not specified",
        mass_spectrum=False,
    ):
        """
        An object for storage and manipulation of peak data.

        Parameters
        ----------
        self.retention_time: float
        self.start: float
        self.end: float
        self.integral: float
        self.assignment: str
        self.concentration: bool or float
        self.conc_error: bool or float
        self.mass_spectrum: bool or MassSpectrum object.
        self.parent_peak_collection: str
        """

        self.retention_time = position
        self.start = start
        self.end = end
        self.integral = integral
        self.assignment = "Unknown"
        self.concentration = False
        self.conc_error = False
        self.mass_spectrum = mass_spectrum
        self.parent_peak_collection = parent

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

        Modifies PeakCollectionElement object attributes.

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

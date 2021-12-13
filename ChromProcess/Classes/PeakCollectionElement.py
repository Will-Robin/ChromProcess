
class PeakCollectionElement:
    '''
    Information on individual peaks
    '''
    def __init__(self, position, integral, start, end, parent = 'not specified',
                    mass_spectrum = False):
        '''
        Parameters
        ----------
        position, integral, start, end: float
        parent: str
        mass_spectrum: bool
        '''

        self.retention_time = position
        self.start = start
        self.end = end
        self.integral = integral
        self.assignment = 'Unknown'
        self.concentration = False
        self.conc_error = False
        self.mass_spectrum = mass_spectrum
        self.parent_peak_collection = parent

    def inspect_peak(self):
        print('retention_time',self.retention_time)
        print('start',self.start)
        print('end',self.end)
        print('integral',self.integral)
        print()

    def reference_integral_to_IS(self, IS_integral):
        '''
        Parameters
        ----------
        IS_integral: float
        '''
        if IS_integral > 0.0:
            self.integral = self.integral/IS_integral
        else:
            pass

    def assign_peak(self, boundaries):
        '''
        Parameters
        ----------
        boundaries: dict
            {'compound name': [lower bound, upper bound]}
        '''
        from ChromProcess import processing_functions as p_f
        self.assignment = p_f.name_peak(self.retention_time,boundaries)

    def apply_linear_calibration(self, A, B, internal_standard = 1.0):
        '''
        Parameters
        ----------
        A, B, internal_standard: float
        '''

        conversion = lambda x : (x-B)/A
        c1 = conversion(self.integral)
        self.concentration = internal_standard*c1

    def apply_quadratic_calibration(self, A, B, C, internal_standard = 1.0):
        '''
        Parameters
        ----------
        A, B, C, internal_standard: float
        '''
        import numpy as np

        conversion = lambda x : ((-B+np.sqrt((B**2)-(4*A*(C-x))))/(2*A))
        c1 = conversion(self.integral)

        self.concentration = internal_standard*c1
        if np.isnan(self.concentration):
            self.apply_linear_calibration(B,C,
                                          internal_standard = internal_standard)

    def calculate_error(self,calibrations,IS_conc,IS_conc_err):
        '''
        Calculation of the standard error on a concentration estimation from
        th calibration.

        Parameters
        ----------
        calibrations: ChromProcess Instrument_Calibration object
            Contains calibration information.
        Returns
        -------
        None

        Modifies PeakCollectionElement object attributes.
        '''
        import numpy as np

        from ChromProcess import calibration as cal_ops
        from ChromProcess import simple_functions as s_f

        if self.assignment in calibrations.calibration_factors:
            yhat = self.integral
            sy2 = 1e-10

            a = calibrations.calibration_factors[self.assignment]['A']
            b = calibrations.calibration_factors[self.assignment]['B']
            c = calibrations.calibration_factors[self.assignment]['C']

            sa2 = calibrations.calibration_factors[self.assignment]['A_variance']
            sb2 = calibrations.calibration_factors[self.assignment]['B_variance']
            sc2 = calibrations.calibration_factors[self.assignment]['C_variance']

            sab = calibrations.calibration_factors[self.assignment]['AB_covariance']
            sac = calibrations.calibration_factors[self.assignment]['AC_covariance']
            sbc = calibrations.calibration_factors[self.assignment]['BC_covariance']

            err = cal_ops.QuadraticPredictionSE(yhat, sy2,
                                                 a, b, c,
                                                 sa2, sb2, sc2,
                                                 sab, sac, sbc)
            err = np.nan_to_num(err)
            val = ((-b+np.sqrt((b**2)-(4*a*(c-self.integral))))/(2*a))

            err = IS_conc*val*s_f.mult_div_error_prop([val, IS_conc],
                                          [err, IS_conc_err])

            self.conc_error = np.nan_to_num(err)

    def dilution_correction(self,factor,factor_error):
        '''
        Parameters
        ----------
        factor: float
        factor_error: float
        '''
        from ChromProcess import simple_functions as s_f
        err=s_f.mult_div_error_prop([self.concentration, factor],
                                    [self.conc_error, factor_error])

        corr_conc = self.concentration*factor
        err*=corr_conc
        self.concentration = corr_conc
        self.conc_error = err


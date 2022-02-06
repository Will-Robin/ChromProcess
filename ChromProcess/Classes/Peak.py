import numpy as np

class Peak:
    def __init__(self,retention_time, indices):
        '''
        Creates a Peak object using a retention time
        and the indices of the places in the data's
        parent chromatogram from which the time and
        signal of the peak can be obtained.

        Parameters
        ----------
        retention_time: float
        indices: list

        Attributes
        ----------
        self.retention_time: float
        self.indices: array or list
        self.integral: float
        self.height: float
        self.ion_chromatograms: dict
        self.ion_integrals: dict
        self.mass_spectrum: list
        self.deconvolution: list
        '''

        self.retention_time = retention_time

        self.indices = indices

        self.integral = 0.0

        self.height = 0.0

        self.ion_chromatograms = {}

        self.ion_integrals = {}

        self.mass_spectrum = []

        self.deconvolution = []

    def get_integral(self, chromatogram, baseline_subtract = False):
        '''
        Get the integral of the peak using a chromatogram. Note that an
        arbitray chromatogram can be passed to this method, meaning it is not
        secure. The baseline substraction substracts a baseliner interpolated
        linearly between the start and the end of the peak.

        Parameters
        ----------
        chromatogram: ChromProcess Chromatogram object
        baseline_subtract: bool

        Returns:
        self.integral: float
            Integral of the peak.
        '''

        time = chromatogram.time[self.indices]
        signal = chromatogram.signal[self.indices]

        if baseline_subtract:
            time_bound = [time[0], time[-1]]
            signal_bound = [signal[0], signal[-1]]
            linterp = np.interp(time, time_bound, signal_bound)
            self.integral = ( np.trapz(signal - linterp, x = time) )
        else:
            self.integral = ( np.trapz(signal, x = time) )

        return self.integral

    def get_height(self, chromatogram):
        '''
        Get the height of the peak.

        Parameters
        ----------
        chromatogram: ChromProcess Chromatogram object

        Returns
        -------
        self.height: float
            Height of the peak.
        '''

        self.height = chromatogram.signal[self.indices]

        return self.height

    def get_mass_spectrum(self, chromatogram):
        '''
        Get the mass spectrum at the apex of the peak.

        Parameters
        ----------
        chromatogram: ChromProcess Chromatogram object
            Parent chromatogram of the peak.

        Return
        ------
        self.mass_spectrum
            Mass spectrum.
        '''

        if len(chromatogram.scan_indices) == 0:
            # no mass spectral information in the chromatogram
            pass
        else:
            time = chromatogram.time

            ind = np.where(time == self.retention_time)[0]

            start = chromatogram.scan_indices[ind][0]
            end = start + chromatogram.point_counts[ind][0]

            self.mass_spectrum = [
                                np.round(chromatogram.mass_values[start:end],2),
                                chromatogram.mass_intensity[start:end]
                                ]

        return self.mass_spectrum



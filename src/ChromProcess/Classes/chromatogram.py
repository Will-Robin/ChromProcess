import numpy as np

from .peak import Peak


class Chromatogram:
    """
    A class for storing chromatographic data.
    """

    def __init__(self) -> None:
        """
        Initialise an empty chromatogram.

        Attributes
        ----------
        self.filename: str
            Name of the file containing the source data.
        self.x_unit: str
            Name of the x units (prefer 'unit name/ SI symbol')
        self.y_unit: str
            Name of the y units (prefer 'unit name/ SI symbol')
        self.time: numpy.ndarray[np.float64]
            Time axis for the chromatogram.
        self.signal: numpy.ndarray[np.float64]
            Signal of the chromatogram.
        self.peaks: dict
            Container for the peaks derived from the chromatogram.
        self.deconvoluted_peaks: dict
            Container for the peaks deconvoluted from the chromatogram.
        self.mz_values: numpy.ndarray[np.int]
            Container for the m/z values from mass spectra.
        self.mz_intensity: numpy.ndarray[np.int]
            Container for the intensity values from mass spectra.
        self.scan_indices: numpy.ndarray[np.int]
            Container for the scan indices from mass spectra.
        self.point_counts: numpy.ndarray[np.int]
            Container for the point counts from mass spectra.
        self.internal_standard: Peak
            The internal standard peak.
        """

        self.filename: str = ""
        self.x_unit: str = ""
        self.y_unit: str = ""

        self.time: np.ndarray = np.array([])
        self.signal: np.ndarray = np.array([])

        self.peaks: dict[float, Peak] = dict()
        self.deconvoluted_peaks: dict[float, Peak] = dict()

        self.mz_values: np.ndarray = np.array([])
        self.mz_intensity: np.ndarray = np.array([])
        self.scan_indices: np.ndarray = np.array([])
        self.point_counts: np.ndarray = np.array([])

        self.internal_standard: Peak = Peak(0.0, 0.0, 0.0)

    def __str__(self) -> str:
        str_repr = f"""
        file name: {self.filename}
        x unit: {self.x_unit}
        y unit: {self.y_unit}
        time: {len(self.time)} data points
        signal: {len(self.signal)} data points
        peaks: {len(self.peaks)}
        deconvoluted peaks: {len(self.deconvoluted_peaks)}
        internal_standard: {self.internal_standard.retention_time}
        mass spectral infomation: {"yes" if len(self.mz_values) > 0 else "no"}
        """
        return str_repr

    def add_peaks(self, peaks: list[Peak]) -> None:
        """
        Add peaks to a chromatogram.

        Parameters
        ----------
        peaks: list[Peak]

        Returns
        ------
        None
        """

        for peak in peaks:
            rt = peak.retention_time
            indices = np.where((self.time >= peak.start) & (self.time <= peak.end))[0]
            peak.indices = indices.tolist()

            if peak.deconvolution_params:
                self.deconvoluted_peaks[rt] = peak
            else:
                self.peaks[rt] = peak

    def set_internal_standard(self, peak: Peak) -> None:
        """
        Set the internal standard of the chromatogram.

        Parameters
        ----------
        peak: Peak

        Returns
        -------

        """
        self.internal_standard = peak

    def get_mass_spectrum(self, time: float) -> tuple[np.ndarray, np.ndarray]:
        """
        Get the mass spectrum at a given time point in the chromatogram.

        Parameters
        ----------
        time: float
            Time point for the mass spectrum in the chromatogram.

        Returns
        -------
        m_z: numpy.ndarray[np.float64]
            m/z value
        intensity: numpy.ndarray[np.float64]
            Ion counts.
        """

        inds = np.where(self.time == time)[0]

        scan_inds = self.scan_indices[inds][0]
        p_counts = self.point_counts[inds][0]

        intensity = self.mz_intensity[scan_inds : scan_inds + p_counts]
        m_z = np.round(self.mz_values[scan_inds : scan_inds + p_counts], 2)

        return m_z, intensity

    def get_ion_chromatogram(
        self, clusters: list[list[float]]
    ) -> dict[float, np.ndarray]:
        """
        Get all ion chromatograms from the Chromatogram using
        pre-defined clusters of m/z values to bin signals.

        Parameters
        ----------
        clusters: list
            Clusters of m/z values. Mass values which are
            together in the list values will be combined in
            the output.

        Returns
        -------
        ion_dict: dict
            Dict of ion chromatograms
        """

        ion_dict = {sum(c) / len(c): np.zeros(len(self.time)) for c in clusters}
        if len(self.scan_indices) != 0:
            scan_brackets = []

            for s in range(0, len(self.scan_indices) - 1):
                scan_brackets.append([self.scan_indices[s], self.scan_indices[s + 1]])

            for s, bracket in enumerate(scan_brackets):
                st_bracket = bracket[0]
                end_bracket = bracket[1]
                inten = self.mz_intensity[st_bracket:end_bracket]
                masses = self.mz_values[st_bracket:end_bracket]

                for m in range(0, len(masses)):
                    for _, c in enumerate(clusters):
                        if masses[m] in c:
                            ion_dict[sum(c) / len(c)][s] = inten[m]
                            break

        return ion_dict

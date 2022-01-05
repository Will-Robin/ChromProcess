import numpy as np

from ChromProcess import Classes

class Chromatogram:
    '''
    A class for storing chromatographic data.
    '''
    def __init__(self):
        '''
        Initialise an empty chromatogram.

        Attributes
        ----------
        self.filename: str
            Name of the file containing the source data.
        self.x_unit: str
            Name of the x units (prefer 'unit name/ SI symbol')
        self.y_unit: str
            Name of the y units (prefer 'unit name/ SI symbol')
        self.time: numpy array
            Time axis for the chromatogram.
        self.signal: numpy array
            Signal of the chromatogram.
        self.peaks: dict
            Container for the peaks derived from the chromatogram.
        self.mz_values: list/numpy array
            Container for the m/z values from mass spectra.
        self.mz_intensity: list/numpy array
            Container for the intensity values from mass spectra.
        self.scan_indices: list/numpy array
            Container for the scan indices from mass spectra.
        self.point_counts: list/numpy array
            Container for the point counts from mass spectra.
        self.internal_standard: ChromProcess.Classes.Peak
            The internal standard peak.
        '''

        self.filename = ''
        self.x_unit = ''
        self.y_unit = ''

        self.time = []
        self.signal = []

        self.peaks = {}

        self.mz_values = []
        self.mz_intensity = []
        self.scan_indices = []
        self.point_counts = []

        self.internal_standard = Classes.Peak(0.0, [])

    def get_mass_spectrum(self, time):
        '''
        Get the mass spectrim at a given time point in the chromatogram.

        time: float
        mass: numpy array
        intensity: numpy array
        '''

        inds = np.where(self.time == time)[0]

        scan_inds = self.scan_indices[inds][0]
        p_counts = self.point_counts[inds][0]

        intensity = self.mz_intensity[scan_inds:scan_inds+p_counts]
        m_z = np.round(self.mz_values[scan_inds:scan_inds+p_counts], 2)

        return m_z, intensity

    def ion_chromatogram(self, clusters):
        '''
        Get all ion chromatograms from the Chromatogram using
        pre-defined clusters of m/z values to bin signals.

        clusters: dict
            Clusters of m/z values. Mass values which are
            together in the list values will be combined in
            the output.
        ion_chromatograms: dict
            Dict of ion chromatograms
        '''

        ion_dict = {}
        if len(self.scan_indices) != 0:
            ion_dict = {
                        np.average(c): np.zeros(len(self.time))
                                                 for c in clusters
                        }

            scan_brackets = []

            for s in range(0,len(self.scan_indices)-1):
                scan_brackets.append(
                                    [
                                    self.scan_indices[s],
                                    self.scan_indices[s+1]
                                    ]
                                )

            for s,bracket in enumerate(scan_brackets):
                st_bracket  = bracket[0]
                end_bracket = bracket[1]
                inten = self.mz_intensity[st_bracket:end_bracket]
                masses = self.mz_values[st_bracket:end_bracket]

                for m in range(0,len(masses)):
                    for _,c in enumerate(clusters):
                        if masses[m] in c:
                            ion_dict[np.average(c)][s] = inten[m]
                            break

        return ion_dict


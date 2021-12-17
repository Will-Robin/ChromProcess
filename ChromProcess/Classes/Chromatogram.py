import sys
import numpy as np

from ChromProcess import Classes

class Chromatogram:
    '''
    A class for storing chromatographic data.
    '''
    def __init__(self):
        '''
        Initialise an empty chromatogram.

        Parameters
        ----------
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

        self.internal_reference = Classes.Peak(0.0, [])

    def write_to_csv(self, filename = 'chromatogram.csv'):
        '''
        Write chromatgram to a .csv file.
        '''

        time   = self.time
        signal = self.signal

        with open(filename, 'w') as f:
            f.write(f'{self.x_unit},{self.y_unit}\n')
            for x in range(0,len(time)):
                f.write(f'{time[x]},{signal[x]}\n')

    def write_peak_collection_text(self, header_text = ''):
        '''
        Create the text for a peak collection based on the Peak objects in the
        chromatogram.

        Parameters
        ----------
        header_text: string 
        '''

        peak_collection_string = ''

        if header_text != '':
            peak_collection_string += header_text

        peak_collection_string += 'IS_retention_time/ min,'
        peak_collection_string += 'IS_integral,IS_peak start/ min,'
        peak_collection_string += 'IS_peak end/ min\n'

        IS_entry = self.internal_reference.peak_table_entry_text(self)

        peak_collection_string += IS_entry

        peak_collection_string += "Retention_time/ min,"
        peak_collection_string += 'integral,peak start/ min,'
        peak_collection_string += "peak end/ min\n"

        for p in self.peaks:
            peak = self.peaks[p]
            peak_collection_string += peak.peak_table_entry_text(self)

        return peak_collection_string

    def write_peak_collection(
                            self, 
                            filename = 'peak_collection.csv',
                            header_text = ""
                            ):
        '''
        For writing peak integrals from a chromatogram to a .csv file.

        Parameters
        ----------
        filename: str or pathlib Path
            Name for the file
        header_text: str
           Text to place at the top of the peak table. 
        Returns
        -------
        None
        '''

        output_text = self.write_peak_collection_text()

        with open(filename, "w") as f:
            f.write(output_text, header_text = header_text)

    def get_mass_spectrum(self, time):
        '''
        Get the mass spectrim at a given time point in the chromatogram.

        time: float
        mass: numpy array
        intensity: numpy array
        '''
        import numpy as np

        inds = np.where(self.time == time)[0]

        scan_inds = self.scan_indices[inds][0]
        p_counts = self.point_counts[inds][0]

        intensity = self.mass_intensity[scan_inds:scan_inds+p_counts]
        mass = np.round(self.mass_values[scan_inds:scan_inds+p_counts], 2)

        return mass, intensity

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
                inten = self.mass_intensity[st_bracket:end_bracket]
                masses = self.mass_values[st_bracket:end_bracket]

                for m in range(0,len(masses)):
                    for _,c in enumerate(clusters):
                        if masses[m] in c:
                            ion_dict[np.average(c)][s] = inten[m]
                            break

        return ion_dict


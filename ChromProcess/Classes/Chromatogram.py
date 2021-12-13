import sys

class Chromatogram:
    def __init__(self, file,
                mass_spec = False,
                channel_select = '360nm'):
        '''
        Initialises a chromatogram object from a file

        Parameters
        ----------
        file: str or pathlib Path
            Path to file.
        mass_spec: bool
            If GC-MS data, whether to load mass spectra information (True) or
            not (False).
        channel_select: str
            For HPLC data, which detector channel to select.
        '''
        import numpy as np
        from pathlib import Path

        if isinstance(file, Path):
            self.initialised_path = file
            self.filename = file.stem.split('.')[0]
            self.filetype = file.suffix.strip('.')
        elif isinstance(file, str):
            # convert the path as a string into a Path object.
            self.initialised_path = Path(file)
            self.filename = self.initialised_path.stem.split('.')[0]
            self.filetype = self.initialised_path.suffix.strip('.')
        else: sys.exit(
                    '''ChromProcess.Classes.Data.Chromatogram:
                        expected file arg to be string or pathlib Path.'''
                        )

        self.peaks = {}

        self.mass_spectra = False

        self.mass_values = []
        self.mass_intensity = []
        self.scan_indices = []
        self.point_counts = []
        self.internal_reference = False

        # act according to the suffix of the file type passed in.
        if self.filetype == 'txt':
            data  = self.get_data_Shimadzu_HPLC(
                                    self.initialised_path,
                                    self.get_info_Shimadzu_HPLC(file)
                                    )

            try:
                chan_ind = data["Wavelength"].index(channel_select)
                self.time = np.array(data["Time"][chan_ind])
                self.signal = np.array(data["Signal"][chan_ind])
            except:
                wvl = data["Wavelength"]
                print(f"{channel_select} not found, try: {wvl}")
                sys.exit()

            self.c_type = 'HPLC'
            self.mass_spectra = False

        elif self.filetype == 'cdf':
            self.time   = self.get_data_cdf_GCMS(
                                    self.initialised_path,
                                    'scan_acquisition_time'
                                    )/60 # converted to minutes
            self.signal = self.get_data_cdf_GCMS(
                                                self.initialised_path,
                                                'total_intensity'
                                                )
            self.c_type = 'GCMS'

            if mass_spec == True:
                self.MS_Load()

        elif self.filetype == 'csv':
            self.time, self.signal = self.load_from_csv(
                                                self.initialised_path
                                                    )
            self.c_type = 'from_csv'

        elif self.filetype == 'CSV':
            self.time, self.signal = self.load_from_csv(
                                                self.initialised_path
                                                    )
            self.c_type = 'from_csv'
        else:
            f_type = self.filetype
            print(f'Unexpected file type ({f_type}). File not loaded')

    def get_data_Shimadzu_HPLC(self, file, dat_dict):
        '''
        A specific parser for .txt files exported from Shimadzu LabSolutions
        software.

        Extracts data from the chromatogram file into a dictionary.

        Parameters
        ----------
        file: str or pathlib Path
            name of chromatogram file
            (ASCII exported from Shimadzu LC software)
        dat_dict: dict
            dictionary of information from chromatogram file.

        Returns
        -------
        Updated dat_dict with chromatogram data added from the file

        '''
        import numpy as np

        for b in range(0,len(dat_dict['Data set start'])):
            with open(file, 'r') as f:
                time = []
                signal = []
                for c,line in enumerate(f):
                    if c > dat_dict['Data set start'][b]:
                        ex = line.strip("\n")
                        inp = ex.split("\t")
                        inp = [e.replace(',','.') for e in inp]
                        time.append(float(inp[0]))
                        signal.append(float(inp[1]))
                        end_time = dat_dict['End Time'][b]
                        if round(float(inp[0]),3) == end_time:
                            dat_dict['Time'].append(np.array(time))
                            dat_dict['Signal'].append(np.array(signal))
                            break
        return dat_dict

    def get_info_Shimadzu_HPLC(self, file):
        '''
        Extracts key information from .txt files exported from Shimadzu LabSolutions
        software.

        Parameters
        ----------
        file: str
             name of chromatogram file (ASCII exported from
             Shimadzu LC software)
        Returns
        -------
        dat_dict: dict
            dictionary of information extracted from file

        '''

        dat_dict = {"Data set start"   : [],
                    "Start Time"       : [],
                    "End Time"         : [],
                    "Intensity Units"  : [],
                    "Time"             : [],
                    "Signal"           : [],
                    "Wavelength"       : []}

        with open(file, 'r') as f:
            for c,line in enumerate(f):
                if 'Start Time'in line:
                    ex =  line.split('\t')
                    ex = [e.replace(',','.') for e in ex]
                    start_time = float(ex[1].strip('\n'))
                    dat_dict['Start Time'].append(start_time)
                if 'R.Time (min)' in line:
                    dat_dict["Data set start"].append(c+1)
                if 'End Time'in line:
                    ex =  line.split('\t')
                    ex = [e.replace(',','.') for e in ex]
                    end_time = float(ex[1].strip('\n'))
                    dat_dict['End Time'].append(end_time)
                if 'Intensity Units'in line:
                    ex =  line.split('\t')
                    ex = [e.replace(',','.') for e in ex]
                    intensity = ex[1].strip('\n')
                    dat_dict['Intensity Units'].append(intensity)
                if 'Wavelength(nm)' in line:
                    ex =  line.split('\t')
                    ex = [e.replace(',','.') for e in ex]
                    wavelength = ex[1].strip('\n')
                    dat_dict['Wavelength'].append(wavelength)

        return dat_dict

    def get_data_cdf_GCMS(self, file, key):
        '''
        Extracts data from a .cdf file using the Dataset function
        from the netCDF4 library.

        Parameters
        ----------
        file: str
            file name of the GCMS .cdf file
        key: str
            key to a set of data in the .cdf file

        Returns
        -------
        output: numpy array
            array filled with values of selected variable

        Notes
        -----
        For a list of contents of the .cdf file, use cdf_info function.
        '''
        from netCDF4 import Dataset

        f = Dataset(file, "r")
        f.set_auto_mask(False)
        output = f.variables[key][:]
        f.close()
        return output

    def MS_Load(self):

        import numpy as np

        file = self.initialised_path

        self.mass_spectra = True
        # Measured intensities
        self.mass_intensity = self.get_data_cdf_GCMS(file, "intensity_values")
        # Measured masses
        self.mass_values = np.round(
                                    self.get_data_cdf_GCMS(
                                                    file,
                                                    "mass_values"
                                                    ),
                                    3)
        # Scan index is starting index
        self.scan_indices = self.get_data_cdf_GCMS(file, "scan_index")
        # Point count is number of elements to read
        self.point_counts = self.get_data_cdf_GCMS(file, "point_count")

    def load_from_csv(self,file, read_start_index = 1):
        '''
        For loading a chromatogram from a .csv file
        Parameters
        ----------
        file: str
            Path to file
        Returns
        -------
        time: numpy array
            time axis
        signal: numpy array
            signal axis
        '''
        import numpy as np

        data = []
        with open(file, 'r') as f:

            for c,line in enumerate(f):
                if c < read_start_index:
                    pass
                else:
                    line_as_list = line.strip('\n').split(',')
                    insertion = [float(x) for x in line_as_list]
                    data.append(insertion)

        data = np.array(data)

        time = data[:,0]
        signal = data[:,1]

        return time, signal

    def write_to_csv(self, filename = ''):
        '''
        Write chromatgram to a .csv file.
        '''
        import numpy as np

        if filename == '':
            filename = self.filename

        time   = self.time
        signal = self.signal
        chrom_out = np.vstack((time,signal))
        chrom_out = chrom_out.T

        with open(filename, 'w') as f:
            f.write('time/ min, signal/ total ion counts\n')
            for x in range(0,len(chrom_out)):
                f.write(f'{chrom_out[x,0]},{chrom_out[x,1]}\n')

    def write_peak_collection_text(self, value = 0, series_unit = "None"):
        '''
        Create the text for a peak collection based on the Peak objects in the
        chromatogram.
        series_unit: string (will be coverted to string).
        value: float or string (will be coverted to string).
        peak_collection_string: string
        '''
        peak_collection_string = ''

        peak_collection_string += f"{series_unit},{value}\n"
        peak_collection_string += '(IS_retention_time/ min,IS_integral,IS_peak start/ min,IS_peak end/ min\n'

        if self.internal_reference:
            st_ind = self.internal_reference.indices[0]
            end_ind = self.internal_reference.indices[-1]
            IS_lower_bound = self.time[st_ind]
            IS_upper_bound = self.time[end_ind]
            IS_RT = self.internal_reference.retention_time
            IS_integral = self.internal_reference.integral
        else:
            IS_RT, IS_integral, IS_lower_bound, IS_upper_bound = 'None', 'None', 'None', 'None'

        peak_collection_string += f"{IS_RT},{IS_integral},{IS_lower_bound},{IS_upper_bound}\n"

        peak_collection_string += "Retention_time/ min,integral,peak start/ min,peak end/ min\n"

        for p in self.peaks:
            st_ind = self.peaks[p].indices[0]
            end_ind = self.peaks[p].indices[-1]

            peak_start = self.time[st_ind]
            peak_end = self.time[end_ind]

            rtn_time = self.peaks[p].retention_time
            integral = self.peaks[p].integral

            peak_collection_string += f"{rtn_time},{integral},{peak_start},{peak_end}\n"

        return peak_collection_string

    def write_peak_collection(self, filename = 'Peak_collection',
                        value = 0.0,
                        series_unit = "None"):
        '''
        For writing peak integrals from a chromatogram to a .csv file.

        Parameters
        ----------
        filename: str or pathlib Path
            Name for the file
        value: float, str, bool
            Value assigned to the variable assigned to the chromatogram.
        series_unit: float, str, bool
            Units of the variable assigned to the chromatogram.
        Returns
        -------
        None
        '''

        with open(filename, "w") as f:
            f.write("{},{}\n".format(series_unit,value))
            f.write('IS_retention_time/ min,IS_integral,IS_peak start/ min,IS_peak end/ min\n')
            if self.internal_reference:
                st_ind = self.internal_reference.indices[0]
                end_ind = self.internal_reference.indices[-1]
                IS_lower_bound = self.time[st_ind]
                IS_upper_bound = self.time[end_ind]
                IS_RT = self.internal_reference.retention_time
                IS_integral = self.internal_reference.integral
            else:
                IS_RT, IS_integral, IS_lower_bound, IS_upper_bound = 'None', 'None', 'None', 'None'

            f.write("{},{},{},{}\n".format(
                                            IS_RT,
                                            IS_integral,
                                            IS_lower_bound,
                                            IS_upper_bound
                                        )
                    )
            f.write("Retention_time/ min,integral,peak start/ min,peak end/ min\n")

            for p in self.peaks:
                st_ind = self.peaks[p].indices[0]
                end_ind = self.peaks[p].indices[-1]
                f.write("{},{},{},{}\n".format(
                                                self.peaks[p].retention_time,
                                                self.peaks[p].integral,
                                                self.time[st_ind],
                                                self.time[end_ind]
                                            )
                        )
    def write_peak_mass_spectra_text(self):
        '''
        Writes the text for a table of mass spectra derived from peaks.
        ms_text: string
        '''

        ms_text = ''

        for p in self.peaks:
            rtn_time = self.peaks[p].retention_time
            ms_text += f'Peak retention time, {rtn_time}\n'
            if self.peaks[p].mass_spectrum:
                ms = self.peaks[p].mass_spectrum

                ms_text += 'm/z,'
                ms_text += ','.join(ms[0])
                ms_text += '\n'

                ms_text += 'relative abundance,'
                rel_inten = [x/max(ms[1]) for x in ms[1]]
                ms_text += ','.join(rel_inten)
                ms_text += '\n'
            else:
                ms_text += 'm/z,empty\n'
                ms_text += 'relative abundance,empty\n'

        return ms_text

    def write_peak_mass_spectra(self, filename = ''):
        '''
        For writing mass spectra of chromatogram peaks to a .csv file.

        Parameters
        ----------

        filename: str
            Name for the file

        Returns
        -------
        None
        '''

        with open(filename, 'w') as f:

            for p in self.peaks:
                if self.peaks[p].mass_spectrum:
                    ms = self.peaks[p].mass_spectrum

                    f.write('Peak retention time, {}\n'.format(
                                            self.peaks[p].retention_time
                                                    )
                            )
                    f.write('m/z,')
                    [f.write('{},'.format(x)) for x in ms[0]]
                    f.write('\n')
                    f.write('relative abundance,')
                    [f.write('{},'.format(x/max(ms[1]))) for x in ms[1]]
                    f.write('\n')
                else:
                    f.write('Peak retention time, {}\n'.format(
                                            self.peaks[p].retention_time
                                                    )
                            )
                    f.write('m/z,empty\n')
                    f.write('relative abundance,empty\n')

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
        import numpy as np

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


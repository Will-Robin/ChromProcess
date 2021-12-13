import sys

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
        '''

        self.retention_time = retention_time

        self.indices = indices

        self.integral = False

        self.height = False

        self.ion_chromatograms = {}

        self.ion_integrals = {}

        self.mass_spectrum = False

        self.character = "Unknown"

        self.deconvolution = None

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
        '''

        import numpy as np

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
        '''

        self.height = chromatogram.signal[self.indices]

    def get_mass_spectrum(self, chromatogram):
        '''
        Get the mass spectrum at the apex of the peak.

        Parameters
        ----------
        chromatogram: ChromProcess Chromatogram object
        '''
        import numpy as np

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

class PeakCollection:
    def __init__(self, file = ''):

        from ChromProcess import Classes


        self.filename = 'not specified'
        self.series_value = None
        self.series_unit = 'not specified'
        self.internal_standard = Classes.PeakCollectionElement(0,0,0,0)
        self.peaks = [Classes.PeakCollectionElement(0,0,0,0)]
        self.assignment = 'not specified'
        self.mass_spectra = []
        self.initial_IS_pos = 0.0
        self.assigned_compounds = []

        if file == '':
            pass
        else:
            self.read_from_file(file)

    def read_from_file(self, file):
        '''
        Read in information to the object from a file.
        '''

        from pathlib import Path
        from ChromProcess import Classes

        fname = Path('')
        if isinstance(file, str):
            fname = Path(file)
        elif isinstance(file, Path):
            fname = file
        else:
            sys.exit('''PeakCollection requires file kwarg to be string or
            pathlib Path.''')

        self.filename = fname.name

        read_line = lambda line: [float(x) for x in line.strip('\n').split(",") if x != '']
        peaks = []

        IS_line_num = -1
        IS = Classes.PeakCollectionElement(0.0, 1, 0.0, 0.0)

        value = 0.0
        variable = ''
        with open(file, "r") as f:
            for c,line in enumerate(f):
                if 'None' in line:
                    pass
                elif c == 0:
                    read = [x for x in line.strip('\n').split(',') if x != '']
                    variable = read[0]
                    value = float(read[1])
                elif 'IS_' in line:
                    IS_line_num = c + 1

                elif c == IS_line_num:
                    if 'None' in line:
                        pass
                    else:
                        read = read_line(line)

                        IS = Classes.PeakCollectionElement(round(read[0],3),
                                       read[1],
                                       round(read[2],3),
                                       round(read[3],3),
                                       parent = self.filename.split('.')[0])
                elif c < 4:
                    pass
                else:
                    rd = read_line(line)
                    peaks.append(
                    Classes.PeakCollectionElement(round(rd[0],3), rd[1],
                                                  round(rd[2],3),
                                                  round(rd[3],3),
                                                  parent = self.filename.split('.')[0]))

        self.series_value = value
        self.series_unit = variable
        self.internal_standard = IS
        self.peaks = peaks
        self.assignment = 'not specified'
        self.mass_spectra = []
        self.initial_IS_pos = IS.retention_time

    def inspect_peaks(self):
        for p in self.peaks:
            p.inspect_peak()

    def remove_peaks_below_threshold(self,threshold):
        '''
        Parameters
        ----------
        threshold: float
        '''
        del_idx = []
        for c,pk in enumerate(self.peaks):
            if pk.integral < threshold:
                del_idx.append(c)

        self.peaks = [v for i,v in enumerate(self.peaks) if i not in del_idx]

    def align_peaks_to_IS(self, IS_set = 0.0):
        '''
        Parameters
        ----------
            IS_set: float
        '''

        is_rt = self.internal_standard.retention_time

        for p in self.peaks:
            p.retention_time = p.retention_time - is_rt + IS_set
            p.start = p.start - is_rt + IS_set
            p.end = p.end - is_rt + IS_set

        for m in self.mass_spectra:
            m.retention_time = m.retention_time - is_rt + IS_set

        self.internal_standard.start = (self.internal_standard.start - is_rt +
                                                                        IS_set)
        self.internal_standard.end = (self.internal_standard.end - is_rt +
                                                                        IS_set)
        self.internal_standard.retention_time = (
                                           self.internal_standard.retention_time
                                           - is_rt + IS_set)
    def add_mass_spectra(self, ms_list):
        '''
        Add mass spectrum information into PeakCollectionElement objects.

        ms_list: list of ChromProcess MassSpectrum objects
            Mass spectra to be added.
        '''
        ms_dict = {}
        for m in ms_list:
            ms_dict[m.retention_time] = m

        peak_dict = {}
        for pk in self.peaks:
            peak_dict[pk.retention_time] = pk

        for p in peak_dict:
            if p in ms_dict:
                peak_dict[p].mass_spectrum = ms_dict[p]

    def get_peak_positions(self):
        import numpy as np
        return np.array([p.retention_time for p in self.peaks])

    def reference_integrals_to_IS(self):
        IS_integral = self.internal_standard.integral
        for p in self.peaks:
            p.reference_integral_to_IS(IS_integral)

    def assign_peaks(self, boundaries):
        '''
        Parameters
        ----------
        boundaries: dict
            {'compound_name': [lower bound, upper bound]}
        '''
        for pk in self.peaks:
            pk.assign_peak(boundaries)

    def apply_calibrations_to_peaks(self, calibrations, IS_conc):
        '''
        Applies calibrations

        Parameters
        ----------
        calibrations: ChromProcess Instrument_Calibration object
            Container for calibration information
        IS_conc: float
            internal standard concentration

        '''

        if IS_conc == 0.0:
            # results in division by 1 during conversion
            IS_conc = 1.0
        for pk in self.peaks:
            if pk.assignment in calibrations.calibration_factors:
                CF = calibrations.calibration_factors[pk.assignment]
                if calibrations.calibration_model == 'linear':
                    pk.apply_linear_calibration(CF['B'], CF['C'],
                                                internal_standard = IS_conc)
                elif calibrations.calibration_model == 'quadratic':
                    pk.apply_quadratic_calibration(CF['A'], CF['B'], CF['C'],
                                                   internal_standard = IS_conc)

    def calculate_conc_errors(self, calibrations, IS_conc, IS_conc_err):
        '''
        Calculation of the standard error on a concentration estimation from
        th calibration.

        Parameters
        ----------
        calibrations: ChromProcess Instrument_Calibration object
            Contains calibration information.
        IS_conc: float
        IS_conc_err: float

        Returns
        -------
        None

        Modifies PeakCollectionElement object attributes.
        '''

        for pk in self.peaks:
            pk.calculate_error(calibrations,IS_conc,IS_conc_err)

    def dilution_correct_peaks(self, analysis):
        '''
        Parameters
        ----------
        analysis: ChromProcess Analysis_Information object
        '''

        dil = analysis.dilution_factor
        dil_err = analysis.dilution_factor_error
        for pk in self.peaks:
            if pk.concentration:
                pk.dilution_correction(dil, dil_err)

    def get_all_assigned_compounds(self):
        assigns = []
        for pk in self.peaks:
            if pk.assignment == 'none':
                pass
            else:
                assigns.append(pk.assignment)

        assigns = list(set(assigns))
        self.assigned_compounds = sorted(assigns, key = lambda x:x.count('C'))


    def write_to_file(self, directory = ''):
        '''
        Write the object to a structured file.

        Parameters
        ----------
        directory: str or pathlib Path
        '''

        from pathlib import Path

        fname = Path(self.filename)
        if isinstance(directory, str):
            if directory == '':
                fname = Path(f"{self.filename}")
        elif isinstance(directory, Path):
            fname = directory/f"{self.filename}"

        IS_header = 'IS_retention_time/ min,IS_integral,IS_peak start/ min,IS_peak end/ min\n'
        pk_header = "Retention_time/ min,integral,peak start/ min,peak end/ min\n"

        with open(fname, "w") as f:
            f.write("{},{}\n".format(self.series_unit,self.series_value))
            f.write(IS_header)
            IS_rt = self.internal_standard.retention_time
            IS_integ = self.internal_standard.integral
            strt = self.internal_standard.start
            end = self.internal_standard.end
            f.write("{},{},{},{}\n".format(IS_rt, IS_integ,strt,end))
            f.write(pk_header)
            for p in self.peaks:
                rt = p.retention_time
                integ = p.integral
                start = p.start
                end = p.end
                f.write("{},{},{},{}\n".format(rt, integ, start, end))

class PeakCollectionSeries:
    def __init__(self, peak_collections, name = 'not specified',
                 conditions = {}):
        '''
        An object which wraps multiple PeakCollection objects.

        Parameters
        ----------
        peak_collections: List of ChromProcess PeakCollection objects
        name: str
        conditions: dict
        '''

        self.name = name
        self.peak_collections = peak_collections
        self.series_values = [pt.series_value for pt in peak_collections]
        self.series_unit = peak_collections[0].series_unit
        self.conditions = conditions
        self.clusters = []
        self.cluster_assignments = []
        self.integral_series = []
        self.concentration_series = []
        self.conc_err_series = []
        self.series_assigned_compounds = []

    def remove_peaks_below_threshold(self,threshold):
        '''
        Parameters
        ----------
        threshold: float (from 0.0 to 1.0)
        '''
        for pc in self.peak_collections:
            pc.remove_peaks_below_threshold(threshold)

    def align_peaks_to_IS(self, IS_set):
        '''
        Parameters
        ----------
        IS_set: float
        '''
        for pc in self.peak_collections:
            pc.align_peaks_to_IS(IS_set = IS_set)

    def get_peak_positions(self):

        import numpy as np

        peak_pos = np.array([])
        for pc in self.peak_collections:
            peak_pos = np.hstack((peak_pos, pc.get_peak_positions()))

        i = np.argsort(peak_pos)

        return peak_pos[i]

    def get_peak_clusters(self, bound = 0.1):
        '''
        bound: float
        '''

        from ChromProcess import simple_functions as sf

        peaks = self.get_peak_positions()

        clusts = []
        for c in sf.cluster(peaks, bound = bound):
            clusts.append(c)

        self.clusters = clusts

    def assign_peaks(self, boundaries):
        '''
        Parameters
        ----------
        boundaries: dict
            {'compound_name', [lower bound, upper bound]}
        '''

        for pc in self.peak_collections:
            pc.assign_peaks(boundaries)

    def assign_clusters(self,boundaries):
        '''
        Parameters
        ----------
        boundaries: dict
            {'compound_name', [lower bound, upper bound]}
        '''
        import numpy as np
        from ChromProcess import processing_functions as p_f

        for c in self.clusters:
            pos = np.mean(c)
            clust_name = p_f.name_peak(pos, boundaries)
            self.cluster_assignments.append(clust_name)

    def get_all_assigned_compounds(self):
        assigns = []
        for pc in self.peak_collections:
            if len(pc.assigned_compounds) == 0:
                pc.get_all_assigned_compounds()
            assigns.extend( pc.assigned_compounds )

        assigns = list(set(assigns))
        self.series_assigned_compounds = sorted(assigns, key = lambda x:x.count('C'))

    def reference_integrals_to_IS(self):
        for pc in self.peak_collections:
            pc.reference_integrals_to_IS()

    def apply_calibrations(self, conditions, calibrations):
        '''
        Parameters
        ----------
        conditions: ChromProcess Analysis_Information object
            container for analysis information

        calibrations: ChromProcess Instrument_Calibration object
            Contains calibration information.
        '''

        IS_conc = conditions.internal_ref_concentration

        for pc in self.peak_collections:
            pc.apply_calibrations_to_peaks(calibrations, IS_conc)

    def calculate_conc_errors(self, calib, conditions):
        '''
        Parameters
        ----------
        calib: ChromProcess Instrument_Calibration object
        '''
        IS_conc = conditions.internal_ref_concentration
        IS_conc_err = conditions.internal_ref_concentration_error

        for pc in self.peak_collections:
            pc.calculate_conc_errors(calib,IS_conc,IS_conc_err)

    def apply_peak_dilution_factors(self,analysis):
        '''
        Parameters
        ----------
        analysis: ChromProcess Analysis_Information object
        '''
        for pc in self.peak_collections:
            pc.dilution_correct_peaks(analysis)

    def make_integral_series(self, cluster_bound = 0.0):

        '''
        Parameters
        ----------
        cluster_bound: float
        '''

        import numpy as np

        if len(self.clusters) == 0:
            self.get_peak_clusters(bound = cluster_bound)

        series_courses = np.zeros((len(self.series_values), len(self.clusters)))

        for c1,pc in enumerate(self.peak_collections):
            for c2,clust in enumerate(self.clusters):
                for pk in pc.peaks:
                    if pk.retention_time == pc.internal_standard.retention_time:
                        continue
                    if pk.retention_time in clust and pk.integral:
                        series_courses[c1,c2] += pk.integral

        self.integral_series = series_courses.T

    def make_concentration_series(self, cluster_bound = 0.0):
        '''
        Parameters
        ----------
        cluster_bound: float
        '''

        import numpy as np

        if len(self.clusters) == 0:
            self.get_peak_clusters(bound = cluster_bound)

        series_courses = np.zeros((len(self.series_values), len(self.clusters)))
        error_courses = np.zeros((len(self.series_values), len(self.clusters)))

        clust_assigns = []
        for c1,pc in enumerate(self.peak_collections):
            for c2,clust in enumerate(self.clusters):
                loc_assigns = []
                for pk in pc.peaks:
                    if pk.retention_time == pc.internal_standard.retention_time:
                        continue
                    if pk.retention_time in clust and pk.concentration:
                        loc_assigns.append(pk.assignment)
                        series_courses[c1,c2] += pk.concentration
                        error_courses[c1,c2] += pk.conc_error

                clust_assigns.append(loc_assigns)

        self.concentration_series = series_courses.T
        self.conc_err_series = error_courses.T

    def concentration_traces_as_dict(self):
        '''
        Parameters
        ----------
        name_conversions: dict
        '''

        import numpy as np
        from ChromProcess import simple_functions as s_f

        if len(self.cluster_assignments) == 0:
            get_name = lambda _: ''
        else:
            get_name = lambda x: self.cluster_assignments[x].split(' ')[0]

        conc_dict = {}
        for x in range(0,len(self.concentration_series)):
            name = get_name(x)
            pos = np.round(np.mean(self.clusters[x]),3)
            if name == '':
                pass
            elif not s_f.isfloat(name):
                conc_dict[f'{name}/ M ({pos})'] = self.concentration_series[x]
            else:
                pass

        return conc_dict

    def integral_traces_as_dict(self):

        import numpy as np

        if len(self.cluster_assignments) == 0:
            get_name = lambda _: ''
        else:
            get_name = lambda x: self.cluster_assignments[x].split(' ')[0]

        integral_dict = {}
        for x in range(0,len(self.integral_series)):
            name = get_name(x)
            if name != '':
                pos = np.mean(self.clusters[x])
                token = f'{name} ({np.round(pos,3)})'
            else:
                cluster_average = np.mean(self.clusters[x])
                val = np.round(cluster_average, 3)
                token = f'{val} ({val})'

            integral_dict[token] = self.integral_series[x]

        return integral_dict

    def concentration_error_traces_dict(self):

        import numpy as np
        from ChromProcess import simple_functions as s_f

        err_dict = {}

        if len(self.cluster_assignments) == 0:
            get_name = lambda _: ''
        else:
            get_name = lambda x: self.cluster_assignments[x].split(' ')[0]

        for x in range(0,len(self.conc_err_series)):
            name = get_name(x)
            pos = np.round(np.mean(self.clusters[x]),3)
            if name == '':
                pass
            elif not s_f.isfloat(name):
                err_dict[f'{name}/ M ({pos})'] = self.conc_err_series[x]
            else:
                pass

        return err_dict

    def write_conditions_header(self, outfile, information):
        '''
        Parameters
        ----------
        outfile: Python file object

        information: ChromProcess Instrument_Calibration object
        '''

        # writing experiment conditions to file
        outfile.write("Dataset,{}\n".format(self.name))
        outfile.write("start_conditions\n")
        for c in self.conditions:
            outfile.write("{},".format(c))
            [outfile.write("{},".format(x)) for x in self.conditions[c]]
            outfile.write("\n")
        outfile.write("end_conditions\n")
        # writing analysis details
        outfile.write("start_analysis_details\n")
        outfile.write('Instrument, {}\n'.format(information.instrument))
        outfile.write("Chromatography_method,{},{}\n".format(information.type, information.method))
        outfile.write("Derivatisation_method,{}\n".format(information.derivatisation))
        outfile.write("Calibrations_file,{}\n".format(information.filename.name))
        outfile.write('Calibration_model,{}\n'.format(information.calibration_model))
        outfile.write("end_analysis_details\n")

    def write_concentrations_to_file(self,
                                    filename, information,
                                    ):
        '''
        Parameters
        ----------
        filename: name for file including path

        information: ChromProcess Instrument_Calibration object
        '''

        import numpy as np
        from pathlib import Path

        if isinstance(filename, str):
            filename = filename
        elif isinstance(filename, Path):
            filename = str(filename)

        out_type = 'concentration_report'
        fname = '{}_{}_{}.csv'.format(filename, information.type, out_type)

        with open(fname, 'w') as outfile:
            # writing experiment conditions to file
            self.write_conditions_header(outfile,information)

            # writing data
            conc_traces = self.concentration_traces_as_dict()
            sorted_keys = sorted([*conc_traces], key = lambda x:x.count('C'))

            outfile.write("start_data\n")

            p_header = [self.series_unit]
            out = np.array([self.series_values])

            for s in sorted_keys:
                p_header.append(s)
                out = np.vstack((out,conc_traces[s]))

            out = out.T
            [outfile.write("{},".format(x)) for x in p_header]

            outfile.write("\n")

            for x in range(0,len(out)):
                for y in range(0,len(out[x])):
                    outfile.write("{},".format(out[x,y]))
                outfile.write("\n")

            outfile.write("end_data\n")

    def write_integrals_to_file(self,
                                filename,
                                information,
                                ):
        '''
        Parameters
        ----------
        filename: name for file including path

        information: ChromProcess Instrument_Calibration object
        '''

        import numpy as np
        from pathlib import Path

        name = ''
        if isinstance(filename, str):
            name = filename
        elif isinstance(filename, Path):
            name = str(filename)

        out_type = 'integral_report'
        fname = '{}_{}_{}.csv'.format(name, information.type, out_type)

        with open(fname, 'w') as outfile:
            # writing experiment conditions to file
            self.write_conditions_header(outfile,information)

            # writing data
            integral_traces = self.integral_traces_as_dict()

            outfile.write("start_data\n")

            p_header = [self.series_unit]
            out = np.array([self.series_values])

            for s in [*integral_traces]:
                p_header.append(s)
                out = np.vstack((out,integral_traces[s]))

            out = out.T
            [outfile.write("{},".format(x)) for x in p_header]

            outfile.write("\n")

            for x in range(0,len(out)):
                for y in range(0,len(out[x])):
                    outfile.write("{},".format(out[x,y]))
                outfile.write("\n")

            outfile.write("end_data\n")

    def create_DataReport_base(self, information):
        '''
        Parameters
        ----------
        information: ChromProcess Instrument_Calibration object
        '''

        from ChromProcess import Classes

        data_report = Classes.DataReport()
        data_report.experiment_code = self.name
        data_report.conditions = self.conditions
        data_report.analysis_details['Instrument'] = information.instrument
        chrom_method = [information.type, information.method]
        data_report.analysis_details['Chromatography_method'] = chrom_method
        deriv_info = information.derivatisation
        data_report.analysis_details['Derivatisation_method'] = deriv_info
        data_report.analysis_details['Calibrations_file'] = information.filename.name
        data_report.analysis_details['Calibration_model'] =information.calibration_model

        data_report.series_values = self.series_values
        data_report.series_unit = self.series_unit

        return data_report

    def create_conc_DataReport(self, information):
        '''
        Parameters
        ----------
        information: ChromProcess Instrument_Calibration object
        '''

        out_type = 'concentration_report'
        fname = '{}_{}_{}.csv'.format(self.name, information.type, out_type)

        data_report = self.create_DataReport_base(information)
        data_report.data = self.concentration_traces_as_dict()
        data_report.errors = self.concentration_error_traces_dict()
        data_report.filename = fname

        return data_report

    def create_integral_DataReport(self, information):
        '''
        Parameters
        ----------
        information: ChromProcess Instrument_Calibration object
        '''

        out_type = 'integral_report'
        fname = '{}_{}_{}.csv'.format(self.name, information.type, out_type)

        data_report = self.create_DataReport_base(information)
        data_report.data = self.integral_traces_as_dict()
        data_report.errors = {}
        data_report.filename = fname

        return data_report

class MassSpectrum:
    def __init__(self, fname, mz, inten, pos = None):
        '''
        Parameters
        ----------
        fname: str or pathlib Path
        mz: array
        inten: array
        pos: None or float
        '''

        from pathlib import Path

        filename = fname
        if isinstance(fname, str):
            filename = Path(fname)
        elif isinstance(fname, Path):
            pass
        else:
            sys.exit('''''')

        self.filename = filename
        self.mz = mz
        self.relative_abundances = inten
        self.retention_time = pos

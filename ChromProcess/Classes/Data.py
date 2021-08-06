from ChromProcess import file_import
import numpy as np
import os

class Peak:
    def __init__(self,retention_time, indices):

        self.retention_time = retention_time

        self.indices = indices

        self.integral = False

        self.ion_chromatograms = {}

        self.ion_integrals = {}

        self.mass_spectrum = False

        self.character = "Unknown"

        self.deconvolution = None

    def get_integral(self, chromatogram, baseline_subtract = False):
        import numpy as np

        time = chromatogram.time[self.indices]
        signal = chromatogram.signal[self.indices]

        if baseline_subtract:

            linterp = np.interp(time,[time[0],time[-1]],[signal[0],signal[-1]])
            self.integral = ( np.trapz(signal, x = time) ) - linterp
        else:
            self.integral = ( np.trapz(signal, x = time) )

        return self.integral

class Chromatogram:
    def __init__(self, file, mass_spec = False, channel_select = '360nm'):
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
        if type(file) != str:
            self.initialised_path = file
            self.filename = file.stem.strip(file.suffix)
            self.filetype = file.suffix.strip('.')
        else:
            self.initialised_path = file
            self.filename = str(file[:-4])
            self.filetype = file.split('.')[-1]

        self.peaks = {} # will become a dict of dicts: {region:{peak:{Peak_indices: [], Peak_start_indices: [], Peak_end_indices: []}}}

        self.mass_spectra = False

        self.mass_values = []
        self.mass_intensity = []
        self.scan_indices = []
        self.point_counts = []
        self.internal_reference = False

        if self.filetype == 'txt':
            data  = self.get_data_Shimadzu_HPLC(self.initialised_path,
                                     self.get_info_Shimadzu_HPLC(file))

            try:
                chan_ind = data["Wavelength"].index(channel_select)
                self.time = np.array(data["Time"][chan_ind])
                self.signal = np.array(data["Signal"][chan_ind])
            except:
                print(f"{channel_select} not found, try: ", data["Wavelength"])
                quit()

            self.c_type = 'HPLC'
            self.mass_spectra = False

        elif self.filetype == 'cdf':
            self.time   = self.get_data_cdf_GCMS(self.initialised_path,
                                                        'scan_acquisition_time')/60 # converted to minutes
            self.signal = self.get_data_cdf_GCMS(self.initialised_path,
                                                        'total_intensity')
            self.c_type = 'GCMS'

            if mass_spec == True:
                self.MS_Load()

        elif self.filetype == 'csv':
            self.time, self.signal = self.load_from_csv(self.initialised_path)
            self.c_type = 'from_csv'

        else:
            print('Unexpected file type ({}). File not loaded'.format(self.filetype))

    def get_data_Shimadzu_HPLC(self, file,dat_dict):
        '''
        Extracts data from the chromatogram into a dictionary.

        Parameters
        ----------
        file: str
            name of chromatogram file (ASCII exported from Shimadzu LC software)
        dat_dict: dict
            dictionary of information from chromatogram file.

        Returns
        -------
        Updated dat_dict with chromatogram data added from the file

        '''

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
                        if round(float(inp[0]),3) == dat_dict['End Time'][b]:
                            dat_dict['Time'].append(np.array(time))
                            dat_dict['Signal'].append(np.array(signal))
                            break
        return dat_dict

    def get_info_Shimadzu_HPLC(self, file):
        '''
        Extracts key information from the exported chromatography file (see dat_dict).

        Parameters
        ----------
        file: str
             name of chromatogram file (ASCII exported from Shimadzu LC software)
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
                    dat_dict['Start Time'].append(float(ex[1].strip('\n')))
                if 'R.Time (min)' in line:
                    dat_dict["Data set start"].append(c+1)
                if 'End Time'in line:
                    ex =  line.split('\t')
                    ex = [e.replace(',','.') for e in ex]
                    dat_dict['End Time'].append(float(ex[1].strip('\n')))
                if 'Intensity Units'in line:
                    ex =  line.split('\t')
                    ex = [e.replace(',','.') for e in ex]
                    dat_dict['Intensity Units'].append(ex[1].strip('\n'))
                if 'Wavelength(nm)' in line:
                    ex =  line.split('\t')
                    ex = [e.replace(',','.') for e in ex]
                    dat_dict['Wavelength'].append(ex[1].strip('\n'))

        return dat_dict

    def get_data_cdf_GCMS(self, f, key):
        '''
        Extracts data from a .cdf file using the Dataset function
        from the netCDF4 library.
        Parameters
        ----------
        f: str
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

        f = Dataset(f, "r")
        f.set_auto_mask(False)
        output = f.variables[key][:]
        f.close()
        return output

    def MS_Load(self):
        file = self.initialised_path

        self.mass_spectra = True
        # Measured intensities
        self.mass_intensity = self.get_data_cdf_GCMS(file, "intensity_values")
        # Measured masses
        self.mass_values = np.round(self.get_data_cdf_GCMS(file, "mass_values"),3)
        # Scan index is starting index
        self.scan_indices = self.get_data_cdf_GCMS(file, "scan_index")
        # Point count is number of elements to read
        self.point_counts = self.get_data_cdf_GCMS(file, "point_count")

    def load_from_csv(self,file):
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
        data = []
        with open(file, 'r') as f:

            for c,line in enumerate(f):
                if c == 0:
                    pass
                else:
                    data.append([float(x) for x in line.strip('\n').split(',')])

        data = np.array(data)

        time = data[:,0]
        signal = data[:,1]

        return time, signal

    def write_to_csv(self, filename = ''):
        '''
        Write chromatgram to a .csv file.
        '''
        if filename == '':
            filename = self.filename

        time   = self.time
        signal = self.signal
        chrom_out = np.vstack((time,signal))
        chrom_out = chrom_out.T

        with open(f'{filename}.csv', 'w') as f:
            f.write('time/ min, signal/ total ion counts\n')
            for x in range(0,len(chrom_out)):
                f.write(f'{chrom_out[x,0]},{chrom_out[x,1]}\n')

    def write_peak_table(self, filename = 'Peak_table',
                        value = 0.0,
                        series_unit = "None"):
        '''
        For writing peak integrals from a chromatogram to a .csv file.

        Parameters
        ----------
        filename: str
            Name for the file
        value: float, str, bool
            Value assigned to the variable assigned to the chromatogram.
        series_unit: float, str, bool
            Units of the variable assigned to the chromatogram.
        Returns
        -------
        None
        '''
        with open("{}.csv".format(filename), "w") as f:
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
        with open('{}.csv'.format(filename), 'w') as f:

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
        if len(self.scan_indices) == 0:
            return {}
        else:

            ion_dict = {np.average(c):np.zeros(len(self.time)) for c in clusters}
            cluster_dict = {np.average(c):c for c in clusters}

            scan_brackets = []

            for s in range(0,len(self.scan_indices)-1):
                scan_brackets.append([self.scan_indices[s],chromatogram.scan_indices[s+1]])

            for s in range(0,len(scan_brackets)):
                inten = self.mass_intensity[scan_brackets[s][0]:scan_brackets[s][1]]
                masses = self.mass_values[scan_brackets[s][0]:scan_brackets[s][1]]

                for m in range(0,len(masses)):
                    for c in clusters:
                        if masses[m] in c:
                            ion_dict[np.average(c)][s] = inten[m]
                            break

            return ion_dict

class Chromatogram_Series:
    def __init__(self,chromatogram_list, information_file):
        '''
        Initialises a series of chromatograms from a list of chromatogram
        objects and an information file.
        '''

        self.chromatograms = chromatogram_list
        self.regions = False

        '''Information provided in file'''
        self.internal_ref_concentration = 1
        self.internal_ref_region = False

        with open(information_file, 'r') as f:
            for line in f:
                if "Dataset" in line:
                    ins = line.strip("\n")
                    self.set_name = ins.split(",")[1]
                if "dilution_factor" in line:
                    ins = line.strip("\n")
                    self.dilution = float(ins.split(",")[1])
                if "series_values" in line:
                    ins = line.strip("\n")
                    spl = ins.split(",")
                    self.x_series =  [float(x) for x in spl[1:] if x != ""]
                if "series_unit" in line:
                    ins = line.strip("\n")
                    self.x_name = ins.split(",")[1]
                if "series_regions" in line:
                    ins = line.strip("\n")
                    reg = [float(x) for x in ins.split(",")[1:] if x != ""]
                    self.regions = [reg[x:x+2] for x in range(0, len(reg), 2)]
                if "internal_ref_region" in line:
                    ins = line.strip("\n")
                    reg = [float(x) for x in ins.split(",")[1:] if x != ""]
                    self.internal_ref_region = [reg[x:x+2] for x in range(0, len(reg), 2) if x != ""][0]
                if "internal_ref_concentration" in line:
                    ins = line.strip("\n")
                    self.internal_ref_concentration = float(ins.split(",")[1])

        if len(self.chromatograms) != len(self.x_series):
            print("The number of chromatograms is {} while the length of the series supplied is {}".format(len(self.chromatograms),len(self.x_series)))
            quit()

        for c,t in zip(self.chromatograms,self.x_series):
             setattr(c,"timepoint", t)

        self.chromatograms.sort(key=lambda x: x.timepoint)
        self.x_series.sort()

        if self.regions == False:
            self.regions = [[chromatogram_list[0].time[0],chromatogram_list[0].time[-1]]]

        '''Read conditions'''
        condset = []
        readstate = False
        with open(information_file, "r") as f:
            for c,line in enumerate(f):
                if "start_conditions" in line:
                    readstate = True
                    line = next(f)
                if "end_conditions" in line:
                    readstate = False
                if readstate:
                    newline = line.strip("\n")
                    condset.append([x for x in newline.split(",") if x != ""])
        c_out = {}
        for c in condset:
            c_out[c[0]] = [float(x) for x in c[1:]]

        self.conditions = c_out

        '''Information to be derived from data'''
        self.integral_series = {}
        self.peak_series = {}
        self.ion_series = {}
        self.internal_ref_integrals = []
        self.conc_series = {}
        self.internal_ref_heights = []
        self.deconvoluted_series = {}

class PeakCollectionElement:
    '''
    Information on individual peakss
    '''
    def __init__(self, position, integral, start, end, parent = 'not specified',
                    mass_spectrum = False):
        '''
        Parameters
        ----------
        position, integral, start, end: float
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
        from ChromProcess import calibration_operations as cal_ops
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
        from ChromProcess import Classes

        if type(file) == str:
            file = Path(file)

        self.filename = file.name

        read_line = lambda line: [float(x) for x in line.strip('\n').split(",") if x != '']
        peaks = []
        IS_pos = 0.0
        IS_integral = 0.0
        IS_bound = [0.0,0.0]
        IS_line_num = -1
        IS = Classes.PeakCollectionElement(0.0, 1, 0.0, 0.0)

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

        is_rt = self.internal_standard.retention_time

        for p in self.peaks:
            p.retention_time = p.retention_time - is_rt + IS_set
            p.start = p.start - is_rt + IS_set
            p.end = p.end - is_rt + IS_set

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
        for pk in self.peaks:
            if pk.assignment in calibrations.calibration_factors:
                CF = calibrations.calibration_factors[pk.assignment]
                if calibrations.calibration_model == 'linear':
                    pk.apply_linear_calibration(CF['C'], CF['B'],
                                                internal_standard = IS_conc)
                elif calibrations.calibration_model == 'quadratic':
                    pk.apply_quadratic_calibration(CF['A'], CF['B'], CF['C'],
                                                   internal_standard = IS_conc)

    def calculate_conc_errors(self, calibrations,IS_conc,IS_conc_err):
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

        if directory == '':
            fname = "{}".format(self.filename)
        else:
            fname = directory/"{}".format(self.filename)

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

        self.name = name
        self.peak_collections = peak_collections
        self.series_values = [pt.series_value for pt in peak_collections]
        self.series_unit = peak_collections[0].series_unit
        self.conditions = conditions
        self.clusters = []

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
        from ChromProcess import processing_functions as p_f
        self.cluster_assignments = []
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

    def make_integral_series(self, cluster_bound = 0.025):
        from ChromProcess import simple_functions as s_f
        from ChromProcess import processing_functions as p_f

        if len(self.clusters) == 0:
            self.get_peak_clusters(self, bound = cluster_bound)

        series_courses = np.zeros((len(self.series_values), len(self.clusters)))

        for c1,pc in enumerate(self.peak_collections):
            for c2,clust in enumerate(self.clusters):
                for pk in pc.peaks:
                    if pk.retention_time == pc.internal_standard.retention_time:
                        continue
                    if pk.retention_time in clust and pk.integral:
                        series_courses[c1,c2] += pk.integral

        self.integral_series = series_courses.T

    def make_concentration_series(self, cluster_bound = 0.025):
        from ChromProcess import simple_functions as s_f
        from ChromProcess import processing_functions as p_f

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
        from ChromProcess import info_params

        conc_dict = {}
        for x in range(0,len(self.concentration_series)):
            name = self.cluster_assignments[x]
            pos = round(np.mean(self.clusters[x]),3)
            if name in info_params.canonical_SMILES:
                smiles = info_params.canonical_SMILES[name.split(' ')[0]]
                conc_dict['{}/ M ({})'.format(smiles, pos)] = self.concentration_series[x]

        return conc_dict

    def integral_traces_as_dict(self):
        from ChromProcess import info_params
        integral_dict = {}
        for x in range(0,len(self.integral_series)):
            name = self.cluster_assignments[x]
            if name in info_params.canonical_SMILES:
                smiles = info_params.canonical_SMILES[name.split(' ')[0]]
                pos = np.mean(self.clusters[x])
                token = smiles + ' ({})'.format(np.round(pos,3))
            else:
                token = name

            integral_dict[token] = self.integral_series[x]

        return integral_dict

    def concentration_error_traces_dict(self):
        from ChromProcess import info_params

        err_dict = {}
        for x in range(0,len(self.conc_err_series)):
            name = self.cluster_assignments[x]
            pos = round(np.mean(self.clusters[x]),3)
            if name in info_params.canonical_SMILES:
                smiles = info_params.canonical_SMILES[name.split(' ')[0]]
                err_dict['{}/ M ({})'.format(smiles, pos)] = self.conc_err_series[x]

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

    def write_concentrations_to_file(self, filename, information):
        '''
        Parameters
        ----------
        filename: name for file including path

        information: ChromProcess Instrument_Calibration object
        '''
        import numpy as np
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

    def write_integrals_to_file(self, filename, information):
        '''
        Parameters
        ----------
        filename: name for file including path

        information: ChromProcess Instrument_Calibration object
        '''
        import numpy as np
        out_type = 'integral_report'
        fname = '{}_{}_{}.csv'.format(filename, information.type, out_type)

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

    def create_conc_DataReport(self,information):
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
        self.filename = fname
        self.mz = mz
        self.relative_abundances = inten
        self.retention_time = pos

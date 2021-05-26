from ChromProcess import file_import
import numpy as np

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

        self.internal_reference = False

        if self.filetype == 'txt':
            data  = file_import.get_data_Shimadzu_HPLC(self.initialised_path,
                                     file_import.get_info_Shimadzu_HPLC(file))

            try:
                chan_ind = data["Wavelength"].index(channel_select)
                self.time = np.array(data["Time"][chan_ind])
                self.signal = np.array(data["Signal"][chan_ind])
            except:
                print("360nm not found, try: ", data["Wavelength"])
                quit()

            self.c_type = 'HPLC'
            self.mass_spectra = False

        if self.filetype == 'cdf':
            self.time   = file_import.get_data_cdf_GCMS(self.initialised_path,
                                                        'scan_acquisition_time')/60 # converted to minutes
            self.signal = file_import.get_data_cdf_GCMS(self.initialised_path,
                                                        'total_intensity')
            self.c_type = 'GCMS'

            if mass_spec == True:
                self.MS_Load()

            else:

                self.mass_spectra = False

                self.mass_values = []
                self.mass_intensity = []
                self.scan_indices = []
                self.point_counts = []
        else:
            print('Unexpected file type ({}), check file in folder. Quitting.'.format(self.filetype))

    def MS_Load(self):
        file = self.initialised_path

        self.mass_spectra = True
        # Measured intensities
        self.mass_intensity = file_import.get_data_cdf_GCMS(file, "intensity_values")
        # Measured masses
        self.mass_values = np.round(file_import.get_data_cdf_GCMS(file, "mass_values"),3)
        # Scan index is starting index
        self.scan_indices = file_import.get_data_cdf_GCMS(file, "scan_index")
        # Point count is number of elements to read
        self.point_counts = file_import.get_data_cdf_GCMS(file, "point_count")


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

class Information:
    def __init__(self, information_file):
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

class Analysis_Information:
    '''
    Stores information about analysis procedure for chromatography data.
    '''
    def __init__(self, regions= [[0,1]], IS_region = [0,1], use_MS = False,
                 analysis_type = 'not specified', i_thres_MS = 500,
                 peak_pick_thres = 0.1, exp_name = 'not specified',
                 dilution_factor = 1.0, dil_err = 0.0,
                 IS_conc = 1.0, IS_err = 0.0,
                 information_file = ''):

        self.experiment_code = exp_name
        self.analysis_type = analysis_type
        self.regions = regions
        self.internal_ref_region = IS_region
        self.use_MS = use_MS
        self.MS_cutoff = i_thres_MS
        self.peak_pick_threshold = peak_pick_thres
        self.dilution_factor = dilution_factor
        self.dilution_factor_error = dil_err
        self.internal_ref_concentration = IS_conc
        self.internal_ref_concentration_error = IS_err

        if information_file == '':
            pass
        else:
            self.read_from_file(information_file)

    def write_to_file(self, directory = ''):
        '''
        directory: str or pathlib Path
        '''
        if directory == '':
            fname = '{}_analysis_details.csv'.format(self.experiment_code)
        else:
            fname = directory/'{}_analysis_details.csv'.format(self.experiment_code)

        with open(fname, 'w') as f:
            f.write('Dataset,{}\n'.format(self.experiment_code))
            f.write('Method,{}\n'.format(self.analysis_type))
            f.write('regions,')
            for r in self.regions:
                for pos in r:
                    f.write('{},'.format(pos))
            f.write('\n')

            f.write('internal_reference_region,')
            for r in self.internal_ref_region:
                f.write('{},'.format(r))
            f.write('\n')
            f.write('extract_mass_spectra,{}\n'.format(str(self.use_MS)))
            f.write('mass_spectra_filter,{}\n'.format(self.MS_cutoff))
            f.write('peak_pick_threshold,{}\n'.format(self.peak_pick_threshold))
            f.write('dilution_factor,{}\n'.format(self.dilution_factor))
            f.write('dilution_factor_error,{}\n'.format(self.dilution_factor_error))
            f.write('internal_ref_concentration,{}\n'.format(self.internal_ref_concentration))
            f.write('internal_ref_concentration_error,{}\n'.format(self.internal_ref_concentration_error))

    def read_from_file(self,fname):
        rdlin = lambda x : [e for e in x.strip('\n').split(',') if e != '']
        with open(fname, 'r') as f:
            for line in f:
                if 'Dataset' in line:
                    ins = rdlin(line)
                    self.experiment_code = ins[1]

                if 'Method' in line:
                    ins = rdlin(line)
                    self.analysis_type = ins[1]

                if 'regions' in line:
                    ins = rdlin(line)
                    reg = [float(x) for x in ins[1:]]
                    self.regions = [reg[x:x+2] for x in range(0, len(reg), 2)]

                if 'internal_reference_region' in line:
                    ins = rdlin(line)
                    reg = [float(x) for x in ins[1:]]
                    self.internal_ref_region = reg

                if 'extract_mass_spectra' in line:
                    ins = rdlin(line)
                    if ins[1] == 'True':
                        self.use_MS = True
                    if ins[1] == 'TRUE':
                        self.use_MS = True
                    if ins[1] == 'False':
                        self.use_MS = False
                    if ins[1] == 'FALSE':
                        self.use_MS = False

                if 'mass_spectra_filter' in line:
                    ins = rdlin(line)
                    self.MS_cutoff = float(ins[1])

                if 'peak_pick_threshold' in line:
                    ins = rdlin(line)
                    self.peak_pick_threshold = float(ins[1])

                if 'dilution_factor,' in line:
                    ins = rdlin(line)
                    self.dilution_factor = float(ins[1])

                if 'dilution_factor_error' in line:
                    ins = rdlin(line)
                    self.dilution_factor_error = float(ins[1])

                if "internal_ref_concentration," in line:
                    ins = rdlin(line)
                    self.internal_ref_concentration = float(ins[1])

                if 'internal_ref_concentration_error' in line:
                    ins = rdlin(line)
                    self.internal_ref_concentration_error = float(ins[1])

class Experiment_Conditions:
    def __init__(self, information_file = ''):

        self.experiment_code = 'not specified'
        self.series_values = []
        self.series_unit = 'not specified'
        self.conditions = {}

        if information_file == '':
            pass
        else:
            self.read_from_file(information_file)

    def write_to_file(self, directory = ''):
        if directory == '':
            fname = '{}_conditions.csv'.format(self.experiment_code)
        else:
            fname = directory/'{}_conditions.csv'.format(self.experiment_code)

        with open(fname, 'w') as f:
            f.write('Dataset,{}\n'.format(self.experiment_code))
            f.write('start_experiment_information\n')

            f.write('series_values,')
            for s in self.series_values:
                f.write('{},'.format(s))

            f.write('\n')
            f.write('series_unit,{}\n'.format(self.series_unit))

            f.write('end_experiment_information\n')
            f.write("start_conditions\n")
            for c in self.conditions:
                f.write("{},".format(c))
                [f.write("{},".format(x)) for x in self.conditions[c]]
                f.write("\n")
            f.write("end_conditions\n")

    def read_from_file(self,information_file):
        with open(information_file, 'r') as f:
            for line in f:
                if "Dataset" in line:
                    ins = line.strip("\n")
                    self.experiment_code = ins.split(",")[1]
                if "series_values" in line:
                    ins = line.strip("\n")
                    spl = ins.split(",")
                    self.series_values =  [float(x) for x in spl[1:] if x != ""]
                if "series_unit" in line:
                    ins = line.strip("\n")
                    self.series_unit = ins.split(",")[1]

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

class Calibration_File:
    '''
    Contains information from a calibration file.
    '''
    def __init__(self,file, type = "None"):

        if type == "None":
            print("Specify calibration type (GCMS or HPLC)")
            print("quitting")
            quit()

        self.calibrations = {}
        self.boundaries = {}
        self.type = type
        info = []
        with open(file, "r") as f:
            readstate = False
            for line in f:
                if "Date_{}".format(type) in line:
                    ins = line.strip("\n").split(",")
                    self.date = ins[1]
                if "Method_file_{}".format(type) in line:
                    ins = line.strip("\n").split(",")
                    self.method = ins[1]
                if "Derivatisation_details_{}".format(type) in line:
                    ins = line.strip("\n").split(",")
                    self.derivatisation = ins[1]

                if "start_calibration" in line:
                    readstate = True
                elif "compound_name" in line:
                    header = line.strip("\n").split(",")
                elif "end_calibration" in line:
                    break
                elif readstate:
                    info.append(line.strip("\n").split(","))

        info = [list(x) for x in zip(*info)]

        comp_ind = header.index('compound_name')

        if type == "GCMS":
            A_ind = header.index('A(GCMS)')
            B_ind = header.index('B(GCMS)')
            C_ind = header.index('C(GCMS)')
            lower_ind = header.index('lower(GCMS)')
            upper_ind = header.index('upper(GCMS)')

        elif type == "HPLC":
            A_ind = header.index('A(HPLC)')
            B_ind = header.index('B(HPLC)')
            C_ind = header.index('C(HPLC)')
            lower_ind = header.index('lower(HPLC)')
            upper_ind = header.index('upper(HPLC)')

        bound_tuples = [(l,u) for l,u in zip(info[lower_ind],info[upper_ind])]
        self.boundaries = {k:[float(v[0]),float(v[1])] for k,v in zip(info[comp_ind],bound_tuples) if v != ("None","None")}

        calib_tuples = [(a,b,c) for a,b,c in zip(info[A_ind],info[B_ind],info[C_ind])]
        self.calibrations = {k:{"A": float(v[0]), "B": float(v[1]), "C": float(v[2])} for k,v in zip(info[comp_ind],calib_tuples) if v != ("None","None","None")}

    def modify_boundaries(self, modified_bounds):

        temp_dict = {}
        for m in modified_bounds:
            if len(modified_bounds[m]) != 2:
                pass
            else:
                temp_dict[m] = modified_bounds[m]

        for i in self.boundaries:
            if i in temp_dict:
                pass
            else:
                temp_dict[i] = self.boundaries[i]

        self.boundaries = {}
        for t in temp_dict:
            self.boundaries[t] = temp_dict[t]


class CalibrationAllocations:
    '''
    Contains information about which calibration file is assigned to which data
    set.
    '''
    def __init__(self,filename):

        self.experiments = []
        self.GCMS_allocations = {}
        self.HPLC_allocations = {}

        if filename == '':
            pass
        else:
            self.import_file(filename)

    def import_file(self,filename):
        rdln = lambda x : [e for e in x.strip('\n').split(',') if e != '']
        assigns = []
        with open(filename, 'r') as f:
            for c, line in enumerate(f):
                if c == 0:
                    header = rdln(line)
                    GCMS_idx = header.index('GCMS_Calibration_file')
                    HPLC_idx = header.index('HPLC_Calibration_file')
                else:
                    ins = rdln(line)
                    self.experiments.append(ins[0])
                    self.GCMS_allocations[ins[0]] = ins[GCMS_idx]
                    self.HPLC_allocations[ins[0]] = ins[HPLC_idx]


class DataPath:
    '''
    Contains information paths to data files.
    '''
    def __init__(self, exp_code, data_type, path):
        self.experiment_code = exp_code
        self.data_type = data_type
        self.path = path

class DataPaths:
    '''
    Container for storing paths to data.
    '''

    def __init__(self, filename):
        from pathlib import Path
        from ChromProcess import Classes

        self.experiment_codes = []
        self.paths = []
        self.exp_code_path = {}
        with open(filename, 'r') as f:
            for c,line in enumerate(f):
                if c == 0:
                    header = line.strip('\n').split(',')
                else:
                    ins = line.strip('\n').split(',')
                    self.experiment_codes.append(ins[0])
                    stem = Path(ins[1])/ins[2]
                    directory = stem/ins[3]
                    dat_path = Classes.DataPath(ins[0], ins[2], directory)
                    self.paths.append(dat_path)
                    self.exp_code_path[ins[0]] = dat_path

class PeakCollectionElement:
    '''
    Information on individual peakss
    '''
    def __init__(self, position, integral, start, end, mass_spectrum = False):
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
        read_line = lambda line: [float(x) for x in line.strip('\n').split(",") if x != '']
        peaks = []
        IS_pos = 0.0
        IS_integral = 0.0
        IS_bound = [0.0,0.0]

        with open(file, "r") as f:
            for c,line in enumerate(f):
                if 'None' in line:
                    pass
                elif c == 0:
                    read = [x for x in line.strip('\n').split(',') if x != '']
                    variable = read[0]
                    value = float(read[1])
                elif c == 2:
                    read = read_line(line)
                    IS_pos = read[0]
                    IS_integral = read[1]
                    IS_bound = read[2:]
                elif c < 4:
                    pass
                else:
                    rd = read_line(line)
                    peaks.append(Classes.PeakCollectionElement(round(rd[0],3), rd[1],
                                                               round(rd[2],3),
                                                               round(rd[3],3)))

        IS = Classes.PeakCollectionElement(round(IS_pos,3), IS_integral,
                                           round(IS_bound[0],3),
                                           round(IS_bound[1],3))

        self.filename = file
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
        threshold: float (from 0.0 to 1.0)
        '''
        del_idx = []
        for c,pk in enumerate(self.peaks):
            if pk.integral < threshold:
                del_idx.append(c)

        for d in sorted(del_idx, reverse=True):
            del self.peaks[d]

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

class DataReport:
    def __init__(self, file = ''):
        self.filename = 'not specified'
        self.experiment_code = 'not specified'
        self.conditions = {}
        self.analysis_details = {}
        self.series_values = np.array([])
        self.series_unit = 'not specified'
        self.data = {}

        if file == '':
            pass
        else:
            self.read_from_file(file)

    def read_from_file(self, file):
        '''
        Parameters
        ----------
        file: pathlib Path or str
            path to file.
        '''

        if type(file) == str:
            from pathlib import Path
            file_n = Path(file)
        else:
            file_n = file

        spl_lin = lambda x : [e for e in x.strip('\n').split(',') if e != '']

        self.filename = file.stem

        with open(file, 'r') as f:
            for line in f:
                ins = spl_lin(line)
                if 'Dataset' in line:
                    self.experiment_code = ins[1]

        condset = self.import_file_section(file, "start_conditions",
                                           "end_conditions")

        c_out = {}
        for c in condset:
            self.conditions[c[0]] = [float(x) for x in c[1:]]

        dataset = self.import_file_section(file, "start_data", "end_data")

        e = [list(i) for i in zip(*dataset)]
        d_out = {}
        for s in e:
            d_out[s[0]] = np.array([0 if x == 'nan' else float(x) for x in s[1:]])

        self.series_unit = dataset[0][0]
        self.series_values = d_out[self.series_unit]

        del d_out[self.series_unit]

        self.data = d_out

        readstate = False
        analysis = self.import_file_section(file, "start_analysis_details",
                                            "end_analysis_details")
        for a in analysis:
            self.analysis_details[a[0]] = [x for x in a[1:]]

    def import_file_section(self, file, start_token, end_token):
        '''
        Parameters
        ----------
        file: str or pathlib Path
            path to file

        start_token: str
            String in line to start reading file from.
        end_token:
            String in line to end reading file from.
        '''

        spl_lin = lambda x : [e for e in x.strip('\n').split(',') if e != '']
        readstate = False
        c_set = []
        with open(file, 'r') as f:
            for c,line in enumerate(f):
                if start_token in line:
                    readstate = True
                    line = next(f)
                if end_token in line:
                    readstate = False
                if readstate:
                    newline = spl_lin(line)
                    c_set.append(newline)

        return c_set

    def write_conditions_header(self, outfile):
        '''
        Parameters
        ----------
        outfile: Python file object
        '''
        # writing experiment conditions to file
        outfile.write("Dataset,{}\n".format(self.experiment_code))
        outfile.write("start_conditions\n")
        for c in self.conditions:
            outfile.write("{},".format(c))
            [outfile.write("{},".format(x)) for x in self.conditions[c]]
            outfile.write("\n")
        outfile.write("end_conditions\n")
        # writing analysis details
        outfile.write("start_analysis_details\n")
        for ad in self.analysis_details:
            outfile.write('{},'.format(ad))
            if type(self.analysis_details[ad]) == str:
                outfile.write('{},'.format(self.analysis_details[ad]))
            else:
                [outfile.write('{},'.format(x)) for x in self.analysis_details[ad]]
            outfile.write('\n')
        outfile.write("end_analysis_details\n")

    def write_to_file(self, filename = '', path = None):
        '''
        Parameters
        ----------
        filename: str
            name for file
        path: pathlib Path object
            Path to folder for file storage.
        '''

        import numpy as np
        an_type = self.analysis_details['Chromatography_method'][0]

        if filename == '':
            filename = self.filename
        elif not filename.endswith('csv'):
            filename = filename + 'csv'
        if path == None:
            fname = filename
        else:
            fname = path/filename

        with open(fname, 'w') as outfile:
            # writing experiment conditions to file
            self.write_conditions_header(outfile)
            # writing data
            sorted_keys = sorted([*self.data], key = lambda x:x.count('C'))

            outfile.write("start_data\n")

            p_header = [self.series_unit]
            out = np.array([self.series_values])

            for s in sorted_keys:
                p_header.append(s)
                out = np.vstack((out,self.data[s]))

            out = out.T
            [outfile.write("{},".format(x)) for x in p_header]

            outfile.write("\n")

            for x in range(0,len(out)):
                for y in range(0,len(out[x])):
                    outfile.write("{},".format(out[x,y]))
                outfile.write("\n")

            outfile.write("end_data\n")

            outfile.write('start_errors\n')
            if len(self.errors)> 0:
                err_out = np.array([self.series_values])

                for s in sorted_keys:
                    err_out = np.vstack((err_out,self.errors[s]))

                err_out = err_out.T
                [outfile.write("{},".format(x)) for x in p_header]
                outfile.write("\n")
                for x in range(0,len(err_out)):
                    for y in range(0,len(err_out[x])):
                        outfile.write("{},".format(err_out[x,y]))
                    outfile.write("\n")

            outfile.write('end_errors\n')

    def find_repeat_data_entries(self):
        entries = []
        repeat_entries = []
        for d in self.data:
            token = d.split(' ')[0]
            if token in entries:
                repeat_entries.append(token)
            entries.append(token)

        return list(set(repeat_entries))

    def remove_repeat_entries(self):
        import numpy as np
        # deleting duplicate entries: taking the entry with the higher signal using the
        # signal sum as a discriminant.
        repeat_entries = self.find_repeat_data_entries()

        for r in repeat_entries:
            compare_keys = []
            for d in self.data:
                if r in d:
                    compare_keys.append(d)

            checkline = np.zeros(len(compare_keys))
            for c,comp in enumerate(compare_keys):
                checkline[c] = np.sum(self.data[comp])

            i_min = np.argmin(checkline)

            del self.data[compare_keys[i_min]]

    def remove_specific_entries(self,remove_list):
        '''
        Parameters
        ----------
        remove_list: list
            List of entries to remove from self.data
        '''
        for r in remove_list:
            del self.data[r]

    def remove_entries_below_threshold(self, threshold):
        '''
        Parameters
        ----------
        threshold: float
            threshold below which entries will be removed.
        '''
        import numpy as np
        # remove entries whose concentrations/integrals do not cross a defined boundary
        del_list = []
        for d in self.data:
            if np.amax(self.data[d]) < threshold:
                del_list.append(d)

        self.remove_specific_entries(del_list)

class MassSpectrum:
    def __init__(self, fname, mz, inten, pos = None):
        self.filename = fname
        self.mz = mz
        self.relative_abundances = inten
        self.retention_time = pos

class Instrument_Calibration:
    def __init__(self, file = ''):

        self.filename = file
        self.calibration_factors = {}
        self.boundaries = {}
        self.type = 'not specified'
        self.date = 'not specified'
        self.method = 'not specified'
        self.derivatisation = 'not specified'
        self.calibration_model = 'not specified'
        self.instrument = 'not specified'
        self.boundaries = {}
        self.calibration_factors = {}
        self.internal_standard_position = False
        self.modified_bounds = {}

        if file == '':
            pass
        else:
            self.read_calibration_file(file)

    def read_calibration_file(self, fname):
        rdlin = lambda x : [e for e in x.strip('\n').split(',') if e != '']
        info = []
        with open(fname, "r") as f:
            readstate = False
            for line in f:
                ins = rdlin(line)
                if 'Technique' in line:
                    self.type = ins[1]
                if 'Instrument' in line:
                    self.instrument = ins[1]
                if "Date" in line:
                    self.date = ins[1]
                if "Method_file" in line:
                    self.method = ins[1]
                if "Derivatisation_details" in line:
                    self.derivatisation = ';'.join(ins[1:])
                if 'calibration_model' in line:
                    self.calibration_model = ins[1]
                if 'internal_standard_position' in line:
                    self.internal_standard_position = float(ins[1])
                if "Instrument" in line:
                    self.instrument = ins[1]
                if "start_calibration" in line:
                    readstate = True
                elif "compound_name" in line:
                    header = line.strip("\n").split(",")
                elif "end_calibration" in line:
                    break
                elif readstate:
                    info.append(ins)

        calib_header = ['A','A_variance','B','B_variance','C', 'C_variance',
                        'AB_covariance', 'AC_covariance',
                        'BC_covariance']

        comp_ind = header.index('compound_name')
        lower_ind = header.index('lower')
        upper_ind = header.index('upper')

        for v in info:
            if 'None' not in v:
                self.calibration_factors[v[comp_ind]] = {}
                for c in calib_header:
                    idx = header.index(c)
                    self.calibration_factors[v[comp_ind]][c] = float(v[idx])

        trans_info = [list(x) for x in zip(*info)]

        bound_tuples = [(l,u) for l,u in zip(trans_info[lower_ind],trans_info[upper_ind])]
        self.boundaries = {k:[float(v[0]),float(v[1])] for k,v in zip(trans_info[comp_ind],bound_tuples) if v != ("None","None")}

    def get_info(self):
        return {'date':self.date, 'method': self.method,
                'derivatisation details': self.derivatisation,
                'calibration model': self.calibration_model,
                'calibration type':self.type}

    def modify_boundaries(self, modified_bounds):
        '''
        Parameters
        ----------
        modified_bounds:
            Dictionary of modifications to make to self.boundaries
        '''

        self.modified_bounds = modified_bounds

        if len(modified_bounds) == 0:
            pass
        else:
            temp_dict = {}
            for m in modified_bounds:
                if len(modified_bounds[m]) == 2:
                    temp_dict[m] = modified_bounds[m]

            for i in self.boundaries:
                if i in temp_dict:
                    pass
                else:
                    temp_dict[i] = self.boundaries[i]

            self.boundaries = {}
            for t in temp_dict:
                self.boundaries[t] = temp_dict[t]

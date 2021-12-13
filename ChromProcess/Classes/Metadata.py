import numpy as np
from ChromProcess.simple_functions import isfloat

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
        with open(information_file, 'r', encoding = 'latin-1') as f:
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
        with open(information_file, "r", encoding = 'latin-1') as f:
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
            c_out[c[0]] = []
            for x in c[1:]:
                if isfloat(x):
                    c_out[c[0]].append(float(x))
                else:
                    c_out[c[0]].append(x)

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
        header = []
        upper_ind = 0
        lower_ind = 0
        A_ind = 0
        B_ind = 0
        C_ind = 0
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
        self.GC_allocations = {}

        if filename == '':
            pass
        else:
            self.import_file(filename)

    def import_file(self,filename):
        rdln = lambda x : [e for e in x.strip('\n').split(',') if e != '']
        GCMS_idx = 0
        HPLC_idx = 0
        GC_idx = 0
        with open(filename, 'r') as f:
            for c, line in enumerate(f):
                if c == 0:
                    header = rdln(line)
                    GCMS_idx = header.index('GCMS_Calibration_file')
                    HPLC_idx = header.index('HPLC_Calibration_file')
                    GC_idx = header.index('GC_Calibration_file')
                else:
                    ins = rdln(line)
                    self.experiments.append(ins[0])
                    self.GCMS_allocations[ins[0]] = ins[GCMS_idx]
                    self.HPLC_allocations[ins[0]] = ins[HPLC_idx]
                    self.GC_allocations[ins[0]] = ins[GC_idx]

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
                    pass
                else:
                    ins = line.strip('\n').split(',')
                    self.experiment_codes.append(ins[0])
                    stem = Path(ins[1])/ins[2]
                    directory = stem/ins[3]
                    dat_path = Classes.DataPath(ins[0], ins[2], directory)
                    self.paths.append(dat_path)
                    self.exp_code_path[ins[0]] = dat_path

class DataReport:
    def __init__(self, file = ''):
        self.filename = 'not specified'
        self.experiment_code = 'not specified'
        self.conditions = {}
        self.analysis_details = {}
        self.series_values = np.array([])
        self.errors = np.array([])
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

        spl_lin = lambda x : [e for e in x.strip('\n').split(',') if e != '']

        self.filename = file.stem

        with open(file, 'r') as f:
            for line in f:
                ins = spl_lin(line)
                if 'Dataset' in line:
                    self.experiment_code = ins[1]

        condset = self.import_file_section(file, "start_conditions",
                                           "end_conditions")

        for c in condset:
            self.conditions[c[0]] = []
            for x in c[1:]:
                if isfloat(x):
                    self.conditions[c[0]].append(float(x))
                else:
                    self.conditions[c[0]].append(x)

        dataset = self.import_file_section(file, "start_data", "end_data")

        e = [list(i) for i in zip(*dataset)]
        d_out = {}
        for s in e:
            d_out[s[0]] = np.array([0 if x == 'nan' else float(x) for x in s[1:]])

        self.series_unit = dataset[0][0]
        self.series_values = d_out[self.series_unit]

        del d_out[self.series_unit]

        self.data = d_out

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
            for _,line in enumerate(f):
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

class Instrument_Calibration:
    def __init__(self, file = 'output.csv'):
        '''
        Parameters
        ----------
        filename: str or pathlib Path
            name for file
        '''

        from pathlib import Path

        if isinstance(file, str):
            file = Path(file)

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
        header = []
        upper_ind = 0
        lower_ind = 0
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
        self.boundaries = {
                             k: [
                                float(v[0]),
                                float(v[1])
                                ] for k,v in zip(trans_info[comp_ind],bound_tuples)
                                                                 if v != ("None","None")}

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

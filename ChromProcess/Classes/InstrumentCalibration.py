
class InstrumentCalibration:
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
        '''
        Parameters
        ----------
        fname: str or pathlib Path
            name for file
        '''

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
        modified_bounds: dict
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

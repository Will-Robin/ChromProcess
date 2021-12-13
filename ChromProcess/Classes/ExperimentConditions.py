
class ExperimentConditions:
    def __init__(self, information_file = ''):
        '''
        Parameters
        ----------
        information_file: str or pathlib Path
        '''

        self.experiment_code = 'not specified'
        self.series_values = []
        self.series_unit = 'not specified'
        self.conditions = {}

        if information_file == '':
            pass
        else:
            self.read_from_file(information_file)

    def write_to_file(self, directory = ''):
        '''
        Parameters
        ----------
        directory: str or pathlib Path
        '''

        from pathlib import Path

        fname = self.experiment_code
        if isinstance(directory, str):
            if directory == '':
                fname = '{}_conditions.csv'.format(self.experiment_code)
            else:
                fname = f'{directory}/{self.experiment_code}_conditions.csv'

        elif isinstance(directory, Path):
            fname = directory/f'{self.experiment_code}_conditions.csv'
        else:
            pass

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
        '''
        Parameters
        ----------
        information_file: str or pathlib Path
        '''

        from ChromProcess import simple_functions as s_f

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

        # Read conditions
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
                if s_f.isfloat(x):
                    c_out[c[0]].append(float(x))
                else:
                    c_out[c[0]].append(x)

        self.conditions = c_out


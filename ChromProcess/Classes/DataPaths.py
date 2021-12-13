
class DataPaths:
    '''
    Container for storing paths to data.
    '''

    def __init__(self, filename):
        '''
        Parameters
        ----------
        filename: str or pathlib Path
        '''

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


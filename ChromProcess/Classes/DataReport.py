
class DataReport:
    def __init__(self, file = ''):
        '''
        Parameters
        ----------
        file: str or pathlib Path
        '''
        import numpy as np

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

        import numpy as np
        from pathlib import Path
        from ChromProcess import simple_functions as s_f

        if isinstance(file, str):
            file = Path(file)

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
                if s_f.isfloat(x):
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
            self.analysis_details[a[0]] = a[1:]

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
            if isinstance(self.analysis_details[ad], str):
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
        from pathlib import Path

        if isinstance(path, str):
            path = Path(path)

        if filename == '':
            filename = self.filename

        if path is None:
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
                    trace = self.errors[s]
                    err_out = np.vstack((err_out,trace))

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

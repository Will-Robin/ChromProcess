
class AnalysisInformation:
    '''
    Stores information about analysis procedure for chromatography data.
    '''
    def __init__(self, regions= [[0,1]], IS_region = [0,1], use_MS = False,
                 analysis_type = 'not specified', i_thres_MS = 500,
                 peak_pick_thres = 0.1, exp_name = 'not specified',
                 dilution_factor = 1.0, dil_err = 0.0,
                 IS_conc = 1.0, IS_err = 0.0,
                 information_file = ''):
        '''
        Parameters
        ----------
        Stores information about analysis procedure for chromatography data.
        regions: list of lists
        IS_region: list
        use_MS: bool
        analysis_type: str
        i_thres_MS: float
        peak_pick_thres: float
        exp_name: float
        dilution_factor: float
        dil_err: float
        IS_conc: float
        IS_err: float
        information_file: str or pathlib Path
        '''

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

    def read_from_file(self, fname):
        '''
        Parameters
        ----------
        fname: str or pathlib Path
        '''

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
                fname = f'{self.experiment_code}_analysis_details.csv'
            else:
                fname = f'{directory}/{self.experiment_code}_analysis_details.csv'
        elif isinstance(directory, Path):
            fname = directory/f'{self.experiment_code}_analysis_details.csv'.format()
        else:
            pass

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


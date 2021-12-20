
class AnalysisInformation:
    '''
    Stores information about analysis procedure for chromatography data.
    '''
    def __init__(self):
        '''
        Attributes
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

        self.experiment_code = ''
        self.analysis_type = ''
        self.regions = []
        self.internal_ref_region = []
        self.use_MS = False
        self.MS_cutoff = 0.0
        self.peak_pick_threshold = 0.0
        self.dilution_factor = 1.0
        self.dilution_factor_error = 0.0
        self.internal_ref_concentration = 0.0
        self.internal_ref_concentration_error = 0.0

    def write_to_file(self, directory = ''):
        '''
        Parameters
        ----------
        directory: str or pathlib Path
        '''
        from pathlib import Path

        exp_code = self.experiment_code
        fname = self.experiment_code
        if isinstance(directory, str):
            if directory == '':
                fname = f'{exp_code}_analysis_details.csv'
            else:
                fname = f'{directory}/{exp_code}_analysis_details.csv'
        elif isinstance(directory, Path):
            fname = directory/f'{exp_code}_analysis_details.csv'
        else:
            pass

        is_conc = self.internal_ref_concentration
        is_conc_err = self.internal_ref_concentration_error
        dil_factor = self.dilution_factor
        dil_factor_err = self.dilution_factor_error
        peak_pick_thresh = self.peak_pick_threshold
        ms_cutoff = self.MS_cutoff
        us_ms = self.use_MS

        with open(fname, 'w') as f:
            f.write(f'Dataset,{exp_code}\n')
            f.write(f'Method,{self.analysis_type}\n')
            f.write(f'regions,')
            for r in self.regions:
                for pos in r:
                    f.write(f'{pos},')
            f.write('\n')

            f.write('internal_reference_region,')
            for r in self.internal_ref_region:
                f.write(f'{r},')
            f.write('\n')
            f.write(f'extract_mass_spectra,{use_ms}\n')
            f.write(f'mass_spectra_filter,{ms_cutoff}\n')
            f.write(f'peak_pick_threshold,{peak_pick_thresh}\n')
            f.write(f'dilution_factor,{dil_factor}\n')
            f.write(f'dilution_factor_error,{dil_factor_err}\n')
            f.write(f'internal_ref_concentration,{is_conc}\n')
            f.write(f'internal_ref_concentration_error,{is_conc_err}\n')


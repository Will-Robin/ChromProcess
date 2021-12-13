import sys

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


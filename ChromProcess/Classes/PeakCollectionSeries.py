import numpy as np
from pathlib import Path

from ChromProcess.Utils.utils import utils
from ChromProcess.Processing.peak import assign_peak
from ChromProcess.Utils.utils.clustering import cluster

class PeakCollectionSeries:
    def __init__(self, peak_collections, name = 'not specified',
                 conditions = {}):
        '''
        An object which wraps multiple PeakCollection objects.

        Parameters
        ----------
        peak_collections: List of ChromProcess PeakCollection objects
        name: str
        conditions: dict
        '''

        self.name = name
        self.peak_collections = peak_collections
        self.series_values = [pt.series_value for pt in peak_collections]
        self.series_unit = peak_collections[0].series_unit
        self.conditions = conditions
        self.clusters = []
        self.cluster_assignments = []
        self.integral_series = []
        self.concentration_series = []
        self.conc_err_series = []
        self.series_assigned_compounds = []

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

        peak_pos = np.array([])
        for pc in self.peak_collections:
            peak_pos = np.hstack((peak_pos, pc.get_peak_positions()))

        i = np.argsort(peak_pos)

        return peak_pos[i]

    def get_peak_clusters(self, bound = 0.1):
        '''
        bound: float
        '''

        peaks = self.get_peak_positions()

        clusts = []
        for c in cluster(peaks, bound = bound):
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

        for c in self.clusters:
            pos = np.mean(c)
            clust_name = assign_peak.assign_retention_time(pos, boundaries)
            self.cluster_assignments.append(clust_name)

    def get_all_assigned_compounds(self):
        assigns = []
        for pc in self.peak_collections:
            if len(pc.assigned_compounds) == 0:
                pc.get_all_assigned_compounds()
            assigns.extend( pc.assigned_compounds )

        assigns = list(set(assigns))
        count_C = lambda x:x.count('C')

        self.series_assigned_compounds = sorted(assigns, key = count_C)

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

        IS_conc = conditions.internal_standard_concentration

        for pc in self.peak_collections:
            pc.apply_calibrations_to_peaks(calibrations, IS_conc)

    def calculate_conc_errors(self, calib, conditions):
        '''
        Parameters
        ----------
        calib: ChromProcess Instrument_Calibration object
        '''
        IS_conc = conditions.internal_standard_concentration
        IS_conc_err = conditions.internal_standard_concentration_error

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

    def make_integral_series(self, cluster_bound = 0.0):

        '''
        Parameters
        ----------
        cluster_bound: float
        '''

        if len(self.clusters) == 0:
            self.get_peak_clusters(bound = cluster_bound)

        series_courses = np.zeros((len(self.series_values), len(self.clusters)))

        for c1,pc in enumerate(self.peak_collections):
            for c2,clust in enumerate(self.clusters):
                for pk in pc.peaks:
                    if pk.retention_time == pc.internal_standard.retention_time:
                        continue
                    if pk.retention_time in clust and pk.integral:
                        series_courses[c1,c2] += pk.integral

        self.integral_series = series_courses.T

    def make_concentration_series(self, cluster_bound = 0.0):
        '''
        Parameters
        ----------
        cluster_bound: float
        '''

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
        '''
        Parameters
        ----------
        name_conversions: dict
        '''

        if len(self.cluster_assignments) == 0:
            get_name = lambda _: ''
        else:
            get_name = lambda x: self.cluster_assignments[x].split(' ')[0]

        conc_dict = {}
        for x in range(0,len(self.concentration_series)):
            name = get_name(x)
            pos = np.round(np.mean(self.clusters[x]),3)
            if name == '':
                pass
            elif not utils.is_float(name):
                conc_dict[f'{name}/ M ({pos})'] = self.concentration_series[x]
            else:
                pass

        return conc_dict

    def integral_traces_as_dict(self):

        if len(self.cluster_assignments) == 0:
            get_name = lambda _: ''
        else:
            get_name = lambda x: self.cluster_assignments[x].split(' ')[0]

        integral_dict = {}
        for x in range(0,len(self.integral_series)):
            name = get_name(x)
            if name != '':
                pos = np.mean(self.clusters[x])
                token = f'{name} ({np.round(pos,3)})'
            else:
                cluster_average = np.mean(self.clusters[x])
                val = np.round(cluster_average, 3)
                token = f'{val} ({val})'

            integral_dict[token] = self.integral_series[x]

        return integral_dict

    def concentration_error_traces_dict(self):

        err_dict = {}

        if len(self.cluster_assignments) == 0:
            get_name = lambda _: ''
        else:
            get_name = lambda x: self.cluster_assignments[x].split(' ')[0]

        for x in range(0,len(self.conc_err_series)):
            name = get_name(x)
            pos = np.round(np.mean(self.clusters[x]),3)
            if name == '':
                pass
            elif not utils.is_float(name):
                err_dict[f'{name}/ M ({pos})'] = self.conc_err_series[x]
            else:
                pass

        return err_dict

    def write_conditions_header(self, outfile, info):
        '''
        Parameters
        ----------
        outfile: Python file object

        information: ChromProcess Instrument_Calibration object
        '''

        # writing experiment conditions to file
        outfile.write(f"Dataset,{self.name}\n")
        outfile.write("start_conditions\n")
        for c in self.conditions:
            outfile.write(f"{c},")
            [outfile.write(f"{x},") for x in self.conditions[c]]
            outfile.write("\n")
        outfile.write("end_conditions\n")

        # writing analysis details
        outfile.write("start_analysis_details\n")
        outfile.write(f'Instrument, {info.instrument}\n')
        outfile.write(f"Chromatography_method,{info.type},{info.method}\n")
        outfile.write(f"Derivatisation_method,{info.derivatisation}\n")
        outfile.write(f"Calibrations_file,{info.filename}\n")
        outfile.write(f'Calibration_model,{info.calibration_model}\n')
        outfile.write(f"end_analysis_details\n")

    def write_concentrations_to_file(self, filename, info):
        '''
        Parameters
        ----------
        filename: name for file including path

        information: ChromProcess Instrument_Calibration object
        '''


        if isinstance(filename, str):
            filename = filename
        elif isinstance(filename, Path):
            filename = str(filename)

        out_type = 'concentration_report'
        fname = '{}_{}_{}.csv'.format(filename, info.type, out_type)

        with open(fname, 'w') as outfile:
            # writing experiment conditions to file
            self.write_conditions_header(outfile,info)

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

    def write_integrals_to_file(self,
                                filename,
                                information,
                                ):
        '''
        Parameters
        ----------
        filename: name for file including path

        information: ChromProcess Instrument_Calibration object
        '''

        import numpy as np
        from pathlib import Path

        name = ''
        if isinstance(filename, str):
            name = filename
        elif isinstance(filename, Path):
            name = str(filename)

        out_type = 'integral_report'
        fname = '{}_{}_{}.csv'.format(name, information.type, out_type)

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

    def create_DataReport_base(self, info):
        '''
        Parameters
        ----------
        information: ChromProcess Instrument_Calibration object
        '''

        from ChromProcess import Classes

        data_report = Classes.DataReport()
        data_report.experiment_code = self.name
        data_report.conditions = self.conditions
        data_report.analysis_details['Instrument'] = info.instrument
        chrom_method = [info.type, info.method]
        data_report.analysis_details['Chromatography_method'] = chrom_method
        deriv_info = info.derivatisation
        data_report.analysis_details['Derivatisation_method'] = deriv_info
        data_report.analysis_details['Calibrations_file'] = info.filename
        data_report.analysis_details['Calibration_model'] =info.calibration_model

        data_report.series_values = self.series_values
        data_report.series_unit = self.series_unit

        return data_report

    def create_conc_DataReport(self, information):
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

import numpy as np

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

    def apply_peak_dilution_factors(self, dilution_factor, error):
        '''
        Parameters
        ----------
        dilution_factor: float
        error: float
        '''
        for pc in self.peak_collections:
            pc.dilution_correct_peaks(dilution_factor, error)

    def create_peak_series(self):
        '''
        Creates arrays for series of peaks using self.clusters to identify 
        similar peaks between chromatograms. If more than one peak from a
        single chromatogram is in a cluster, their values are added together.

        The ordering of the series in the arrays is given by the order of the
        self.clusters.
        '''

        series_length = len(self.series_values)
        number_of_peaks = len(self.clusters)

        integral_series = np.zeros((series_length, number_of_peaks))
        concentration_series = np.zeros((series_length, number_of_peaks))
        error_series = np.zeros((series_length, number_of_peaks))

        for c1,pc in enumerate(self.peak_collections):
            is_rt = pc.internal_standard.retention_time

            for c2,clust in enumerate(self.clusters):
                for pk in pc.peaks:

                    if pk.retention_time == is_rt:
                        continue

                    if pk.retention_time in clust:
                        if pk.integral:
                            integral_series[c1,c2] += pk.integral

                        if pk.concentration:
                            concentration_series[c1,c2] += pk.concentration

                        if pk.conc_error:
                            error_series[c1,c2] += pk.conc_error

        self.integral_series = integral_series.T
        self.concentration_series = concentration_series.T
        self.conc_err_series = error_series.T

    def series_traces_as_dict(self):
        '''
        Create dictionaries of peak series values derived using the
        self.clusters.
        '''
        conc_dict = {}
        err_dict = {}
        integral_dict = {}

        if len(self.cluster_assignments) == 0:
            cluster_names = ['' for _ in self.clusters]
        else:
            cluster_names = [n.split(' ')[0] for n in self.cluster_assignments]

        for c1,pc in enumerate(self.peak_collections):

            for c2,clust in enumerate(self.clusters):
                name = cluster_names[c2]
                average_position = sum(clust)/len(clust)

                for pk in pc.peaks:
                    if pk.retention_time in clust:
                        if pk.integral:
                            #Prevent unnamed compounds to be saved as (value) 
                            #as this will cause problems when saving as csv.
                            if name == '': 
                                token = f'{average_position}'
                            else:
                                token = name + f' ({average_position})'
                            if token not in integral_dict:
                                integral_dict[token] = [0.0 for _ in self.series_values]
                            integral_dict[token][c1] += pk.integral

                        if pk.concentration:
                            token = name + '/ M'
                            if token not in conc_dict:
                                conc_dict[token] = [0.0 for _ in self.series_values]
                            conc_dict[token][c1] += pk.concentration

                        if pk.conc_error:
                            token = name + '/ M'
                            if token not in err_dict:
                                err_dict[token] = [0.0 for _ in self.series_values]
                            err_dict[token][c1] += pk.conc_error

        return conc_dict, err_dict, integral_dict

    def write_data_reports(self, filename, information, cluster_removal_limit = False):
        '''
        Parameters
        ----------
        filename: name for file including path
        information: ChromProcess Analysis_Information object
        '''
        import ChromProcess.Writers as write

        write.peak_collection_series_to_data_report(self, filename, information, cluster_removal_limit = cluster_removal_limit)

import numpy as np
from ChromProcess import Classes
from ChromProcess.Writers import peak_collection_to_csv


class PeakCollection:
    def __init__(self):
        """
        An object for storing and operating upon collections of peaks.

        Attributes
        ----------
        self.filename: str
            Name of the file that the object was created from.
        self.series_value: float
            Value of the point in the series at which the peaks were measured.
        self.series_unit: str
            Unit for the series values.
        self.internal_standard: Classes.Peak
            The internal standard Peak
        self.peaks: list
            List of peaks.
        self.mass_spectra: list
            List of mass spectra.
        self.initial_IS_pos: float
            The position of the internal standard when the object is created
            from a file.
        self.assigned_compounds: list
            A list of names of assigned compounds in the Peak
        """

        self.filename = "not specified"
        self.series_value = 0.0
        self.series_unit = "not specified"
        self.internal_standard = Classes.Peak(0.0, 0.0, 0.0)
        self.peaks = [Classes.Peak(0.0, 0.0, 0.0)]
        self.mass_spectra = []
        self.initial_IS_pos = 0.0
        self.assigned_compounds = []

    def remove_peaks_below_threshold(self, threshold):
        """
        Remove peaks below a certain integral threshold (note that this
        operates on the values held in the Peak integral attributes. If
        they have been normalised to an internal standard, the
        threshold value should probably be lower than for the 'raw'
        integral data).

        Parameters
        ----------
        threshold: float

        Returns
        -------
        None
        """

        del_idx = []
        for c, pk in enumerate(self.peaks):
            if pk.integral < threshold:
                del_idx.append(c)

        self.peaks = [v for i, v in enumerate(self.peaks) if i not in del_idx]

    def align_peaks_to_IS(self, IS_set=0.0):
        """
        Align the peak integral retention times to the internal standard's
        position, which will be set to IS_set.

        Parameters
        ----------
        IS_set: float
            Position to set the internal standard's retention time.

        Returns
        -------
        None
        """

        is_rt = self.internal_standard.retention_time
        is_start = self.internal_standard.start
        is_end = self.internal_standard.end

        if is_rt == 0.0:
            print(
                """Internal standard retention time is set to 0.0.
            No alignment performed."""
            )
        else:
            for p in self.peaks:
                p.retention_time = p.retention_time - is_rt + IS_set
                p.start = p.start - is_rt + IS_set
                p.end = p.end - is_rt + IS_set

            for m in self.mass_spectra:
                m.retention_time = m.retention_time - is_rt + IS_set

            self.internal_standard.start = is_start - is_rt + IS_set
            self.internal_standard.end = is_end - is_rt + IS_set
            self.internal_standard.retention_time = is_rt - is_rt + IS_set

    def add_mass_spectra(self, ms_list):
        """
        Add mass spectrum information into Peak objects.

        Parameters
        ----------
        ms_list: list of ChromProcess MassSpectrum objects
            Mass spectra to be added.

        Returns
        -------
        None
        """

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
        """
        Get the positions of all of the peaks in the PeakCollection.

        Parameters
        ----------

        Returns
        -------
        array: 2D numpy array
        """
        return np.array([p.retention_time for p in self.peaks])

    def reference_integrals_to_IS(self):
        """
        Divide all peak integrals by the integral of the internal standard
        peak.

        Parameters
        ----------

        Returns
        -------
        None
        """

        IS_integral = self.internal_standard.integral
        for p in self.peaks:
            p.reference_integral_to_IS(IS_integral)

    def reference_heights_to_IS(self):
        """
        Divide all peak heights by the height of the internal standard
        peak.

        Parameters
        ----------

        Returns
        -------
        None
        """

        IS_height = self.internal_standard.height
        for p in self.peaks:
            p.reference_height_to_IS(IS_height)

    def assign_peaks(self, boundaries):
        """
        Assign peaks using boundaries.

        Parameters
        ----------
        boundaries: dict
            {'compound_name': [lower bound, upper bound]}

        Returns
        -------
        None
        """

        for pk in self.peaks:
            pk.assign_peak(boundaries)

    def apply_calibrations_to_peaks(self, calibrations, IS_conc):
        """
        Apply calibrations to peals

        Parameters
        ----------
        calibrations: ChromProcess Instrument_Calibration object
            Container for calibration information
        IS_conc: float
            internal standard concentration

        Returns
        -------
        None
        """

        if IS_conc == 0.0:
            # results in division by 1 during conversion
            IS_conc = 1.0
        for pk in self.peaks:
            if pk.assignment in calibrations.calibration_factors:

                CF = calibrations.calibration_factors[pk.assignment]

                if calibrations.calibration_model == "linear":
                    pk.apply_linear_calibration(
                        CF["B"], CF["C"], internal_standard=IS_conc
                    )

                elif calibrations.calibration_model == "quadratic":
                    pk.apply_quadratic_calibration(
                        CF["A"], CF["B"], CF["C"], internal_standard=IS_conc
                    )
                else:
                    print("Calibration type not recognised")

    def calculate_conc_errors(self, calibrations, IS_conc, IS_conc_err):
        """
        Calculation of the standard error on a concentration estimation from
        a calibration calculation.

        Parameters
        ----------
        calibrations: ChromProcess Instrument_Calibration object
            Contains calibration information.
        IS_conc: float
        IS_conc_err: float

        Returns
        -------
        None
        """

        for pk in self.peaks:
            pk.calculate_error(calibrations, IS_conc, IS_conc_err)

    def dilution_correct_peaks(self, dilution_factor, error):
        """
        Divide peaks by dilution_factor, accounting for error.

        Parameters
        ----------
        dilution_factor: float
        error: float

        Returns
        -------
        None
        """

        for pk in self.peaks:
            if pk.concentration:
                pk.dilution_correction(dilution_factor, error)

    def get_all_assigned_compounds(self):
        """
        Get the names of all of the compounds assigned in the PeakCollection.

        Parameters
        ----------

        Returns
        -------
        self.assigned_compounds: list
            List of names of assigned compounds.
        """
        assigns = []
        for pk in self.peaks:
            if pk.assignment == "none":
                pass
            else:
                assigns.append(pk.assignment)

        assigns = list(set(assigns))

        self.assigned_compounds = sorted(assigns, key=lambda x: x.count("C"))

        return self.assigned_compounds[:]

    def write_to_file(self, directory=""):
        """
        Write the object to a structured file.

        Parameters
        ----------
        directory: str or pathlib Path

        Returns
        -------
        None
        """

        peak_collection_to_csv(self, directory=directory)

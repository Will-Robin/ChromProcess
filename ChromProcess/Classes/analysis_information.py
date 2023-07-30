class AnalysisInformation:
    """
    Stores information about analysis procedure for chromatography data.
    """

    def __init__(self):
        """
        Initialisation of the AnalysisInformation object.

        Attributes
        ----------
        self.experiment_code: str
            Experiment code.
        self.analysis_type: str
            The type of chromatographic methods used.
        self.regions: list[lists]
            List of upper and lower bounds describing the regions of a chromatogram.
        self.internal_standard_region: list
            Upper and lower bounds for the internal standard.
        self.use_MS: bool
            Whether to extract mass spectra information for GC-MS data.
        self.MS_cutoff: float
            The lower threshold for inclusion of mass spectral signals.
        self.peak_pick_threshold: float
            The relative cut off for peak picking.
        self.dilution_factor: float
            Sample dilution factor.
        self.dilution_factor_error: float
            Error on dilution factor.
        self.internal_standard_concentration: float
            Concentration of the internal standard.
        self.internal_standard_concentration_error: float
            Error on the concentration of the internal standard.
        """

        self.experiment_code = ""
        self.analysis_type = ""
        self.derivatisation_method = ""
        self.instrument = ""
        self.instrument_method = ""
        self.calibration_model = ""
        self.calibration_file = ""
        self.regions = []
        self.internal_standard_region = []
        self.use_MS = False
        self.MS_cutoff = 0.0
        self.peak_pick_threshold = 0.0
        self.dilution_factor = 1.0
        self.dilution_factor_error = 0.0
        self.internal_standard_concentration = 0.0
        self.internal_standard_concentration_error = 0.0

    def __str__(self):
        str_repr = f"""
        experiment code: {self.experiment_code}
        analysis type: {self.analysis_type}
        derivatisation method: {self.derivatisation_method}
        instrument: {self.instrument}
        instrument method: {self.instrument_method}
        calibration model: {self.calibration_model}
        calibration file: {self.calibration_file}
        regions: {self.regions}
        internal standard region: {self.internal_standard_region}
        use MS: {self.use_MS}
        MS cutoff: {self.MS_cutoff}
        peak pick threshold: {self.peak_pick_threshold}
        dilution factor: {self.dilution_factor}
        dilution factor error: {self.dilution_factor_error}
        IS concentration: {self.internal_standard_concentration}
        IS concentration_error: {self.internal_standard_concentration_error}
        """
        return str_repr

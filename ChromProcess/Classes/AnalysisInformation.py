class AnalysisInformation:
    """
    Stores information about analysis procedure for chromatography data.
    """

    def __init__(self):
        """
        Stores information about analysis procedure for chromatography data.

        Attributes
        ----------
        self.experiment_code: str
            Experiment code.
        self.analysis_type: str
            The type of chromatographic methods used.
        self.regions: list of lists
            List of upper and lower bounds describing the regions of a chromatogram.
        self.internal_standard_region: list
            Upper and lower bounds for the internal standard.
        self.use_MS: bool
            Whether to extract mass spectra information for GC-MS data.
        self.MS_cutoff: float
            The lower threshold for inclusion of mass spectral signals.
        self.peak_pick_threshold: float
            The relative cut off for peak picking.
        self.dilution_factor = 1.0
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

    def write_to_file(self, directory=""):
        """
        Parameters
        ----------
        directory: str or pathlib Path
            Directory in which the file will be saved in.

        Returns
        -------
        None

        """
        from pathlib import Path

        exp_code = self.experiment_code
        fname = self.experiment_code
        if isinstance(directory, str):
            if directory == "":
                fname = f"{exp_code}_analysis_details.csv"
            else:
                fname = f"{directory}/{exp_code}_analysis_details.csv"
        elif isinstance(directory, Path):
            fname = directory / f"{exp_code}_analysis_details.csv"
        else:
            pass

        is_conc = self.internal_standard_concentration
        is_conc_err = self.internal_standard_concentration_error
        dil_factor = self.dilution_factor
        dil_factor_err = self.dilution_factor_error
        peak_pick_thresh = self.peak_pick_threshold
        ms_cutoff = self.MS_cutoff
        use_ms = self.use_MS

        with open(fname, "w") as f:
            f.write(f"Dataset,{exp_code}\n")
            f.write(f"Method,{self.analysis_type}\n")
            f.write(f"regions,")
            for r in self.regions:
                for pos in r:
                    f.write(f"{pos},")
            f.write("\n")

            f.write("internal_standard_region,")
            for r in self.internal_standard_region:
                f.write(f"{r},")
            f.write("\n")
            f.write(f"extract_mass_spectra,{use_ms}\n")
            f.write(f"mass_spectra_filter,{ms_cutoff}\n")
            f.write(f"peak_pick_threshold,{peak_pick_thresh}\n")
            f.write(f"dilution_factor,{dil_factor}\n")
            f.write(f"dilution_factor_error,{dil_factor_err}\n")
            f.write(f"internal_standard_concentration,{is_conc}\n")
            f.write(f"internal_standard_concentration_error,{is_conc_err}\n")

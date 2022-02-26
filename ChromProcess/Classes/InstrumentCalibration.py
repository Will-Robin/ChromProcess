class InstrumentCalibration:
    def __init__(self):
        """
        A container for calibration information.

        Attributes
        ----------
        calibration_factors: dict
        boundaries: dict
        type: str
        date: str
        method: str
        derivatisation: str
        calibration_model: str
        instrument: str
        calibration_factors: dict
        internal_stanard_position: float
        modified_bounds: dict
        """

        self.filename = "not specified"
        self.calibration_factors = {}
        self.boundaries = {}
        self.type = "not specified"
        self.date = "not specified"
        self.method = "not specified"
        self.derivatisation = "not specified"
        self.calibration_model = "not specified"
        self.instrument = "not specified"
        self.calibration_factors = {}
        self.internal_standard_position = 0.0
        self.modified_bounds = {}

    def get_info(self):
        return {
            "date": self.date,
            "method": self.method,
            "derivatisation details": self.derivatisation,
            "calibration model": self.calibration_model,
            "calibration type": self.type,
        }

    def modify_boundaries(self, modified_bounds):
        """
        Modify the assignment boundaries in the calibration.

        Parameters
        ----------
        modified_bounds: dict
            Dictionary of modifications to make to self.boundaries

        Returns
        -------
        None
        """

        self.modified_bounds = modified_bounds

        if len(modified_bounds) == 0:
            pass
        else:
            temp_dict = {}
            for m in modified_bounds:
                if len(modified_bounds[m]) == 2:
                    temp_dict[m] = modified_bounds[m]

            for i in self.boundaries:
                if i in temp_dict:
                    pass
                else:
                    temp_dict[i] = self.boundaries[i]

            self.boundaries = {}
            for t in temp_dict:
                self.boundaries[t] = temp_dict[t]

from pathlib import Path


class ExperimentConditions:
    """
    A container for experimental information.
    """

    def __init__(self):
        """
        Initialise an empty ExperimentConditions object.

        Attributes
        ----------
        self.experiment_code: str
            Code name for the experiment.
        self.series_values: list[Any]
            List of series values for the experiment.
        self.series_unit: str
            Unit of the series values.
        self.conditions: dict
            Experimental conditions.
        """

        self.experiment_code: str = "not specified"
        self.series_values: list[str | float] = []
        self.series_unit: str = "not specified"
        self.conditions: dict[str, float | str | list[float]] = dict()


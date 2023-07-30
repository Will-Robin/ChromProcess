from typing import Any


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
        self.series_values: list[str, float]
            List of series values for the experiment.
        self.series_unit: str
            Unit of the series values.
        self.conditions: dict[str, Any]
            Experimental conditions.
        """

        self.experiment_code: str = "not specified"
        self.series_values: list[str | float] = []
        self.series_unit: str = "not specified"
        self.conditions: dict[str, Any] = dict()

    def __str__(self):
        str_repr = f"""
        experiment_code: {self.experiment_code}
        series_values: {self.series_values}
        series_unit: {self.series_unit}
        conditions: {self.conditions}
        """
        return str_repr

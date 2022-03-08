import tomli
from ChromProcess import Classes

def conditions_from_toml(filename):
    """
    Create a ExperimentConditions object from a .toml file.


    Parameters
    ----------
    filename: str or pathlib Path

    Returns
    -------
    conditions: Classes.ExperimentConditions
    """

    conditions = Classes.ExperimentConditions()

    with open(filename, "r") as f:
        text = f.read()

    info_dict = tomli.loads(text)

    conditions.experiment_code = info_dict["Dataset"]
    conditions.series_values = info_dict["Series_values"]
    conditions.series_unit = info_dict["Series_unit"]

    conds = {}
    for c in info_dict["conditions"]:
        val = info_dict["conditions"][c][0]
        unit = info_dict["conditions"][c][1]

        conds[f"{c}/ {unit}"] = val

    conditions.conditions = conds

    return conditions


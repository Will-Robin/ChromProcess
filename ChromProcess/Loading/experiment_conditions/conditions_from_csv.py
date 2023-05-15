from ChromProcess.Classes import ExperimentConditions
from ChromProcess.Utils.utils import utils


def conditions_from_csv(filename: str) -> ExperimentConditions:
    """
    Create a ExperimentConditions object from a formatted .csv file.

    Parameters
    ----------
    filename: str or pathlib Path

    Returns
    -------
    conditions: ExperimentConditions
    """

    conditions = ExperimentConditions()

    with open(filename, "r", encoding="latin-1") as f:
        for line in f:
            if "Dataset" in line:
                ins = line.strip("\n")
                conditions.experiment_code = ins.split(",")[1]
            if "Series_values" in line:
                ins = line.strip("\n")
                spl = ins.split(",")
                conditions.series_values = [x for x in spl[1:] if x != ""]
            if "Series_unit" in line:
                ins = line.strip("\n")
                conditions.series_unit = ins.split(",")[1]

    # Read conditions
    condset = []
    readstate = False
    with open(filename, "r", encoding="latin-1") as f:
        for c, line in enumerate(f):
            if "start_conditions" in line:
                readstate = True
                line = next(f)
            if "end_conditions" in line:
                readstate = False
            if readstate:
                newline = line.strip("\n")
                condset.append([x for x in newline.split(",") if x != ""])

    c_out: dict[str, list[float | str]] = dict()
    for cndtn in condset:
        cndtn_name = cndtn[0]
        c_out[cndtn_name] = []
        for x in cndtn[1:]:
            if utils.is_float(x):
                c_out[cndtn_name].append(float(x))
            else:
                c_out[cndtn_name].append(x)

    conditions.conditions = c_out

    return conditions

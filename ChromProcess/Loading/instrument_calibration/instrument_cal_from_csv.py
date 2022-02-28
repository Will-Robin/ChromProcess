from pathlib import Path
from ChromProcess import Classes


def instrument_cal_from_csv(filename):
    """
    Create an InstrumentCalibration object from a formatted .csv file.

    Parameters
    ----------
    fname: str or pathlib Path
        name for file

    Returns
    -------
    calibration: Classes.InstrumentCalibration
    """

    if isinstance(filename, str):
        file = Path(filename)
    else:
        file = filename

    calibration = Classes.InstrumentCalibration()
    calibration.filename = file.name

    rdlin = lambda x: [e for e in x.strip("\n").split(",") if e != ""]
    info = []
    header = []
    upper_ind = 0
    lower_ind = 0
    with open(filename, "r") as f:
        readstate = False
        for line in f:
            ins = rdlin(line)
            if "Technique" in line:
                calibration.type = ins[1]
            if "Instrument" in line:
                calibration.instrument = ins[1]
            if "Date" in line:
                calibration.date = ins[1]
            if "Method_file" in line:
                calibration.method = ins[1]
            if "Derivatisation_details" in line:
                calibration.derivatisation = ";".join(ins[1:])
            if "calibration_model" in line:
                calibration.calibration_model = ins[1]
            if "internal_standard_position" in line:
                calibration.internal_standard_position = float(ins[1])
            if "Instrument" in line:
                calibration.instrument = ins[1]
            if "start_calibration" in line:
                readstate = True
            elif "compound_name" in line:
                header = line.strip("\n").split(",")
            elif "end_calibration" in line:
                break
            elif readstate:
                info.append(ins)

    calib_header = [
        "A",
        "A_variance",
        "B",
        "B_variance",
        "C",
        "C_variance",
        "AB_covariance",
        "AC_covariance",
        "BC_covariance",
    ]

    comp_ind = header.index("compound_name")
    lower_ind = header.index("lower")
    upper_ind = header.index("upper")

    for v in info:
        if "None" not in v:
            calibration.calibration_factors[v[comp_ind]] = {}
            for c in calib_header:
                idx = header.index(c)
                calibration.calibration_factors[v[comp_ind]][c] = float(v[idx])

    trans_info = [list(x) for x in zip(*info)]

    bound_tuples = [
        (l, u) for l, u in zip(trans_info[lower_ind], trans_info[upper_ind])
    ]
    calibration.boundaries = {
        k: [float(v[0]), float(v[1])]
        for k, v in zip(trans_info[comp_ind], bound_tuples)
        if v != ("None", "None")
    }

    return calibration

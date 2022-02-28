import os


def local_assignments_from_csv(file):
    """
    Read a comma-separated values file containing assignments of compounds to
    retention time boundaries.

    File format:
        compound_1_name,lower bound (float), upper bound (float),
        compound_2_name,lower bound (float), upper bound (float),
        etc.

    Parameters
    ----------
    file: str
        Path to file

    Returns
    -------
    modified_bounds: dict
        modified_bounds['compound_1_name'] = [lower: float, upper: float]
    """

    modified_bounds = {}

    if not os.path.exists(file):
        print(
            f"""ChromProcess.file_import.read_local_assignments:
        Local assignments file at {file} not found"""
        )
        print(
            """ChromProcess.file_import.read_local_assignments:
        Returning empty dict"""
        )
        return modified_bounds
    else:
        with open(file, "r") as f:
            for line in f:
                spl = line.strip("\n").split(",")
                spl = [x for x in spl if x != ""]
                if len(spl) == 0:
                    pass
                elif "set_IS_pos" in line:
                    modified_bounds[spl[0]] = [float(spl[1])]
                else:
                    modified_bounds[spl[0]] = [float(spl[1]), float(spl[2])]

    return modified_bounds

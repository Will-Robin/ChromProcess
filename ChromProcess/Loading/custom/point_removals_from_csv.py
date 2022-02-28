def read_point_removals_file(fname):
    """
    For reading a file containing datasets and point indices to remove in the
    data sets.

    Parameters
    ----------
    fname: str
        Path to the file

    Returns
    -------
    point_removals: dict
        dictionary of point removal indices. keys are experiment codes,
        list of ints are items.
    """

    point_removals = {}
    with open(fname, "r") as f:
        for c, line in enumerate(f):
            if c == 0:
                pass
            else:
                ins = line.strip("\n").split(",")
                point_removals[ins[0]] = [int(x) for x in ins[1:] if x != ""]

    return point_removals

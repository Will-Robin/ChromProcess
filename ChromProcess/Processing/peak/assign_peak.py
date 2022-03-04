def assign_retention_time(retention_time, boundaries):
    """
    Possible TODO: Re-implement as decision tree. The current method works, but
    relies upon Python 3 dictionary ordering and potentially some wasted
    iteration. However, performance is not currently an issue.

    Takes a peak (retention times of peak) and assigns peak name based on a
    dictionary of boundaries. The assignment priority is based on the
    interation order of the dict (Python 3 dict)

    Parameters
    ----------
    retention_time: float
        retention time of a peak
    boundaries: dict
        dictionary of boundaries for peak assignments

    Returns
    -------
    peak_name: str
        name of the peak. If the peak has no assignment, the peak retention time
        is returned
    """

    # Default assignment to rounded retention time.
    peak_name = str(round(retention_time, 3))
    for b in boundaries:
        if boundaries[b][0] < retention_time < boundaries[b][1]:
            peak_name = b
            break

    return peak_name

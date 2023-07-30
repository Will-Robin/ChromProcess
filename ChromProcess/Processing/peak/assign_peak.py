from ChromProcess.Classes import Peak


def assign_retention_time(peak: Peak, boundaries: dict[str, list[float]]) -> str:
    """
    Takes a peak (retention times of peak) and assigns peak name based on a
    dictionary of boundaries. The assignment priority is based on the
    iteration order of the dict (Python 3 dict).

    Parameters
    ----------
    retention_time: Peak
        Peak to be assigned.
    boundaries: dict
        dictionary of boundaries for peak assignments

    Returns
    -------
    peak_name: str
        Name for the peak. If the peak has no assignment, the peak retention
        time is returned, rounded to 3 decimal places.
    """

    # Default assignment to rounded retention time.
    peak_name = str(round(peak.retention_time, 3))
    for b in boundaries:
        if boundaries[b][0] < peak.retention_time < boundaries[b][1]:
            peak_name = b
            break

    return peak_name

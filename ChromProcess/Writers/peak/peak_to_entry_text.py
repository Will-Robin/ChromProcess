def peak_to_entry_text(peak, chromatogram):
    """
    For writing the peak's row in a peak table.

    Parameters
    ----------
    peak: Classes.Peak
    chromatogram: Classes.Chromatogram

    Returns
    -------
    entry: str
        Row for the peak table output
    """

    st_ind = 0
    end_ind = 0
    peak_start = 0.0
    peak_end = 0.0
    if len(peak.indices) > 0:
        st_ind = peak.indices[0]
        end_ind = peak.indices[-1]
        peak_start = chromatogram.time[st_ind]
        peak_end = chromatogram.time[end_ind]

    rtn_time = peak.retention_time
    integral = peak.integral

    entry = f"{rtn_time},{integral},{peak_start},{peak_end}\n"

    return entry

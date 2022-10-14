def peak_to_entry_text(peak, chromatogram):
    """
    For writing the peak's row in a peak table.

    Parameters
    ----------
    peak: ChromProcess.Classes.Peak
    chromatogram: ChromProcess.Classes.Chromatogram

    Returns
    -------
    entry: str
        Row for the peak table output
    """

    # Default values
    st_ind = 0
    end_ind = 0
    peak_start = 0.0
    peak_end = 0.0
    peak_height = 0.0

    if len(peak.indices) > 0:
        st_ind = peak.indices[0]
        end_ind = peak.indices[-1]

        peak_start = chromatogram.time[st_ind]
        peak_end = chromatogram.time[end_ind]
        peak_height = peak.height

    rtn_time = peak.retention_time
    integral = peak.integral

    entry_list = [rtn_time, integral, peak_start, peak_end, peak_height, "\n"]

    entry = ",".join([str(x) for x in entry_list])

    return entry

from ChromProcess import Classes


def analysis_from_csv_legacy(fname):
    """
    Create and AnalysisInformation object from a formatted csv file.
    Written for an older version of the file format.

    Parameters
    ----------
    fname: str or pathlib Path

    Returns
    -------
    analysis: Classes.AnalysisInformation
    """

    analysis = Classes.AnalysisInformation()

    rdlin = lambda x: [e for e in x.strip("\n").split(",") if e != ""]

    with open(fname, "r") as file:
        lines = file.readlines()

    for line in lines:
        ins = rdlin(line)
        if "Dataset" in line:
            analysis.experiment_code = ins[1]

        if "Method" in line and not "Instrument" in line:
            analysis.analysis_type = ins[1]

        if "regions" in line:
            reg = [float(x) for x in ins[1:]]
            analysis.regions = [reg[x : x + 2] for x in range(0, len(reg), 2)]

        if "internal_reference_region" in line:
            reg = [float(x) for x in ins[1:]]
            analysis.internal_standard_region = reg

        if "extract_mass_spectra" in line:
            use_ms = ins[1].lower()
            if use_ms == "true":
                analysis.use_MS = True
            if use_ms == "false":
                analysis.use_MS = False

        if "mass_spectra_filter" in line:
            analysis.MS_cutoff = float(ins[1])

        if "peak_pick_threshold" in line:
            analysis.peak_pick_threshold = float(ins[1])

        if "dilution_factor," in line:
            analysis.dilution_factor = float(ins[1])

        if "dilution_factor_error" in line:
            analysis.dilution_factor_error = float(ins[1])

        if "internal_ref_concentration," in line:
            analysis.internal_standard_concentration = float(ins[1])

        if "internal_ref_concentration_error" in line:
            analysis.internal_standard_concentration_error = float(ins[1])

        if "instrument" in line and not "method" in line:
            analysis.instrument = ins[1]

        if "instrument_method" in line:
            analysis.instrument_method = ins[1]

        if "derivatisation_method" in line:
            analysis.derivatisation_method = ins[1]

        if "calibration_model" in line:
            analysis.calibration_model = ins[1]

        if "calibration_file" in line:
            analysis.calibration_file = ins[1]

    return analysis

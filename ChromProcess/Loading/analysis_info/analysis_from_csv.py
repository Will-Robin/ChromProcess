from ChromProcess.Classes import AnalysisInformation

def rdlin(x):
    return [e for e in x.strip("\n").split(",") if e != ""]

def analysis_from_csv(fname: str) -> AnalysisInformation:
    """
    Create and AnalysisInformation object from a formatted csv file.

    Parameters
    ----------
    fname: str

    Returns
    -------
    analysis: AnalysisInformation
    """

    analysis = AnalysisInformation()

    with open(fname, "r") as file:
        lines = file.readlines()

    for line in lines:
        ins = rdlin(line)
        if "Dataset" in line:
            analysis.experiment_code = ins[1]

        if "Method" in line and "Instrument" not in line:
            analysis.analysis_type = ins[1]

        if "Regions" in line:
            reg = [float(x) for x in ins[1:]]
            analysis.regions = [reg[x : x + 2] for x in range(0, len(reg), 2)]

        if "Internal_standard_region" in line:
            reg = [float(x) for x in ins[1:]]
            analysis.internal_standard_region = reg

        if "Extract_mass_spectra" in line:
            use_ms = ins[1].lower()
            if use_ms == "true":
                analysis.use_MS = True
            if use_ms == "false":
                analysis.use_MS = False

        if "Mass_spectra_filter" in line:
            analysis.MS_cutoff = float(ins[1])

        if "Peak_pick_threshold" in line:
            analysis.peak_pick_threshold = float(ins[1])

        if "Dilution_factor," in line:
            analysis.dilution_factor = float(ins[1])

        if "Dilution_factor_error" in line:
            analysis.dilution_factor_error = float(ins[1])

        if "Internal_standard_concentration," in line:
            analysis.internal_standard_concentration = float(ins[1])

        if "Internal_standard_concentration_error" in line:
            analysis.internal_standard_concentration_error = float(ins[1])

        if "Instrument" in line and "method" not in line:
            analysis.instrument = ins[1]

        if "Instrument_method" in line:
            analysis.instrument_method = ins[1]

        if "Derivatisation_method" in line:
            analysis.derivatisation_method = ins[1]

        if "Calibration_model" in line:
            analysis.calibration_model = ins[1]

        if "Calibration_file" in line:
            analysis.calibration_file = ins[1]

    return analysis

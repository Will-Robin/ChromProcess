import tomli
from ChromProcess import Classes


def analysis_from_toml(fname):
    """
    Create and AnalysisInformation object from a formatted csv file.

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
        text = file.read()

    analysis_dict = tomli.loads(text)

    analysis.experiment_code = analysis_dict["Dataset"]

    method_info = analysis_dict["Method_information"]
    analysis.analysis_type = method_info["Method"]

    analysis.instrument = method_info["Instrument"]
    analysis.instrument_method = method_info["Instrument_method"]
    analysis.derivatisation_method = method_info["Derivatisation_method"]
    analysis.calibration_model = method_info["Calibration_model"]
    analysis.calibration_file = method_info["Calibration_file"]

    regions = analysis_dict["Chromatogram_Regions"]
    analysis.regions = [regions[x] for x in regions]

    is_reg = analysis_dict["Internal_standard_region"]
    analysis.internal_standard_region = is_reg["internal_standard_1"]

    ms_opts = analysis_dict["Mass_spectra_options"]
    use_ms = ms_opts["Extract_mass_spectra"].lower()

    if use_ms == "true":
        analysis.use_MS = True
    if use_ms == "false":
        analysis.use_MS = False

    analysis.MS_cutoff = ms_opts["Mass_spectra_filter"]

    peak_picking_options = analysis_dict["Peak_picking_options"]
    analysis.peak_pick_threshold = peak_picking_options["Peak_pick_threshold"]

    sample_info = analysis_dict["Sample_information"]
    analysis.dilution_factor = sample_info["Dilution_factor"]
    analysis.dilution_factor_error = sample_info["Dilution_factor_error"]
    analysis.internal_standard_concentration = sample_info["Internal_standard_concentration"]
    analysis.internal_standard_concentration_error = sample_info["Internal_standard_concentration_error"]

    return analysis

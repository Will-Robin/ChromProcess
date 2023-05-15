from ChromProcess.Classes import ExperimentConditions
from ChromProcess.Classes import AnalysisInformation


def write_conditions_header(
    name: str, conditions: ExperimentConditions, information: AnalysisInformation
) -> str:
    """
    Write the text for a formatted file header detailing experiment
    information.

    Parameters
    ----------
    name: str
    conditions: ExperimentConditions
    information: AnalysisInformation

    Returns
    -------
    header_text: str
    """

    header_text = ""

    header_text += f"Dataset,{name}\n"
    header_text += "start_conditions\n"
    for c in conditions.conditions:
        header_text += f"{c},"

        if type(conditions.conditions[c]) == list:
            for x in conditions.conditions[c]:
                header_text += f"{x},"
        else:
            header_text += f"{conditions.conditions[c]},"

        header_text += "\n"

    header_text += "end_conditions\n"

    header_text += "start_analysis_details\n"
    header_text += f"Instrument,{information.instrument}\n"
    header_text += f"Chromatography_type,{information.analysis_type}\n"
    header_text += f"Chromatography_method,{information.instrument_method}\n"
    deriv_method = information.derivatisation_method.replace(",", ";")
    deriv_list = deriv_method.split("\n")
    write_deriv = " ".join(deriv_list)
    header_text += f"Derivatisation_method,{write_deriv}\n"
    header_text += f"Calibrations_file,{information.calibration_file}\n"
    header_text += f"Calibration_model,{information.calibration_model}\n"
    header_text += f"end_analysis_details\n"

    return header_text

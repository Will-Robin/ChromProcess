
def write_conditions_header(name, conditions, information):
    '''
    Parameters
    ----------
    name: str
    conditions: Classes.ExperimentConditions
    information: ChromProcess AnalysisInformation object

    Returns
    -------
    header_text: str
    '''

    header_text = ''
    # writing experiment conditions to file
    header_text += f"Dataset,{name}\n"
    header_text += "start_conditions\n"
    for c in conditions:
        header_text += f"{c},"

        for x in conditions[c]:
            header_text += f"{x},"

        header_text += "\n"
    header_text += "end_conditions\n"

    # writing analysis details
    header_text += "start_analysis_details\n"
    header_text += f'Instrument, {information.instrument}\n'
    header_text += f"Chromatography_method,{information.analysis_type},{information.instrument_method}\n"
    header_text += f"Derivatisation_method,{information.derivatisation_method}\n"
    header_text += f"Calibrations_file,{information.calibration_file}\n"
    header_text += f'Calibration_model,{information.calibration_model}\n'
    header_text += f"end_analysis_details\n"

    return header_text


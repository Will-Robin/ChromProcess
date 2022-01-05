
def write_conditions_text(data_report):
    '''
    Parameters
    ----------
    data_report: Classes.DataReport
    '''

    header_string = ''

    # writing experiment conditions to file
    header_string += f"Dataset,{data_report.experiment_code}\n"
    header_string += "start_conditions\n"
    for c in data_report.conditions:
        header_string += f"{c},"
        for x in data_report.conditions[c]:
            header_string += f"{x},"

        header_string += "\n"
    header_string += "end_conditions\n"

    # writing analysis details
    header_string += "start_analysis_details\n"
    for ad in data_report.analysis_details:
        header_string += f"{ad},"
        if isinstance(data_report.analysis_details[ad], str):
            header_string += f"{data_report.analysis_details[ad]},"
        else:
            for x in data_report.analysis_details[ad]:
                header_string += f"{x},"
        header_string += "\n"

    header_string += "end_analysis_details\n"

    return header_string

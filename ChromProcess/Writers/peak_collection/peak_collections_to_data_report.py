from pathlib import Path
from ChromProcess.Classes import PeakCollection
from ChromProcess.Classes import AnalysisInformation

from ChromProcess.Utils.utils import utils
from ChromProcess.Writers.general import write_header


def series_traces_as_dict(peak_collections: list[PeakCollection]):
    """
    Create dictionaries of peak series values derived using the
    self.clusters.

    Parameters
    ----------

    Returns
    -------
    conc_dict: dict
    err_dict: dict
    integral_dict: dict
    """

    conc_dict = dict()
    err_dict = dict()
    integral_dict = dict()
    height_dict = dict()

    if len(self.cluster_assignments) == 0:
        cluster_names = ["" for _ in self.clusters]
    else:
        cluster_names = [n.split(" ")[0] for n in self.cluster_assignments]

    for c1, pc in enumerate(self.peak_collections):

        for c2, clust in enumerate(self.clusters):
            if len(clust) == 0:
                continue

            name = cluster_names[c2]
            average_position = round(sum(clust) / len(clust), 3)

            for pk in pc.peaks:
                if pk.retention_time in clust:
                    if pk.integral:
                        token = name + f" [{average_position}]"
                        if token not in integral_dict:
                            integral_dict[token] = [0.0 for _ in self.series_values]
                        integral_dict[token][c1] += pk.integral

                    if pk.height:
                        token = name + f" [{average_position}]"
                        if token not in height_dict:
                            height_dict[token] = [0.0 for _ in self.series_values]
                        height_dict[token][c1] += pk.height

                    if pk.concentration:
                        token = name + "/ M"
                        if token not in conc_dict:
                            conc_dict[token] = [0.0 for _ in self.series_values]
                        conc_dict[token][c1] += pk.concentration

                    if pk.conc_error:
                        token = name + "/ M"
                        if token not in err_dict:
                            err_dict[token] = [0.0 for _ in self.series_values]
                        err_dict[token][c1] += pk.conc_error

    return conc_dict, err_dict, integral_dict, height_dict


def peak_collection_series_conc_report_text(
    peak_collection: list[PeakCollection], information: AnalysisInformation
) -> str:
    """
    Write a peak collection series as formatted data report text.

    Parameters
    ----------
    peak_collection_series: list[PeakCollection]
    information: AnalysisInformation

    Returns
    -------
    conc_report_text: str
    """

    # create output dictionaries
    conc_dict, err_dict, _, _ = series_traces_as_dict()

    # create spreadsheet-like output
    conc_header, conc_grid = utils.peak_dict_to_spreadsheet(
        conc_dict,
        peak_collection_series.series_values,
        peak_collection_series.series_unit,
    )

    peak_err_header, err_grid = utils.peak_dict_to_spreadsheet(
        err_dict,
        peak_collection_series.series_values,
        peak_collection_series.series_unit,
    )

    header_text = write_header.write_conditions_header(
        peak_collection_series.name, peak_collection_series.conditions, information
    )

    conc_report_text = header_text

    conc_report_text += "start_data\n"

    for x in conc_header:
        conc_report_text += f"{x},"

    conc_report_text += "\n"

    for x in range(0, len(conc_grid)):
        for y in range(0, len(conc_grid[x])):
            val = conc_grid[x][y]
            conc_report_text += f"{val},"
        conc_report_text += "\n"

    conc_report_text += "end_data\n"

    conc_report_text += "start_errors\n"

    for x in peak_err_header:
        conc_report_text += f"{x},"

    conc_report_text += "\n"

    for x in range(0, len(err_grid)):
        for y in range(0, len(err_grid[x])):
            val = err_grid[x][y]
            conc_report_text += f"{val},"
        conc_report_text += "\n"

    conc_report_text += "end_errors\n"

    return conc_report_text


def peak_collection_series_integral_report_text(peak_collection_series, information):
    """
    Write a peak collection series as formatted data report text.

    Parameters
    ----------
    peak_collection_series: Classes.PeakCollectionSeries
    information: ChromProcess Analysis_Information object

    Returns
    -------
    integral_report_text: str
    """

    _, _, integral_dict, _ = peak_collection_series.series_traces_as_dict()

    header_text = write_header.write_conditions_header(
        peak_collection_series.name, peak_collection_series.conditions, information
    )

    peak_integral_header, integ_grid = utils.peak_dict_to_spreadsheet(
        integral_dict,
        peak_collection_series.series_values,
        peak_collection_series.series_unit,
    )

    integral_report_text = header_text

    integral_report_text += "start_data\n"

    for x in peak_integral_header:
        integral_report_text += f"{x},"
    integral_report_text += "\n"

    for x in range(0, len(integ_grid)):
        for y in range(0, len(integ_grid[x])):
            val = integ_grid[x][y]
            integral_report_text += f"{val},"
        integral_report_text += "\n"

    integral_report_text += "end_data\n"

    return integral_report_text


def peak_collection_series_to_data_report(
    peak_collection_series, filename, information
):
    """
    Write a peak collection series as a formatted data report file.

    Parameters
    ----------
    peak_collection_series: ChromProcess.Classes.PeakCollectionSeries
    filename: name for file including path
    information: ChromProces.Classes.AnalysisInformation

    Returns
    -------
    None
    """

    if isinstance(filename, str):
        filename = filename
    elif isinstance(filename, Path):
        filename = str(filename)

    analysis_type = information.analysis_type

    conc_fname = f"{filename}_{analysis_type}_concentration_report.csv"
    integral_fname = f"{filename}_{analysis_type}_integral_report.csv"
    height_fname = f"{filename}_{analysis_type}_peak_height_report.csv"

    # create output dictionaries
    (
        conc_dict,
        err_dict,
        integral_dict,
        height_dict,
    ) = peak_collection_series.series_traces_as_dict()

    # create spreadsheet-like output
    # Peak concentrations
    conc_header, conc_grid = utils.peak_dict_to_spreadsheet(
        conc_dict,
        peak_collection_series.series_values,
        peak_collection_series.series_unit,
    )

    # Peak integrals
    peak_integral_header, integ_grid = utils.peak_dict_to_spreadsheet(
        integral_dict,
        peak_collection_series.series_values,
        peak_collection_series.series_unit,
    )

    # Peak heights
    peak_height_header, height_grid = utils.peak_dict_to_spreadsheet(
        height_dict,
        peak_collection_series.series_values,
        peak_collection_series.series_unit,
    )

    # Peak errors
    peak_err_header, err_grid = utils.peak_dict_to_spreadsheet(
        err_dict,
        peak_collection_series.series_values,
        peak_collection_series.series_unit,
    )

    # Create a header containing experimental conditions
    header_text = write_header.write_conditions_header(
        peak_collection_series.name, peak_collection_series.conditions, information
    )

    # Write concentration report to file
    conc_rep_text = peak_collection_series_conc_report_text(
        peak_collection_series, information
    )
    with open(conc_fname, "w") as outfile:
        outfile.write(conc_rep_text)

    # Write integral report to file
    i_data_rep_text = peak_collection_series_integral_report_text(
        peak_collection_series, information
    )
    with open(integral_fname, "w") as outfile:
        outfile.write(i_data_rep_text)

        outfile.write(header_text)

        outfile.write("start_data\n")

        [outfile.write("{},".format(x)) for x in peak_integral_header]

        outfile.write("\n")
        for x in range(0, len(integ_grid)):
            for y in range(0, len(integ_grid[x])):
                val = integ_grid[x][y]
                outfile.write(f"{val},")
            outfile.write("\n")

        outfile.write("end_data\n")

    # Write height report to file
    with open(height_fname, "w") as outfile:

        outfile.write(header_text)

        outfile.write("start_data\n")

        [outfile.write("{},".format(x)) for x in peak_height_header]

        outfile.write("\n")
        for x in range(0, len(height_grid)):
            for y in range(0, len(height_grid[x])):
                val = height_grid[x][y]
                outfile.write(f"{val},")
            outfile.write("\n")

        outfile.write("end_data\n")

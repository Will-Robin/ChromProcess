"""
Creates a template directory for an experiment

Consider using cookiecutter (https://github.com/cookiecutter/cookiecutter).
"""
import os
import sys
from pathlib import Path

# The home directory is default.
default_parent_directory = Path.home()

exp_code = "experiment_code"
default_folder_name = exp_code

# defaults for the conditions file
default_series_values = ["1", "2", "3", "4", "5", "6", "7", "8", "9"]
default_series_unit = "Sample_number/ n"
example_condition_name = "Condition 1"
example_condition_value = "Condition_1_value"

# defaults for the analysis file
default_analysis_method = "GC"
default_regions = ["0", "1", "2", "3"]
default_internal_standard_region = ["6.7", "6.8"]
mass_spec_default = "FALSE"
default_peak_pick_threshold = "0.1"
default_dilution_factor = "1"
default_dilution_factor_error = "0.0447"
default_internal_standard_concentration = "0.0008"
default_internal_standard_concentration_error = "9.89E-06"

# basis for conditions file
conditions = {
    "Dataset": [exp_code],
    "start_experiment_information": [],
    "Series_values": default_series_values,
    "Series_unit": [default_series_unit],
    "end_experiment_information": [],
    "start_conditions": [],
    "example_condition_name": [example_condition_value],
    "end_conditions": [],
}

analysis_details = {
    "Dataset": [exp_code],
    "Method": [default_analysis_method],
    "Regions": default_regions,
    "Internal_standard_region": default_internal_standard_region,
    "Extract_mass_spectra": [mass_spec_default],
    "Peak_pick_threshold": [default_peak_pick_threshold],
    "Dilution_factor": [default_dilution_factor],
    "Dilution_factor_error": [default_dilution_factor_error],
    "Internal_standard_concentration": [default_internal_standard_concentration],
    "Internal_standard_concentration_error": [
        default_internal_standard_concentration_error
    ],
}

local_assignments = {}


def write_dict_to_file(fname, this_dict):
    """
    fname: str
        Name for the file.
    this_dict: dict
        Dictionary to be written.
    """

    with open(fname, "w") as f:
        for k in this_dict:
            if len(this_dict[k]) == 0:
                info = ""
            else:
                info = ",".join(this_dict[k])

            f.write(f"{k},{info}\n")


def build_project(
    project_name, parent_directory, conditions, analysis_details, local_assignments
):

    """
    project_name: str
        A name for the project.
    """

    # create the project folder and write files into it
    os.makedirs(f"{parent_directory}/{project_name}", exist_ok=True)

    project_root = f"{parent_directory}/{project_name}"

    write_dict_to_file(f"{project_root}/{project_name}_conditions.csv", conditions)

    write_dict_to_file(
        f"{project_root}/{project_name}_analysis_details.csv", analysis_details
    )

    write_dict_to_file(
        f"{project_root}/{project_name}_local_assignments.csv", local_assignments
    )


if __name__ == "__main__":

    folder_name = default_folder_name
    parent_directory = default_parent_directory

    if len(sys.argv) >= 2:
        folder_name = sys.argv[1]

    if len(sys.argv) >= 3:
        parent_directory = sys.argv[2]

    print(f"Creating {folder_name} in {parent_directory}.")

    build_project(
        folder_name,
        parent_directory,
        conditions,
        analysis_details,
        local_assignments,
    )

from ChromProcess.Loading import conditions_from_toml
from ChromProcess.Loading import conditions_from_csv

from ChromProcess.Loading import analysis_from_csv
from ChromProcess.Loading import analysis_from_toml

# Analysis
toml_file = "Example/Analysis/example_analysis_details.toml"
csv_file = "Example/Analysis/example_analysis_details.csv"

analysis_toml = analysis_from_toml(toml_file)
analysis_csv = analysis_from_csv(csv_file)

print("TOML")

print(analysis_toml.experiment_code)
print(analysis_toml.analysis_type)
print(analysis_toml.derivatisation_method)
print(analysis_toml.instrument)
print(analysis_toml.instrument_method)
print(analysis_toml.calibration_model)
print(analysis_toml.calibration_file)
print(analysis_toml.regions)
print(analysis_toml.internal_standard_region)
print(analysis_toml.use_MS)
print(analysis_toml.MS_cutoff)
print(analysis_toml.peak_pick_threshold )
print(analysis_toml.dilution_factor )
print(analysis_toml.dilution_factor_error )
print(analysis_toml.internal_standard_concentration )
print(analysis_toml.internal_standard_concentration_error )

print("CSV")
print(analysis_csv.experiment_code)
print(analysis_csv.analysis_type)
print(analysis_csv.derivatisation_method)
print(analysis_csv.instrument)
print(analysis_csv.instrument_method)
print(analysis_csv.calibration_model)
print(analysis_csv.calibration_file)
print(analysis_csv.regions)
print(analysis_csv.internal_standard_region)
print(analysis_csv.use_MS)
print(analysis_csv.MS_cutoff )
print(analysis_csv.peak_pick_threshold )
print(analysis_csv.dilution_factor )
print(analysis_csv.dilution_factor_error )
print(analysis_csv.internal_standard_concentration )
print(analysis_csv.internal_standard_concentration_error )

# Conditions
toml_file = "Example/ExperimentalData/example_conditions.toml"
csv_file = "Example/ExperimentalData/example_conditions.csv"

conditions_toml = conditions_from_toml(toml_file)
conditions_csv = conditions_from_csv(csv_file)

print("TOML")
print(conditions_toml.conditions)
print(conditions_toml.experiment_code)
print(conditions_toml.series_values)
print(conditions_toml.series_unit)

print("CSV")
print(conditions_csv.conditions)
print(conditions_csv.experiment_code)
print(conditions_csv.series_values)
print(conditions_csv.series_unit)

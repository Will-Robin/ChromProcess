import os
from ChromProcess import Classes
from ChromProcess.Loading import peak_collection_from_csv
from ChromProcess.Loading import analysis_from_csv
from ChromProcess.Loading import conditions_from_csv

experiment_number = 'FRN141'
experiment_folder = r"C:\Users\thijs\Documents\PhD\Data\FRN141"
peak_collection_directory = f'{experiment_folder}\PeakCollections'
conditions_file = f'{experiment_folder}\{experiment_number}_conditions.csv'
analysis_file = f'{experiment_folder}\{experiment_number}_analysis_details.csv'
data_report_directory = f'{experiment_folder}\DataReports'
os.makedirs(data_report_directory, exist_ok=True)

conditions = conditions_from_csv(conditions_file)
analysis = analysis_from_csv(analysis_file)

peak_tables = []
for file in os.listdir(peak_collection_directory):
    if file.endswith('.csv'):
        peak_tables.append(peak_collection_from_csv(f'{peak_collection_directory}/{file}'))

# Create series of peak collections
series = Classes.PeakCollectionSeries(
                                    peak_tables, 
                                    name = f'{experiment_number}',
                                    conditions = conditions.conditions
                                    )
IS_pos = 7.43
series.align_peaks_to_IS(IS_pos)
series.reference_integrals_to_IS()
 # 5% of internal standard integral if integrals are normalised to IS
#series.remove_peaks_below_threshold(peak_removal_limit)
peak_agglomeration_boundary = 0.025 # distance cutoff 
series.get_peak_clusters(bound = peak_agglomeration_boundary)
cluster_removal_limit = 0.008


series.write_data_reports(f'{data_report_directory}/{series.name}', analysis,cluster_removal_limit = cluster_removal_limit) # create arrays for output

import os
import numpy as np
from pathlib import Path
from ChromProcess import file_import
from ChromProcess import file_output
from ChromProcess import peak_operations as peak_ops
from ChromProcess import processing_functions as p_f
from ChromProcess import chromatogram_operations as chrom_ops

'''
A program for loading sets of chromatograms from the source_folder and
extracting basic information from them (chromatogram traces, peak integrals,
mass spectra), which are deposited in the store_folder.

Minimal processing is performed on the chromatograms.
'''
# Parameters used in the analysis
use_MS = True # Whether to extract mass spectra or not
analysis_type = 'GCMS' # Which kind of data are to be analysed
inten_thres_MS = 500 # Threshold for removing background mass spectra
peak_pick_thres = 0.1 # Threshold for peak picking as a proportion of the max
                      # signal in a region.

# Point the program to where the data are
source_folder =  Path(r'C:\Users\willi\Documents\Data\GCMS\FRN\FRN107\A')
# State directory in which to store results
store_folder = Path(r'C:\Users\willi\Documents\Data\GCMS\FRN\FRN107\A')
os.makedirs(store_folder, exist_ok = True)

# Read in the data files
if analysis_type == 'GCMS':
    # load .cdf files from GCMS analysis
    chroms, cond_file = file_import.load_cdf_from_directory(source_folder,
                                                            ms = use_MS)
elif analysis_type == 'HPLC':
    # load .txt files exported from Shimadzu LabSolutions for HPLC analysis
    chroms, cond_file = file_import.directoryLoadShimadzuASCII(source_folder)
else:
    print('analysis_type provided is {}.')
    print('Please choose from GCMS or HPLC. Quitting.'.format(analysis_type))
    quit()

# Get experiment conditions information
conditions = file_import.ImportExperimentConditionsFile(cond_file)

# Pre-process the chromatograms
for c in chroms:
    if analysis_type == 'GCMS':
        # remove low intensity ion chromatograms
        # reconstitue total ion chromatogram
        p_f.MS_intensity_threshold_chromatogram(c, threshold = inten_thres_MS)

# Find integral information
for c in chroms:
    if analysis_type == "GCMS":
        # Get interal reference integrals
        chrom_ops.internalRefIntegral(c, conditions.internal_ref_region)

    # Get peaks in regions of the chromatogram
    for r in conditions.regions:
        chrom_ops.pickPeaksRegion(c, r, threshold = peak_pick_thres)

# Output chromatograms
os.makedirs(store_folder/'Chromatograms', exist_ok = True)
dest_dir = store_folder/'Chromatograms'
for c in chroms:
    file_output.chromatogram_to_csv_GCMS(c, filename = dest_dir/c.filename)

# Output peak table
os.makedirs(store_folder/'PeakTables', exist_ok = True)
dest_dir = store_folder/'PeakTables'
for c,v in zip(chroms, conditions.x_series):
    file_output.write_peak_table(c, filename = dest_dir/c.filename, value = v,
                                 series_unit = conditions.x_name)
if analysis_type == 'GCMS':
    # Output peak mass spectra
    os.makedirs(store_folder/'PeakMassSpectra', exist_ok = True)
    dest_dir = store_folder/'PeakMassSpectra'
    for c in chroms:
        for p in c.peaks:
            peak_ops.peakIonChromatogram(c.peaks[p],c)

        file_output.write_peak_mass_spectra(c, filename = dest_dir/c.filename)

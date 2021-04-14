import sys
sys.path.append(r'/Users/williamrobinson/Documents/Nijmegen/packages')
from ChromProcess import Classes, os, plt, series_builder, plotting, file_output, info_params
from ChromProcess import processing_functions
from pathlib import Path
import numpy as np
from NorthNet import Global_formatting

min_font = Global_formatting.min_font # should be around 1 mm on paper
col_width = Global_formatting.col_width
two_col_width = Global_formatting.two_col_width
full_page = Global_formatting.full_page
modest_display_item = Global_formatting.modest_display_item
composite_figure = Global_formatting.composite_figure
lines = Global_formatting.lines
ticklength = Global_formatting.ticklength
tickpad = Global_formatting.tickpad

directory = Path(r'/Users/williamrobinson/Documents/Nijmegen/Dynamic_environment_project')
save_data_folder = directory/'DataAnalysis'
plot_dir = Path(directory/'Paper_plots/SugarChromOverlay')

regions = {'C4':[7.9,8.4],'C5':[10.25,11.25],'C6':[13.5,15.5],'C7':[18,21]}

os.chdir(r"/Users/williamrobinson/Documents/Nijmegen/Data/2020_calibrations/All_sugar_overlay") # changes the directory to where files are stored

calib_file = r"/Users/williamrobinson/documents/nijmegen/packages/info_files/2020_03_16_Combined_Species_Info.csv"
IS_pos = 6.74 # 6.956/6.925/ 6.74
information = Classes.Calibration_File(calib_file, type = "GCMS")

filelist = os.listdir() # get list of files in directory
filelist.sort() # sort the files

'''Writes the chromatograms to csv files and creates Chromatogram objects from
.cdf files'''
chroms = []
names = []
specs = []
for f in filelist:
    if f.endswith(".cdf"):
        print("loading", f)
        chroms.append(Classes.Chromatogram(f, mass_spec = True))
        names.append(f.split("_")[0])
        specs.append(f[3:-6])
    if "conditions" in f:
        cond_file = f

series = Classes.Chromatogram_Series(chroms, cond_file)
series.regions = [regions[r] for r in regions]
counter = 0
for b in series.regions:

    for c in chroms: # plotting chromatograms for inspection.
        name = c.filename.split('_')[1].lower()

        if name == 'dha':
            name ='dihydroxyacetone'

        if name not in info_params.canonical_SMILES:
            print(name)
            continue

        smiles = info_params.canonical_SMILES[name]
        processing_functions.MS_intensity_threshold_chromatogram(c, threshold = 500)

        region_inds = np.where((c.time > b[0]) & (c.time < b[1]))[0]

        loc_time = c.time[region_inds]
        signal = c.signal[region_inds]

        if np.amax(signal) < 5.5e4:
            pass
        else:
            fig, ax = plt.subplots(figsize = (two_col_width/3,two_col_width/3))
            peak_ind = np.where(signal == np.amax(signal))[0]
            #name = processing_functions.name_peak(loc_time[peak_ind[0]],information.boundaries)
            ax.plot(loc_time, signal/1000000, c = info_params.colour_assignments[smiles])
            #ax.annotate(name, xy = (loc_time[peak_ind],signal[peak_ind]))

            ax.set_xlabel('Retention_time/ min.', fontsize = min_font*3)
            ax.set_ylabel('Total ion count/ 10$^6$', fontsize = min_font*3)
            ax.tick_params(axis='both', which='major', labelsize = min_font*2,
                           length = ticklength*2, pad = tickpad*2)

            ax.set_position([0.2,0.15,0.75,0.75])
            ax.set_title(name, fontsize = min_font*4)

            counter += 1
            plt.savefig(plot_dir/'{}Chromatogram.png'.format(name), dpi = 300)
            plt.close()
print(counter)

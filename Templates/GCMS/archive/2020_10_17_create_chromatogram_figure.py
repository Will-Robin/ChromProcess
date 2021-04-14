import sys
sys.path.append(r'C:\Users\willi\Documents\Packages')
from ChromProcess import Classes, os, plt, series_builder, plotting, file_output, file_import
from ChromProcess import mass_spectra, calibration_functions, info_params
from ChromProcess import processing_functions
import pickle
import numpy as np
import matplotlib.patches as patches
from matplotlib.patches import ConnectionPatch
import matplotlib.transforms as mtransforms
sys.path.append(r'C:\Users\willi\Documents\Packages\Scripts\Data_processing\Clustering\Workflow')
import Global_formatting
from pathlib import Path
from scipy import interpolate


min_font = Global_formatting.min_font # should be around 1 mm on paper
col_width = Global_formatting.col_width
two_col_width = Global_formatting.two_col_width
full_page = Global_formatting.full_page
modest_display_item = Global_formatting.modest_display_item
composite_figure = Global_formatting.composite_figure
lines = Global_formatting.lines/6
ticklength = Global_formatting.ticklength*2
tickpad = Global_formatting.tickpad*2

min_font = min_font
figure_width = 4.7/2.54
figure_height = 5.67/2.54
axis_width = 0.5
axis_height = 0.5
lines = lines*2

chrom_limits = (9,18)
zoom_lu = (10.1,11)

colour_map = plt.cm.cividis
C_chain_regions = {"tetradecane" : [6.5,6.9],
                   "C3": [5.2,6.6],
                   "C4": [7.5,8.8],
                   "C5": [9.3,12],
                   "C6": [12,18]}

C_chain_colors = { "tetradecane" : "#000000",
                   "C3": "#f09c08",
                   "C4": "#2738e7",
                   "C5": "#cb340b",
                   "C6": "#30bd37",
                   "C7": "#592387"}

directory = Path(r"C:\Users\willi\Documents\Data\GCMS\FRN\FRN089\B")

chroms, cond_file = file_import.load_cdf_from_directory(directory, ms = False) # get experimental data and conditions
series = Classes.Chromatogram_Series(chroms, cond_file) # Convert chromatogram list and conditions into a Series object

axw = 0.1
axh= 0.1
fig,ax = plt.subplots(ncols = len(chroms),figsize = (figure_width,figure_height))
positions = np.flip(np.linspace(0.1,0.9,len(chroms)))
for c,ch in enumerate(chroms):
    for reg in C_chain_regions:
        inds = np.where((ch.time>C_chain_regions[reg][0])&(ch.time<C_chain_regions[reg][1]))[0]
        ax[c].plot(ch.time[inds],ch.signal[inds], linewidth = lines, c = C_chain_colors[reg])

    L = positions[c]-axw/2
    B = positions[c]-axh/2
    W = axw
    H = axh
    ax[c].set_position([L,B,W,H])
    ax[c].set_xticklabels([])
    ax[c].set_yticklabels([])
    ax[c].patch.set_alpha(0.0)
    ax[c].spines['right'].set_visible(False)
    ax[c].spines['top'].set_visible(False)

plt.savefig(r'C:\Users\willi\Documents\Dynamic_environment_project\Paper_plots\ExampleChromStack.png', dpi = 800, transparent = True)

plt.close()

t_max = []
t_min = []
for c in chroms:
    #processing_functions.MS_intensity_threshold_chromatogram(c, threshold = 500)
    t_min.append(c.time[0])
    t_max.append(c.time[-1])

maxt = min(t_max)
mint = max(t_min)
heat_stack = np.zeros((len(chroms),13236))
for cn,c in enumerate(chroms):
    f = interpolate.interp1d(c.time,c.signal, kind = "linear")
    xnew = np.linspace(mint,maxt, num = 13236)
    ynew = f(xnew)
    heat_stack[cn] = ynew

fchrommap_inds = np.where((xnew>chrom_limits[0])&(xnew<chrom_limits[1]))[0]

X = xnew[fchrommap_inds]
Y = np.array([x/1000 for x in series.x_series])
zi = heat_stack[:,fchrommap_inds]
xi, yi = np.meshgrid(X, Y)

file = 'FRN089_036.cdf'
home = os.getcwd()
os.chdir(directory)
chrom = Classes.Chromatogram(file, mass_spec = True)
os.chdir(home)
processing_functions.MS_intensity_threshold_chromatogram(chrom, threshold = 500)

fchrom_inds = np.where((chrom.time>chrom_limits[0])&(chrom.time<chrom_limits[1]))[0]
i_n = np.where((chrom.time >  C_chain_regions['tetradecane'][0])&(chrom.time < C_chain_regions['tetradecane'][1]))[0]
IS_height = np.amax(chrom.signal[i_n])

from matplotlib.collections import PolyCollection
fig  = plt.figure(1, figsize = (figure_width,figure_height))

ax = fig.add_subplot(2, 2, 2)
ax.plot(chrom.time,chrom.signal/IS_height, linewidth = lines, c = '#000000', alpha=0.5)
for reg in C_chain_regions:
    inds = np.where((chrom.time>C_chain_regions[reg][0])&(chrom.time<C_chain_regions[reg][1]))[0]
    ax.plot(chrom.time[inds],chrom.signal[inds]/IS_height, linewidth = lines, c = C_chain_colors[reg])
ax.set_ylabel('Signal/ IS height',fontsize = min_font*2)
ax.set_xticklabels([])
ax.set_xlim(np.amin(chrom.time),chrom_limits[1])

ax = fig.add_subplot(2, 2, 3)
inds = np.where((xnew>zoom_lu[0])&(xnew<zoom_lu[1]))[0]
xi, yi = np.meshgrid(xnew[inds], Y)
zi = heat_stack[:,inds]
CS = ax.contourf(xi,yi,zi, 100, cmap=colour_map)
ax.set_xlabel('Retention time/ min.',fontsize = min_font*2)
ax.set_ylabel('time/ ks',fontsize = min_font*2)

cent_WH = [ [(0.56,0.75),(0.8,0.3)],
            [(0.56,0.35),(0.8, axis_height)] ]

ax = fig.axes
for x in range(0,len(ax)):
    params = cent_WH[x]
    W = params[1][0]
    H = params[1][1]
    L = params[0][0]-W/2
    B = params[0][1]-H/2
    ax[x].set_position([L,B,W,H])

for x in range(0,len(ax)):
    ax[x].tick_params(axis='both', which='major', labelsize = min_font, length = ticklength, pad = tickpad)
    ax[x].patch.set_alpha(0.0)
    ax[x].spines['right'].set_visible(False)
    ax[x].spines['top'].set_visible(False)

coord_conn = [[(zoom_lu[0],np.amin(chrom.signal)), (np.amin(xi),np.amax(yi))],
              [(zoom_lu[1],np.amin(chrom.signal)), (np.amax(xi),np.amax(yi))]]

for c in coord_conn:
    con2 = ConnectionPatch(xyA=c[0], coordsA=ax[0].transData,
                           xyB=c[1], coordsB=ax[1].transData,
                           arrowstyle="-", color = '#000000',
                           alpha = 1, linewidth = 1)
    fig.add_artist(con2)

plt.savefig(r'C:\Users\willi\Documents\Dynamic_environment_project\Paper_plots\ExampleChrom.png', dpi = 800, transparent = True)
plt.show()

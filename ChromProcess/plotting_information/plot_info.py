from pathlib import Path

def readTwoColInfo(file, col_idx):
    assignments = {}
    with open(file, 'r') as f:
        for c, line in enumerate(f):
            if c > 0:
                spl = line.strip('\n').split(',')
                assignments[spl[0]] = spl[col_idx]
            else:
                pass
    return assignments

path = Path(__file__)
script_dir = path.parent

plot_params = readTwoColInfo(script_dir/'plotting_parameters.csv', 1)

min_font = float(plot_params['min_font'])
lines = float(plot_params['lines'])
labels = float(plot_params['labels'])
font = float(plot_params['font'])
axiswidth = float(plot_params['axiswidth'])
height = float(plot_params['height'])
width = float(plot_params['width'])
legend_font = float(plot_params['legend_font'])
min_font = float(plot_params['min_font'])
ticklength = float(plot_params['ticklength'])
tickpad = float(plot_params['tickpad'])

Nat_format = readTwoColInfo(script_dir/'Nature_formatting.csv', 1)

col_width = float(Nat_format['col_width'])
two_col_width = float(Nat_format['two_col_width'])
full_page =float(Nat_format['full_page'])
modest_display_item = float(Nat_format['modest_display_item'])
composite_figure = float(Nat_format['composite_figure'])

cluster_grayscale = readTwoColInfo(script_dir/'cluster_colourmap.csv', 2)
cluster_colour_map = readTwoColInfo(script_dir/'cluster_colourmap.csv', 1)

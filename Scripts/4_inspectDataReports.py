import os
from pathlib import Path

from ChromProcess import Classes
import plotly.graph_objects as go

# importing information
system_root = Path('/Users/williamrobinson/Documents/Nijmegen')
storage_stem = system_root/'PrebioticDatabase'
data_folder = storage_stem/'Data/GCMS/FRN'

exp_code = '1A'

root_folder = data_folder/f'{exp_code}/DataReports'

integral_report =  Classes.DataReport(root_folder/f'{exp_code}_GCMS_integral_report.csv')
conc_report =  Classes.DataReport(root_folder/f'{exp_code}_GCMS_concentration_report.csv')

fig = go.Figure()

for c,d in enumerate(integral_report.data):
    fig.add_trace(
    go.Scatter(
                x = integral_report.series_values,
                y = integral_report.data[d],
                line = dict(color = 'rgb(1, 1, 1)', width = 1),
                name = 'data')
    )

fig.layout.xaxis.color = '#000000'
fig.layout.yaxis.color = '#000000'
fig.layout.xaxis.linecolor = '#000000'
fig.layout.yaxis.linecolor = '#000000'

fig.update_layout(title='Data',
                  xaxis_title='[]/[internal standard]',
                  yaxis_title='integral response/ IS respose',
                  plot_bgcolor='rgb(236, 236, 247)',
                  xaxis_showgrid=False, yaxis_showgrid=False)
fig.show()

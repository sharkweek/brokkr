import pandas as pd
import plotly.graph_objects as go
from lasso.dyna import Binout

# get results data
binout = Binout('binout')
energy_fields_m = [x for x in binout.read('matsum') if x[-6:] == 'energy']
energy_fields_g = [x for x in binout.read('glstat') if x[-6:] == 'energy']

# parts[[matsum_id, part_id, title]]
parts = [[0, 1, 'Part 1'],
         [1, 2, 'Part 2'],
         [2, 3, 'Part 3']]

# create and populate dataframe with part data
parts_df =pd.DataFrame(columns=['matsum_id', 'pid', 'title', 'data', 'trace'])
for part in parts:
    for field in energy_fields_m:
        df = pd.DataFrame(
            {'matsum_id': [part[0]],
             'pid': [part[1]],
             'title': [part[2]],
             'data': [field],
             'trace': [go.Scatter(x=binout.read('matsum', 'time'),
                                  y=binout.read('matsum', field).T[part[0]],
                                  name=field[:-7])]}
        )
        parts_df = parts_df.append(df, ignore_index=True)

# add global data to dataframe
for field in energy_fields_g:
    parts_df = parts_df.append(
        pd.DataFrame({
            'matsum_id': [None],
            'pid': [None],
            'title': ["Global"],
            'data': [field],
            'trace': [go.Scatter(x=binout.read('glstat', 'time'),
                                 y=binout.read('glstat', field),
                                 name=field)]
        })
    )

# create figure with traces
layout = go.Layout(title={'text': 'Part Energies',
                          'x': 0.5,
                          'xanchor': 'center'},
                   xaxis={'title': 't (s)',
                          'tickformat': '.0s',
                          'hoverformat': '.3s',
                          'showspikes': True,
                          'spikethickness': 1},
                   yaxis={'title': 'Energy (lbf-in)',
                          'tickformat': '.0s',
                          'hoverformat': '.3s',
                          'showspikes': True,
                          'spikethickness': 1},
                   legend={'title': 'Energy Type'},
                   hovermode='closest',
                   height=500,
                   width=800)

fig = go.Figure(layout=layout)

for trace in parts_df.trace:
    fig.add_trace(trace)

# add dropdown menu to figure
titles = [x for x in parts_df.title.unique()]
buttons = [{'label': x,
            'method': 'update',
            'args': [{'visible': parts_df.title == x}]}
            for x in parts_df.title.unique()]
fig.update_layout(updatemenus=[
    {
        'active': 0,
        'buttons': buttons,
        'direction': "down",
        'x': 0,
        'y': 1.3,
        'showactive': True,
    }
])

# print figure
fig.write_html('part-energies.html')

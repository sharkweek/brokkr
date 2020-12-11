import pandas as pd
import plotly.graph_objects as go
from lasso.dyna import Binout

# get results data
binout = Binout('binout')
velocity_fields = [x for x in binout.read('matsum') if 'rbvelocity' in x]

# Update part list based on the model
# parts[[matsum_id, part_id, title]]
parts = [[0, 1, 'Part 1'],
         [1, 2, 'Part 2'],
         [2, 3, 'Part 3'],
         [3, 4, 'Front Panel'],
         [4, 5, 'Front Ramp'],
         [5, 6, 'Aero Stack'],
         [6, 7, 'Side Buildup'],
         [7, 8, 'Sides'],
         [8, 9, 'Attach Flange'],
         [9, 10, 'Radomes'],
         [10, 11, 'Core 2'],
         [11, 12, 'Core 1'],
         [12, 14, 'Side Buildup Ramp'],
         [13, 15, 'Liner 1'],
         [14, 16, 'Closeout Ribs'],
         [15, 17, 'Liner 2'],
         [16, 18, 'Liner 3']]

# create and populate dataframe with part data
parts_df = pd.DataFrame(
    columns=['matsum_id', 'pid', 'part_title', 'direction', 'trace']
)

for part in parts:
    for field in velocity_fields:
        df = pd.DataFrame(
            {'matsum_id': [part[0]],
             'pid': [part[1]],
             'part_title': [part[2]],
             'direction': [field[0]],
             'trace': [go.Scatter(x=binout.read('matsum', 'time'),
                                  y=binout.read('matsum', field).T[part[0]],
                                  name=field[0])]})
        parts_df = parts_df.append(df, ignore_index=True)

# create resultant velocity
for part in parts:
    # get resultant
    temp_df = pd.DataFrame({'time': binout.read('matsum', 'time')})

    temp_df['resultant'] = 0
    for vel in velocity_fields:
        temp_df['resultant'] += pd.Series(
            binout.read('matsum', vel).T[part[0]]
        )**2

    temp_df['resultant'] = temp_df['resultant']**(0.5)

    # add resultant traces to parts_df
    parts_df = parts_df.append(pd.DataFrame(
        {
            'matsum_id': [part[0]],
            'pid': [part[1]],
            'part_title': [part[2]],
            'direction': ['resultant'],
            'trace': [go.Scatter(x=binout.read('matsum', 'time'),
                                 y=temp_df['resultant'])]
        }
    ), ignore_index=True)

# create figure with traces
layout = go.Layout(title={'text': 'Part Rigid Body Velocities',
                          'x': 0.5,
                          'xanchor': 'center'},
                   xaxis={'title': 't (s)',
                          'tickformat': '.0s',
                          'hoverformat': '.3s',
                          'showspikes': True,
                          'spikethickness': 1},
                   yaxis={'title': 'Velocity (in/s)',
                          'tickformat': '.0s',
                          'hoverformat': '.3s',
                          'showspikes': True,
                          'spikethickness': 1},
                   legend={'title': 'Direction'},
                   hovermode='closest',
                   height=500,
                   width=800)

fig = go.Figure(layout=layout)

for trace in parts_df.trace:
    fig.add_trace(trace)

# add dropdown menu to figure
dropdown_bypart = [{'label': x,
                   'method': 'update',
                   'args': [{'visible': parts_df.part_title == x,
                             'name': parts_df.direction}]}
                  for x in parts_df.part_title.unique()]

dropdown_bydir = [{'label': x,
                  'method': 'update',
                  'args': [{'visible': parts_df.direction == x,
                             'name': parts_df.part_title}]}
                 for x in parts_df.direction.unique()]

fig.update_layout(updatemenus=[{'active': 0,
                                'type': 'dropdown',
                                'buttons': dropdown_bypart,
                                'direction': "down",
                                'x': 0.05,
                                'xanchor': 'left',
                                'y': 1.3,
                                'showactive': True},
                               {'active': 0,
                                'type': 'dropdown',
                                'buttons': dropdown_bydir,
                                'direction': 'down',
                                'x': 0.05,
                                'xanchor': 'left',
                                'y': 1.17,
                                'showactive': True}])

fig.add_annotation({'text': 'By Part: ',
                    'xanchor': 'right',
                    'x': 0.05,
                    'xref': 'paper',
                    'y': 1.28,
                    'yref': 'paper',
                    'showarrow': False})
fig.add_annotation({'text': 'By Direction: ',
                    'xanchor': 'right',
                    'x': 0.05,
                    'xref': 'paper',
                    'y': 1.15,
                    'yref': 'paper',
                    'showarrow': False})

# update initial state of figure
fig.update_traces({'visible': False})

# print figure
fig.write_html(
    binout.read('glstat', 'title').strip().replace(' ', '-').lower()
    + '-part-velocities.html'
)

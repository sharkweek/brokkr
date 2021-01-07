from pandas import DataFrame, Index, Series
from numpy import ndarray
from lasso.dyna import Binout
from plotly.graph_objects import (
    Figure,
    Scatter,
    Layout
)


class Extended_Binout(Binout):

    def glstat(self, *args):
        """`glstat` fields."""

        return self.read('glstat', *args)

    def matsum(self, *args):
        """`matsum` fields."""

        return self.read('matsum', *args)

    def nodfor(self, *args):
        """`nodfor` fields."""

        return self.read('nodfor', *args)

    def rbdout(self, *args):
        """`rbdout` fields."""

        return self.read('rbdout', *args)

    def ssstat(self, *args):
        """`ssstat` fields."""

        return self.read('ssstat', *args)

    def as_df(self, *args):
        """Return data as a pandas DataFrame."""

        data = self.read(*args)

        # validate time-based data
        if type(data) != ndarray:
            raise ValueError("Not a time series")
        elif data.shape[0] != self.read(args[0], 'time').shape[0]:
            raise ValueError("Not a time series")
        else:
            time = Index(self.read(args[0], 'time'), name='time')

        # create dataframe
        if data.ndim > 1:
            df = DataFrame(index=time)
            for i, j in enumerate(self.read(args[0], 'ids')):
                df[str(j)] = data.T[i]
        else:
            df = Series(data, index=time, name=args[-1])

        return df

    @property
    def df(self):
        return self._DF(self)

    class _DF:
        """Inner class containing methods to return data as `DataFrames`."""

        def __init__(self, outer):
            """Reference the outer object for inner methods."""

            self.outer = outer

        @property
        def matsum_legend(self):
            """`matsum` legend as a pandas `DataFrame`"""

            legend = self.outer.read('matsum', 'legend')
            return DataFrame({
                'pid': self.outer.read('matsum', 'ids'),
                'title': [legend[i:i + 80].strip()
                          for i in range(0, len(legend), 80)]
            })

    @property
    def plot(self):
        return self._Plot(self)

    class _Plot:
        """Inner class containing methods to return interactive plots."""

        def __init__(self, outer):
            """Reference the outer object for inner methods."""

            self.outer = outer

        def matsum_energy(self, filename=None):
            """Plot `matsum` energy fields as Plotly figure."""

            energy_fields_m = [
                x for x in self.outer.read('matsum') if x[-6:] == 'energy'
            ]
            energy_fields_g = [
                x for x in self.outer.read('glstat') if x[-6:] == 'energy'
            ]

            # create and populate dataframe with part data
            legend = self.outer.df.matsum_legend
            plot_df = DataFrame(columns=['matsum_id',
                                         'pid',
                                         'part_title',
                                         'type',
                                         'trace'])
            for part in legend.title:
                for field in energy_fields_m:
                    rec = legend[legend.title == part]
                    plot_df = plot_df.append(
                        {"matsum_id": rec.index[0],
                         "pid": rec.pid.iloc[0],
                         "part_title": rec.title.iloc[0]
                                       + " (" + str(rec.pid.iloc[0]) + ")",
                         "type": field[:-7],
                         "trace": Scatter(
                            x=self.outer.read('matsum', 'time'),
                            y=self.outer.read('matsum', field).T[rec.index[0]],
                            name=field[:-7]
                        )},
                        ignore_index=True
                    )

            # add global data to dataframe
            for field in energy_fields_g:
                plot_df = plot_df.append(
                    {'matsum_id': None,
                     'pid': None,
                     'part_title': "Global",
                     'type': field[:-7],
                     'trace': Scatter(x=self.outer.read('glstat', 'time'),
                                      y=self.outer.read('glstat', field),
                                      name=field[:-7])},
                    ignore_index=True
                )

            # create figure with traces
            layout = Layout(title={'text': 'Part Energies',
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

            fig = Figure(layout=layout)

            for trace in plot_df.trace:
                fig.add_trace(trace)

            # add dropdown menu to figure
            dropdown_bypart = [{'label': x,
                                'method': 'update',
                                'args': [{'visible': plot_df.part_title == x,
                                        'name': plot_df.type}]}
                                        for x in plot_df.part_title.unique()]

            dropdown_bydir = [{'label': x,
                               'method': 'update',
                               'args': [{'visible': plot_df.type == x,
                                         'name': plot_df.part_title}]}
                              for x in plot_df.type.unique()]

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
            fig.add_annotation({'text': 'By Type: ',
                                'xanchor': 'right',
                                'x': 0.05,
                                'xref': 'paper',
                                'y': 1.15,
                                'yref': 'paper',
                                'showarrow': False})

            if filename is not None:
                # print figure
                fig.write_html(
                    self.outer.read('glstat',
                                    'title'). strip().replace(' ', '-').lower()
                        + '-part-energies.html',
                    include_mathjax='cdn',
                    default_width=800,
                )
            else:
                return fig

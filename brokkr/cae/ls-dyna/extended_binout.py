from pandas import DataFrame
from lasso.dyna import Binout

class Extended_Binout(Binout):
    def read(self, *args):
        """Extends the `read` functionality of `Binout`"""
        if args[0:2] == ('matsum', 'legend'):
            return self._matsum_legend_df()
        else:
            return super().read(*args)

    def _matsum_legend_df(self):
        """reformat the `matsum` legend as a pandas `DataFrame`"""
        legend = super().read('matsum', 'legend')
        titles = [legend[i:i + 80].strip() for i in range(0, len(legend), 80)]
        parts = [list(x) for x
                 in zip(range(0, len(titles)),
                        super().read('matsum', 'ids'),
                        titles)]
        df = pd.DataFrame({
            'pid': super().read('matsum', 'ids'),
            'title': [legend[i:i + 80].strip()
                      for i in range(0, len(legend), 80)]
        })

        return df

    def matsum_energy_plot(self, filename=None):
        """Plot `matsum` energy fields as Plotly figure."""

        energy_fields_m = [x for x
                           in super().read('matsum') if x[-6:] == 'energy']
        energy_fields_g = [x for x
                           in super().read('glstat') if x[-6:] == 'energy']

        # create and populate dataframe with part data
        legend = self._matsum_legend_df()
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
                     "part_title": rec.title.iloc[0],
                     "type": field[:-7],
                     "trace": Scatter(
                         x=super().read('matsum', 'time'),
                         y=super().read('matsum', field).T[rec.index[0]],
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
                 'trace': Scatter(x=super().read('glstat', 'time'),
                                  y=super().read('glstat', field),
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
                super().read('glstat',
                             'title'). strip().replace(' ', '-').lower()
                + '-part-energies.html',
                include_mathjax='cdn',
                default_width=800,
            )
        else:
            return fig

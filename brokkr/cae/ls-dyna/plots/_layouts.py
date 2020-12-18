"""Default Plotly layouts for LS-DYNA plots"""

from plotly.graph_objects import Layout

DEFAULT_LAYOUT = Layout(title={'x': 0.5,
                               'xanchor': 'center'},
                        xaxis={'tickformat': '.3s'},
                        yaxis={'tickformat': '.3s'},
                        hovermode='x')

MASS_V_TIME_LAYOUT = Layout(title={'text': 'Mass v. Time',
                                   'x': 0.5,
                                   'xanchor': 'center'},
                            xaxis={'title': 'Time',
                                   'tickformat': '.3s'},
                            yaxis={'title': 'Mass',
                                   'tickformat': '.3s'},
                            hovermode='x')

ENERGY_V_TIME_LAYOUT = Layout(title={'text': 'Energy v. Time',
                                     'x': 0.5,
                                     'xanchor': 'center'},
                              xaxis={'title': 'Time',
                                     'tickformat': '.3s'},
                              yaxis={'title': 'Energy',
                                     'tickformat': '.3s'},
                              hovermode='x')

FORCE_V_TIME_LAYOUT = Layout(title={'text': 'Force v. Time',
                                    'x': 0.5,
                                    'xanchor': 'center'},
                             xaxis={'title': 'Time',
                                    'tickformat': '.3s'},
                             yaxis={'title': 'Force',
                                    'tickformat': '.3s'},
                             hovermode='x')

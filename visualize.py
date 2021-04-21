#!/Users/jennysheng/anaconda3/envs/tools/bin/python

import sys
import os
import glob
import numpy as np
import pandas as pd
import gzip
import matplotlib.pyplot as plt
from scipy import stats
from scipy.optimize import curve_fit
import itertools
from collections import Counter
import plotly.express as px
import plotly.graph_objects as go
import plotly.offline as offline
import plotly.subplots as subplots

wt_full = ('MSGFRKMAFPSGKVEGCMVQVTCGTTTLNGLWLDDVVYCPRHVICT'
           'SEDMLNPNYEDLLIRKSNHNFLVQAGNVQLRVIGHSMQNCVLKLKV'
           'DTANPKTPKYKFVRIQPGQTFSVLACYNGSPSGVYQCAMRPNFTIK'
           'GSFLNGSCGSVGFNIDYDCVSFCYMHHMELPTGVHAGTDLEGNFYG'
           'PFVDRQTAQAAGTDTTITVNVLAWLYAAVINGDRWFLNRFTTTLND'
           'FNLVAMKYNYEPLTQDHVDILGPLSAQTGIAVLDMCASLKELLQNG'
           'MNGRTILGSALLEDEFTPFDVVRQCSGVTFQ')
wt_full = [x for x in wt_full]

def make_heatmap(df, x_, y_, wt_, show = True, save = False, **kwarg):
    '''
    Given a list of dataframe generated from
    amino_acid attribute, return an interactive plotly figure.
    _________________
    Input:
    df: dataframe of amino acids at each site.
    x_: x-axis labels
    y_: y-axis labels (list in inverted order)
    wt_: a list of the wildtype residues at each of the amino acid positions
    if save = True add kwarg name: file path for saving figure

    Output: Plotly figure saved at location specified at filename.
    Saves an html and a pdf version.
    '''
    df = df.reindex(y_)
    fig = go.Figure(data = go.Heatmap(z = df,
                                 x = x_, y = df.index, colorscale='RdBu', zmid=0,
                                     zmin = -3))
    #Add marker to denote wildtype
    wt = wt_
    fig.add_scatter(y=wt, x = x_, mode="markers",
                    marker=dict(size=4, color="Black"),
                    name="wt")
    fig.update_layout(xaxis=dict(rangeslider=dict(visible=True),
                             type="linear"), xaxis_showgrid=False,
                             yaxis_showgrid=False)
    fig.layout.font.family = 'Arial'
    fig.update_layout({
        'plot_bgcolor': 'rgba(166, 166, 166, 0.7)',
        'paper_bgcolor': 'rgba(0, 0, 0, 0)',
        })
    if show == True:
        fig.show()
    if save == True:
        offline.plot(fig, filename = kwarg['name']+'.html')
        fig.write_image(kwarg['name']+'.pdf')

def make_heatmap_vert(df, wt_, show = True, save = False, **kwarg):
    '''
    Given a list of dataframe generated from
    amino_acid attribute, return an interactive plotly figure.
    _________________
    Input:
    df: dataframe of amino acids at each site.
    # x_: x-axis labels
    # y_: y-axis labels (list in inverted order)
    wt_: a list of the wildtype residues at each of the amino acid positions
    if save = True add kwarg name: file path for saving figure

    Output: Plotly figure saved at location specified at filename.
    Saves an html and a pdf version.
    '''
    fig = go.Figure(data = go.Heatmap(z = df,
                    x = df.columns, y = [int(x[4:]) for x in list(df.index)],
                    colorscale='RdBu', zmid=0))
    #Add marker to denote wildtype
    wt = wt_
    fig.add_scatter(y=[int(x[4:]) for x in list(df.index)], x = wt_,
                    mode="markers",
                    marker=dict(size=4, color="Black"),
                    name="wt")
    fig.update_layout(
        yaxis=dict(autorange='reversed')
    )
    fig.update_layout(xaxis_showgrid=False,
                      yaxis_showgrid=False,
                      width=600,
                      height=2200,)
    if show == True:
        fig.show()
    if save == True:
        offline.plot(fig, filename = kwarg['name']+'.html')
        fig.write_image(kwarg['name']+'.pdf')

def heatmap_patient(
        df, x_, y_, wt_, position, value, show = True, save = False,**kwarg):
    '''
    Given a list of dataframe generated from
    amino_acid attribute, return an interactive plotly figure.
    Includes sequenced patient clinical variants.
    _________________
    Input:
    df: dataframe of amino acids at each site.
    x_: x-axis labels
    y_: y-axis labels (list in inverted order)
    position: residue position of clinical variant call
    value: amino acid position of clinical variant call
    if save = True add kwarg name: file path for saving figure

    Output: Plotly figure saved at location specified at filename
    Saves an html and a pdf version
    '''
    df = df.reindex(y_)
    fig = go.Figure(data = go.Heatmap(z = df,
                                 x = x_, y = df.index, colorscale='RdBu',
                                 zmid=0))
# Add marker to denote patient variants
    fig.add_scatter(y=value, x = position, mode="markers",
                    marker=dict(size=8, color="Black", symbol='x')
                   )
#     Add marker to denote wildtype
    wt = wt_
    fig.add_scatter(y=wt, x = x_, mode="markers",
                    marker=dict(size=4, color="Black"))


    fig.update_layout({
            'plot_bgcolor': 'rgba(166, 166, 166, 0.7)',
            'paper_bgcolor': 'rgba(0, 0, 0, 0)',
            })
    fig.update_layout(xaxis=dict(rangeslider=dict(visible=True),
                             type="linear"))
    fig.layout.font.family = 'Arial'
    fig.update_layout(showlegend=False, xaxis_showgrid=False,
            yaxis_showgrid=False)

    if show == True:
        fig.show()
    if save == True:
        offline.plot(fig, filename = kwarg['name']+'.html')
        fig.write_image(kwarg['name']+'.pdf')

def heatmap_rsa(
        df, x_, y_, position, grouped_aa,
        rsa, single_muts, show = True, save = False, **kwarg):
    '''
    Given a list of dataframe generated from
    amino_acid attribute, return an interactive plotly figure.
    Includes sequenced patient clinical variants, rsa, and secondary
    structure information.
    _________________
    Input:
    df: dataframe of amino acids at each site.
    x_: x-axis labels
    y_: y-axis labels (list in inverted order)
    wt_: a list of the wildtype residues at each of the amino acid positions
    position: residue position of clinical variant call
    if save = True add kwarg name: file path for saving figure
    rsa: dataframe of rsa information produced by dms-tools2
    single_muts: list of snv as seen in clinical isolates

    Output: Plotly figure saved at location specified at filename
    Saves an html and a pdf version
    '''


    fig = subplots.make_subplots(rows=4, cols=1,
            shared_xaxes=True,
            vertical_spacing=0.02, row_heights=[0.05, 0.05,0.05, 0.95])


    all_res_regroup = df.reindex(grouped_aa)
    fig.add_trace(go.Heatmap(
            z = all_res_regroup,
            x = x_,
            y = all_res_regroup.index,
            colorscale='RdBu', zmid=0,
            zmin = -3),
            row = 4, col = 1)
    fig.update_layout({
            'plot_bgcolor': 'rgba(166, 166, 166, 0.7)',
            'paper_bgcolor': 'rgba(0, 0, 0, 0)',
            },
            xaxis_showgrid=False, yaxis_showgrid=False)

# add RSA info
    fig.add_trace(go.Heatmap(
            z = rsa['RSA'],
            x = rsa['site'], y = ['RSA']*306, colorscale='Blues'),
            row = 2, col = 1)

# add secondary structure info
    fig.add_trace(go.Heatmap(
            z = rsa['SS_num'],
            x = rsa['site'], y = ['SS']*306, colorscale='RdBu'),
            row = 1, col = 1)

# add average score info
    fig.add_trace(go.Heatmap(
            z = df.iloc[::-1].mean(),
            x = x_, y = ['Average']*306, colorscale='RdBu',
            zmid=0,),
            row = 3, col = 1)

# Add marker to denote clinical variants

    fig.add_scatter(
            y=single_muts, x = position, mode="markers",
            marker=dict(size=5, color="Black", symbol='star-triangle-down'),
            row = 4, col = 1
                       )
# Add marker to denote wildtype
    wt = wt_full[1:]
    fig.add_scatter(y=wt, x = x_, mode="markers",
            marker=dict(size=4, color="Black"),
            name="wt", row = 4, col = 1)


    fig.layout.font.family = 'Arial'

    fig.update_layout({
            'plot_bgcolor': 'rgba(166, 166, 166, 0.7)',
            'paper_bgcolor': 'rgba(0, 0, 0, 0)',
            },
            xaxis_showgrid=False, yaxis_showgrid=False)
    fig.update_yaxes(showgrid=False, row=4, col=1)

    if show == True:
        fig.show()
    if save == True:
        offline.plot(fig, filename = kwarg['name']+'.html')
        fig.write_image(kwarg['name']+'.pdf')

def heatmap_rsa_vert(
        df, position, grouped_aa,
        rsa, single_muts, means, show = True, save = False, **kwarg):
    '''
    Given a list of dataframe generated from
    amino_acid attribute, return an interactive plotly figure.
    Includes sequenced patient clinical variants, rsa, and secondary
    structure information.
    _________________
    Input:
    df: dataframe of amino acids at each site.
    x_: x-axis labels
    y_: y-axis labels (list in inverted order)
    wt_: a list of the wildtype residues at each of the amino acid positions
    position: residue position of clinical variant call
    if save = True add kwarg name: file path for saving figure
    rsa: dataframe of rsa information produced by dms-tools2
    single_muts: list of snv as seen in clinical isolates

    Output: Plotly figure saved at location specified at filename
    Saves an html and a pdf version
    '''


    fig = subplots.make_subplots(rows=1, cols=4,
            shared_yaxes=True,
            vertical_spacing=0.02, column_widths=[0.05, 0.05,0.05, 0.95])

    fig.add_trace(go.Heatmap(
            z = df,
            x = df.columns,
            y = [int(x[4:]) for x in list(df.index)],
            colorscale='RdBu', zmid=0,
            zmin = -3),
            row = 1, col = 4)
    fig.update_layout({
            'plot_bgcolor': 'rgba(166, 166, 166, 0.7)',
            'paper_bgcolor': 'rgba(0, 0, 0, 0)',
            },
            xaxis_showgrid=False, yaxis_showgrid=False)

# add RSA info
    fig.add_trace(go.Heatmap(
            z = rsa['RSA'],
            y = rsa['site'], x = ['RSA']*306, colorscale='Blues'),
            row = 1, col = 2)

# add secondary structure info
    fig.add_trace(go.Heatmap(
            z = rsa['SS_num'],
            y = rsa['site'], x = ['SS']*306, colorscale='RdBu'),
            row = 1, col = 1)

# add average score info
    fig.add_trace(go.Heatmap(
            z = means,
            y = [int(x[4:]) for x in list(df.index)],
            x = ['Average']*306, colorscale='RdBu',
            zmid=0,),
            row = 1, col = 3)

# Add marker to denote clinical variants

    fig.add_scatter(
            y=position, x = single_muts, mode="markers",
            marker=dict(size=5, color="Black", symbol='star-triangle-down'),
            row = 1, col = 4
                       )
# Add marker to denote wildtype
    wt = wt_full[1:]
    fig.add_scatter(x=wt, y = [int(x[4:]) for x in list(df.index)],
            mode="markers",
            marker=dict(size=4, color="Black"),
            name="wt", row = 1, col = 4)


    fig.layout.font.family = 'Arial'

    fig.update_layout({
            'plot_bgcolor': 'rgba(166, 166, 166, 0.7)',
            'paper_bgcolor': 'rgba(0, 0, 0, 0)',
            },
            xaxis_showgrid=False,
            yaxis_showgrid=False,
            width=900,
            height=2200,
            yaxis=dict(autorange='reversed'))
    fig.update_layout(
    yaxis = dict(
        tickmode = 'linear',
        dtick = 10
    )
)

    # fig.update_layout(yaxis=dict(autorange='reversed'), row = 1, column = 1)

    fig.update_yaxes(showgrid=False, row=1, col=4)

    if show == True:
        fig.show()
    if save == True:
        offline.plot(fig, filename = kwarg['name']+'.html')
        fig.write_image(kwarg['name']+'.pdf')

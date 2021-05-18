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
from plotly.subplots import make_subplots

wt_full = ('MSGFRKMAFPSGKVEGCMVQVTCGTTTLNGLWLDDVVYCPRHVICT'
           'SEDMLNPNYEDLLIRKSNHNFLVQAGNVQLRVIGHSMQNCVLKLKV'
           'DTANPKTPKYKFVRIQPGQTFSVLACYNGSPSGVYQCAMRPNFTIK'
           'GSFLNGSCGSVGFNIDYDCVSFCYMHHMELPTGVHAGTDLEGNFYG'
           'PFVDRQTAQAAGTDTTITVNVLAWLYAAVINGDRWFLNRFTTTLND'
           'FNLVAMKYNYEPLTQDHVDILGPLSAQTGIAVLDMCASLKELLQNG'
           'MNGRTILGSALLEDEFTPFDVVRQCSGVTFQ')

def original_dist(
        folder, suffix, samples, sets, res_redo,
        all_sets, save = False, **kwarg):
    '''
    Distribution shape of original scores from screen.
    _______________
    Input:
    folder: column name in sample spreadsheet that points to folder
    suffix: suffix of the file name
    if save = True add kwarg name: file path for saving figure
    '''

    fig = make_subplots(
    rows=5, cols=6)

    layout= itertools.product(range(1,6), range(1,7))
    for x, pos in list(zip(sets + res_redo, layout)):
        # old replicates
        if x in sets:
            fchange = pd.read_csv(list(samples[samples['Set']==\
                str(x)][folder])[0]\
                +str(x) + suffix, index_col = [0])
            flat_list = [item for sublist in fchange.values for item in sublist]
            fig.add_trace(go.Histogram(x=flat_list,
                                              xbins=dict(# bins for histogram
                    start=min(flat_list),
                    end=max(flat_list),
                    size=0.25
                ),), row=pos[0], col=pos[1])

        # new replicates single residues
        elif x in res_redo:
            start = list(samples[samples['Set'] == str(x)]['Start range'])[0]
            end = list(samples[samples['Set'] == str(x)]['End range'])[0]
            sites = list(samples[samples['Set'] == str(x)]['Sites'])[0]
            sites = ['Res '+ str(x) for x in sites.split(',')]
            fchange = pd.read_csv(list(samples[samples['Set']==\
                str(x)][folder])[0]\
                +str(x) +suffix, index_col = [0])
            fchange.columns = ['Res '+str(x) for x in list(range(start, end))]
            fchange = fchange[sites]
            flat_list = [item for sublist in fchange.values for item in sublist]
            fig.add_trace(go.Histogram(x=flat_list,
                                              xbins=dict( # bins for histogram
                    start=min(flat_list),
                    end=max(flat_list),
                    size=0.25
                ),), row=pos[0], col=pos[1])

    fig.update_layout(height=700, width=900,
                      title_text=kwarg['title'])
    fig.show()
    if save == True:
        plotly.offline.plot(fig, filename = kwarg['name'])

def transform_sigma(folder, suffix, samples, sets, res_redo, all_sets):
    '''
    Takes the folder and the suffix of the files and
    computes the standard deviations--scaling factors--for
    all sets.
    ____________
    Input:
    folder--name of category in the sample spreadsheet that points to
    folder where data are stored
    suffix--suffix of the file name
    Output:
    list in which each element is paired. First item in pair is
    the set, second is the numerical value for the standard deviation
    --scaling factor--of the set.
    '''
    sigma_list = []
    for x in all_sets:
        fchange = pd.read_csv(list(samples[samples['Set']==str(x)][folder])[0]\
                              + str(x) + suffix, index_col = [0])
        start = list(samples[samples['Set'] == str(x)]['Start range'])[0]
        end = list(samples[samples['Set'] == str(x)]['End range'])[0]
        fchange.columns = ['Res '+str(x) for x in list(range(start, end))]
        # name the columns
        wt_subseq = wt_full[start:end] #find WT residues for the set

        flat_list = np.array([item for sublist in fchange.values for item \
            in sublist])
        mean = flat_list[~np.isnan(flat_list)].mean() # mean of the set
        var = flat_list[~np.isnan(flat_list)].var() # variance of the set
        # set the variance of all set to 1
        var_norm = (flat_list-mean)/np.sqrt(var)+mean
        sigma_list.append([x, np.sqrt(var)])
    return(sigma_list)

def transform_dist(
        folder, suffix, samples,
        sets, res_redo, all_sets, set21, save = False, **kwarg):
    '''
    Distribution shape of original scores from screen.
    ___________________
    Input:
    folder: column name in sample spreadsheet that points to folder
    suffix: suffix of the file name
    if save = True add kwarg name: file path for saving figure

    '''
    # average value of wt
    # average value of stop codon
    # number of residues in set
    mean_stop = {}
    len_set = {}

    fig = make_subplots(
    rows=5, cols=6)
    layout= itertools.product(range(1,6), range(1,7))
    for x, pos in list(zip(sets + res_redo + set21, layout)):
        if x in sets:
            fchange = pd.read_csv(list(samples[samples['Set']==\
                str(x)][folder])[0]\
                + str(x) + suffix, index_col = [0])
            start = list(samples[samples['Set'] == str(x)]['Start range'])[0]
            end = list(samples[samples['Set'] == str(x)]['End range'])[0]
            fchange.columns = ['Res '+str(x) for x in list(range(start, end))]
            # name the columns
            wt_subseq = wt_full[start:end] #find WT residues for the set

            fchange.columns = ['Res '+str(x) for x in list(range(start, end))]

            #set average wt to 0
            cols = fchange.columns
            wt_vals = []
            for row, col in zip(wt_subseq, cols):
                wt_vals.append(fchange.loc[row, col])

            wt_mean = np.mean(wt_vals)
            fchange = fchange - wt_mean
            # add to dict for mean stop
            mean_stop[str(x)] = np.mean(fchange.loc['*'])
            len_set[str(x)] = len(fchange.columns)

            stop_mean = np.mean(fchange.loc['*'])
            scale_factor = -1/stop_mean
            fchange_norm = (fchange - wt_mean)*scale_factor

            flatten_fchange = fchange_norm.values
            flat_list = np.array([item for sublist in flatten_fchange for\
                item in sublist])

            fig.add_trace(go.Histogram(x=flat_list,
                                              xbins=dict( # bins for histogram
                    start=min(flat_list),
                    end=max(flat_list),
                    size=0.25
                ),), row=pos[0], col=pos[1])

        elif x in set21:
            fchange = pd.read_csv(list(samples[samples['Set']==\
                str(x)][folder])[0]\
                + str(x) + suffix, index_col = [0])
            start = list(samples[samples['Set'] == str(x)]['Start range'])[0]
            end = list(samples[samples['Set'] == str(x)]['End range'])[0]
            fchange.columns = ['Res '+str(x) for x in list(range(start, end))]
            # name the columns
            wt_subseq = wt_full[start:end] #find WT residues for the set

            fchange.columns = ['Res '+str(x) for x in list(range(start, end))]

            #set average wt to 0
            cols = fchange.columns[:2]
            wt_vals = []
            for row, col in zip(wt_subseq, cols):
                wt_vals.append(fchange.loc[row, col])

            wt_mean = np.mean(wt_vals)
            fchange = fchange - wt_mean
            # add to dict for mean stop
            mean_stop[str(x)] = np.mean(fchange.loc['*'][:2])
            len_set[str(x)] = 2

            stop_mean = np.mean(fchange.loc['*'][:2])
            scale_factor = -1/stop_mean
            fchange_norm = (fchange - wt_mean)*scale_factor

            flatten_fchange = fchange_norm.values
            flat_list = np.array([item for sublist in flatten_fchange for\
                item in sublist])

            fig.add_trace(go.Histogram(x=flat_list,
                                              xbins=dict( # bins for histogram
                    start=min(flat_list),
                    end=max(flat_list),
                    size=0.25
                ),), row=pos[0], col=pos[1])

        else: # for all individually repeated residues
            set_ind = x.find('R') #identify the R notation for the repeated set
            set_redo = x[:set_ind]
            start = list(samples[samples['Set'] == str(x)]['Start range'])[0]
            end = list(samples[samples['Set'] == str(x)]['End range'])[0]
            sites_ = list(samples[samples['Set'] == str(x)]['Sites'])[0]
            sites = ['Res '+ str(x) for x in sites_.split(',')]
            fchange = pd.read_csv(list(samples[samples['Set']==\
                str(x)][folder])[0]\
                + str(x) + suffix, index_col = [0])
            fchange.columns = ['Res '+str(x) for x in list(range(start, end))]
            fchange = fchange[sites]
            #find WT residues for the set
            wt_subseq = [wt_full[int(ind)] for ind in sites_.split(',')]
            cols = fchange.columns
            wt_vals = []
            for row, col in zip(wt_subseq, cols):
                wt_vals.append(fchange.loc[row, col])

            # Calculate scaling values for slotting in individual residues
            wt_mean = np.mean(wt_vals)
            fchange = fchange-wt_mean
            stop_mean = np.mean(fchange.loc['*'])
            scale_factor = mean_stop[set_redo]/stop_mean
            fchange_norm = fchange*scale_factor

    #         print(x, np.mean(wt_vals), np.var(wt_vals))
            flatten_fchange = fchange_norm.values
            flat_list = np.array([item for sublist in flatten_fchange for \
                item in sublist])

            fig.add_trace(go.Histogram(x=flat_list,
                            xbins=dict( # bins used for histogram
                    start=min(flat_list),
                    end=max(flat_list),
                    size=0.25
                ),), row=pos[0], col=pos[1])

    fig.update_layout(height=700, width=900,
                      title_text=kwarg['title'])
    fig.show()
    if save == True:
        plotly.offline.plot(fig, filename = kwarg['name'])

def transform_dist_sigma(
        folder, suffix, samples,
            sets, res_redo, all_sets, save = False, **kwarg):
    '''
    Distribution shape of original scores from screen. WT from each set
    is set to 0 and variance of each set is set to 1.
    ___________________
    Input:
    folder: column name in sample spreadsheet that points to folder
    suffix: suffix of the file name
    samples: dataframe of sample_spreadsheet with data specs
    sets: complete sets
    res_redo: residues that were individually sequenced
    all_sets: all sets including those that were individually resequenced
    if save = True add kwarg name: file path for saving figure


    '''
    # average value of wt
    # average value of stop codon
    # number of residues in set
    mean_stop = {}
    len_set = {}

    fig = make_subplots(
    rows=5, cols=6)
    layout= itertools.product(range(1,6), range(1,7))
    for x, pos in list(zip(sets + res_redo, layout)):
        if x in sets:
            fchange = pd.read_csv(list(samples[samples['Set']\
                                ==str(x)][folder])[0]\
                                + str(x) + suffix, index_col = [0])
            start = list(samples[samples['Set'] == str(x)]['Start range'])[0]
            end = list(samples[samples['Set'] == str(x)]['End range'])[0]
            # name the columns
            fchange.columns = ['Res '+str(x) for x in list(range(start, end))]
            wt_subseq = wt_full[start:end] #find WT residues for the set

            flat_list = np.array([item for sublist in
                    fchange.values for item in sublist])
            mean = flat_list[~np.isnan(flat_list)].mean() # mean of the set
            var = flat_list[~np.isnan(flat_list)].var() # variance of the set

            # set the variance of all sets to 1
            fchange_norm = (fchange-mean)/np.sqrt(var) + mean
            fchange_norm.columns = ['Res '+str(x) for x
                    in list(range(start, end))] # name the columns

            #set average wt to 0
            cols = fchange_norm.columns
            wt_vals = []
            for row, col in zip(wt_subseq, cols):
                wt_vals.append(fchange_norm.loc[row, col])
            fchange_norm = fchange_norm - np.mean(wt_vals)
            # add to dict for mean stop
            mean_stop[str(x)] = np.mean(fchange_norm.loc['*'])
            len_set[str(x)] = len(fchange_norm.columns)
            flatten_fchange = fchange_norm.values
            flat_list = np.array([item for sublist in
                    flatten_fchange for item in sublist])

            fig.add_trace(go.Histogram(x=flat_list,
                                        xbins=dict( # bins used for histogram
                    start=min(flat_list),
                    end=max(flat_list),
                    size=0.25
                    ),), row=pos[0], col=pos[1])

        else: # for all individually repeated residues
            set_ind = x.find('R') #identify the R notation for the repeated set
            set_redo = x[:set_ind]
            start = list(samples[samples['Set'] == str(x)]['Start range'])[0]
            end = list(samples[samples['Set'] == str(x)]['End range'])[0]
            sites_ = list(samples[samples['Set'] == str(x)]['Sites'])[0]
            sites = ['Res '+ str(x) for x in sites_.split(',')]
            fchange = pd.read_csv(list(samples[samples['Set']\
                                ==str(x)][folder])[0]\
                                + str(x) + suffix, index_col = [0])
            fchange.columns = ['Res '+str(x) for x in list(range(start, end))]
            fchange = fchange[sites]
            wt_subseq = [wt_full[int(ind)] for ind in sites_.split(',')]
            cols = fchange.columns
            wt_vals = []
            for row, col in zip(wt_subseq, cols):
                wt_vals.append(fchange.loc[row, col])

            # Calculate scaling values for slotting in individual residues
            wt_mean = np.mean(wt_vals)
            stop_mean = np.mean(fchange.loc['*'])
            scale_factor = mean_stop[set_redo]/stop_mean
            fchange_norm = (fchange - wt_mean)*scale_factor

    #         print(x, np.mean(wt_vals), np.var(wt_vals))
            flatten_fchange = fchange_norm.values
            flat_list = np.array([item for sublist in\
                    flatten_fchange for item in sublist])

            fig.add_trace(go.Histogram(x=flat_list,
                    xbins=dict( # bins used for histogram
                    start=min(flat_list),
                    end=max(flat_list),
                    size=0.25
                ),), row=pos[0], col=pos[1])

    fig.update_layout(height=700, width=900,
                      title_text=kwarg['title'])
    fig.show()
    if save == True:
        plotly.offline.plot(fig, filename = kwarg['name'])

def transform_dist_mat(folder, suffix, samples, sets, res_redo, all_sets):
    '''
    Distribution shape of original scores from screen.
    ___________________
    Input:
    folder: column name in sample spreadsheet that points to folder
    suffix: suffix of the file name
    if save = True add kwarg name: file path for saving figure

    '''
    # average value of wt
    # average value of stop codon
    # number of residues in set
    len_set = {}
    df_list = []
    mean_stop = {}

    for x in sets + res_redo:
        if x in sets:
            fchange = pd.read_csv(list(samples[samples['Set']==\
                str(x)][folder])[0]\
                + str(x) + suffix, index_col = [0])
            start = list(samples[samples['Set'] == str(x)]['Start range'])[0]
            end = list(samples[samples['Set'] == str(x)]['End range'])[0]
            fchange.columns = ['Res '+str(x) for x in list(range(start, end))]
            wt_subseq = wt_full[start:end] #find WT residues for the set

            flat_list = np.array([item for sublist in fchange.values for\
                item in sublist])
            mean = flat_list[~np.isnan(flat_list)].mean() # mean of the set
            var = flat_list[~np.isnan(flat_list)].var() # variance of the set

            fchange.columns = ['Res '+str(x) for x in list(range(start, end))]

            #set average wt to 0
            cols = fchange.columns
            wt_vals = []
            for row, col in zip(wt_subseq, cols):
                wt_vals.append(fchange.loc[row, col])

            wt_mean = np.mean(wt_vals)
            fchange = fchange - wt_mean
            mean_stop[str(x)] = np.mean(fchange.loc['*']) # add to mean stop
            len_set[str(x)] = len(fchange.columns)

            stop_mean = np.mean(fchange.loc['*'])
            scale_factor = 1/stop_mean
            fchange_norm = (fchange - wt_mean)*scale_factor

            flatten_fchange = fchange_norm.values
            df_list.append(fchange_norm)

        else: # for all individually repeated residues
            set_ind = x.find('R') #identify the R notation for the repeated set
            set_redo = x[:set_ind]
            start = list(samples[samples['Set'] == str(x)]['Start range'])[0]
            end = list(samples[samples['Set'] == str(x)]['End range'])[0]
            sites_ = list(samples[samples['Set'] == str(x)]['Sites'])[0]
            sites = ['Res '+ str(x) for x in sites_.split(',')]
            fchange = pd.read_csv(list(samples[samples['Set']==\
                str(x)][folder])[0]\
                + str(x) + suffix, index_col = [0])
            fchange.columns = ['Res '+str(x) for x in list(range(start, end))]
            fchange = fchange[sites]
            #find WT residues for the set
            wt_subseq = [wt_full[int(ind)] for ind in sites_.split(',')]
            cols = fchange.columns
            wt_vals = []
            for row, col in zip(wt_subseq, cols):
                wt_vals.append(fchange.loc[row, col])

            # Calculate scaling values for slotting in individual residues
            wt_mean = np.mean(wt_vals)
            stop_mean = np.mean(fchange.loc['*'])
            scale_factor = mean_stop[set_redo]/stop_mean
            fchange_norm = (fchange - wt_mean)*scale_factor
            df_list.append(fchange_norm)
    return df_list

def transform_matrix(folder, suffix, samples, sets, res_redo, all_sets, set21):
    '''
    Transform each set so that WT fixed at 0 and stop codon is normalized to
    -1 in each set.
    __________
    Input:
    folder: column name in sample spreadsheet that points to folder
    suffix: suffix of the file name
    samples: sample spreadsheet
    sets: all complete sets
    res_redo: all invividually resequenced sets
    all_sets: all sets
    set21: set21--treated separately because of the C terminus
    '''
    mean_stop = {}
    len_set = {}
    set_list = []
    for file in sets:
        fchange = pd.read_csv(list(samples[samples['Set']==\
            str(file)][folder])[0]\
            +str(file) + suffix, index_col = [0])
        start = list(samples[samples['Set'] == str(file)]['Start range'])[0]
        end = list(samples[samples['Set'] == str(file)]['End range'])[0]
        fchange.columns = ['Res '+str(x) for x in list(range(start, end))]
        wt_subseq = wt_full[start:end] #find WT residues for the set

        flat_list = np.array([item for sublist in fchange.values\
            for item in sublist])
        mean = flat_list[~np.isnan(flat_list)].mean() # mean of the set
        var = flat_list[~np.isnan(flat_list)].var() # variance of the set

        fchange.columns = ['Res '+str(x) for x in list(range(start, end))]

        #set average wt to 0
        cols = fchange.columns
        wt_vals = []
        for row, col in zip(wt_subseq, cols):
            wt_vals.append(fchange.loc[row, col])

        wt_mean = np.mean(wt_vals)
        fchange = fchange - wt_mean
        mean_stop[str(file)] = np.mean(fchange.loc['*'])
        len_set[str(file)] = len(fchange.columns)

        stop_mean = np.mean(fchange.loc['*'])
        scale_factor = -1/stop_mean
        fchange_norm = fchange*scale_factor

        set_list.append(fchange_norm)

    for file in set21:
        fchange = pd.read_csv(list(samples[samples['Set']==\
            str(file)][folder])[0]\
            + str(file) + suffix, index_col = [0])
        start = list(samples[samples['Set'] == str(file)]['Start range'])[0]
        end = list(samples[samples['Set'] == str(file)]['End range'])[0]
        # name the columns
        wt_subseq = wt_full[start:end] #find WT residues for the set

        fchange.columns = ['Res '+str(x) for x in list(range(start, end))]

        #set average wt to 0
        cols = fchange.columns[:2]
        wt_vals = []
        for row, col in zip(wt_subseq, cols):
            wt_vals.append(fchange.loc[row, col])

        wt_mean = np.mean(wt_vals)
        fchange = fchange - wt_mean
        # add to dict for mean stop
        mean_stop[str(file)] = np.mean(fchange.loc['*'][:2])
        len_set[str(file)] = 2

        stop_mean = np.mean(fchange.loc['*'][:2])
        scale_factor = -1/stop_mean
        fchange_norm = (fchange - wt_mean)*scale_factor

        set_list.append(fchange_norm)

    set_list_res = []
    for file in res_redo:
        set_ind = file.find('R') #identify the R notation for the repeated set
        set_redo = file[:set_ind]
        start = list(samples[samples['Set'] == str(file)]['Start range'])[0]
        end = list(samples[samples['Set'] == str(file)]['End range'])[0]
        sites_ = list(samples[samples['Set'] == str(file)]['Sites'])[0]
        sites = ['Res '+ str(x) for x in sites_.split(',')]
        fchange = pd.read_csv(list(samples[samples['Set']==\
            str(file)][folder])[0]\
            + str(file) + suffix, index_col = [0])
        fchange.columns = ['Res '+str(x) for x in list(range(start, end))]
        fchange = fchange[sites]
        wt_subseq = [wt_full[int(ind)] for ind in sites_.split(',')]
        cols = fchange.columns
        wt_vals = []
        for row, col in zip(wt_subseq, cols):
            wt_vals.append(fchange.loc[row, col])

        # Calculate scaling values for slotting in individual residues
        wt_mean = np.mean(wt_vals)
        fchange = fchange-wt_mean
        stop_mean = np.mean(fchange.loc['*'])
        scale_factor = -1/stop_mean
        fchange_norm = fchange *scale_factor


        set_list_res.append(fchange_norm)
    all_residues = pd.concat(set_list, axis = 1)

    all_res_redo = pd.concat(set_list_res, axis = 1)
    all_res_redo = all_res_redo.fillna('NaN')
    all_residues.update(all_res_redo)
    order = ['Res '+str(x) for x in range(1, 307)]
    all_residues = all_residues[order]
    all_residues = all_residues.applymap(lambda x: x if not \
        isinstance(x, str) else np.nan)
    return(all_residues)

def raw_dist(folder, samples, sets, res_redo, all_sets):
    '''
    Returns the mean and standard error of the raw data.
    __________
    Input:
    folder--folder where the datasets are stored. points
    to column in sample spreadsheet (string)
    Output:
    melted dataframe with residue and mutation along with raw
    means and standard deviations
    '''
    fchange_list = []
    redo_list = []
    for x in sets + res_redo:
        # old replicates
        start = list(samples[samples['Set']==str(x)]['Start range'])[0]
        end = list(samples[samples['Set']==str(x)]['End range'])[0]
        sites = list(samples[samples['Set']==str(x)]['Sites'])[0]
        directory = list(samples[samples['Set']==str(x)][folder])[0]
        if x in sets:
            for y in range(start, end):
                fchange = pd.read_csv(directory + '/set' + str(x) + \
                    '_residue' + str(y) + '.csv')
                fchange['residue'] = [y]*len(fchange)
                fchange_list.append(fchange)

        # new replicates single residues

        elif x in res_redo:
            sites = [str(x) for x in sites.split(',')]
            for y in sites:
                fchange = pd.read_csv(directory + '/set' + str(x) + \
                    '_residue' + y + '.csv')
                fchange['residue'] = [y]*len(fchange)
                fchange_list.append(fchange)
                redo_list.append(fchange)

    #list of residues and amino acids along with raw mean and standard error
    error = pd.concat(fchange_list)
    mid = [x[1] for x in error['Translation']]
    error['middle'] = mid
    return(error)

def transform_matrix_sigma(folder,
        suffix, samples, sets, res_redo, all_sets):
    '''
    Transforms the data by set such that all the wildtypes are fixed at zero
    and standard deviation of each set is set to 1.
    __________
    Input:
    folder: column name in sample spreadsheet that points to folder
    suffix: suffix of the file name
    samples: sample spreadsheet
    sets: all complete sets
    res_redo: all invividually resequenced sets
    all_sets: all sets
    '''

    mean_stop = {}
    len_set = {}
    set_list = []
    for file in sets:
        fchange = pd.read_csv(list(samples[samples['Set']
                            ==str(file)][folder])[0]\
                            +str(file) + suffix, index_col = [0])
        start = list(samples[samples['Set'] == str(file)]['Start range'])[0]
        end = list(samples[samples['Set'] == str(file)]['End range'])[0]
        fchange.columns = ['Res '+str(x) for x in list(range(start, end))]
        wt_subseq = wt_full[start:end] #find WT residues for the set

        flat_list = np.array([item for sublist in fchange.values \
                for item in sublist])
        mean = flat_list[~np.isnan(flat_list)].mean() # mean of the set
        var = flat_list[~np.isnan(flat_list)].var() # variance of the set
        # set the variance of all set to 1
        # normalize the set to unit variance
        fchange_norm = (fchange-mean)/np.sqrt(var)+mean
        #new label for columns
        fchange_norm.columns = ['Res '+str(x) for x in list(range(start, end))]
        #figure out average wt values for set and linear transform
        cols = fchange_norm.columns
        wt_vals = []
        for row, col in zip(wt_subseq, cols):
            wt_vals.append(fchange_norm.loc[row, col])
        fchange_norm = fchange_norm - np.mean(wt_vals)
        # add to dict for mean stop
        mean_stop[str(file)] = np.mean(fchange_norm.loc['*'])
        len_set[str(file)] = len(fchange_norm.columns)

        set_list.append(fchange_norm)

    set_list_res = []
    for file in res_redo:
        set_ind = file.find('R') #identify the R notation for the repeated set
        set_redo = file[:set_ind]
        start = list(samples[samples['Set'] == str(file)]['Start range'])[0]
        end = list(samples[samples['Set'] == str(file)]['End range'])[0]
        sites_ = list(samples[samples['Set'] == str(file)]['Sites'])[0]
        sites = ['Res '+ str(x) for x in sites_.split(',')]
        fchange = pd.read_csv(list(samples[samples['Set']\
                            ==str(file)][folder])[0]\
                            + str(file) + suffix, index_col = [0])
        fchange.columns = ['Res '+str(x) for x in list(range(start, end))]
        fchange = fchange[sites]
        #find WT residues for the set
        wt_subseq = [wt_full[int(ind)] for ind in sites_.split(',')]
        cols = fchange.columns
        wt_vals = []
        for row, col in zip(wt_subseq, cols):
            wt_vals.append(fchange.loc[row, col])

        # Calculate scaling values for slotting in individual residues
        wt_mean = np.mean(wt_vals)
        stop_mean = np.mean(fchange.loc['*'])
        scale_factor = mean_stop[set_redo]/(stop_mean-wt_mean)
        fchange_norm = (fchange - wt_mean)*scale_factor

        set_list_res.append(fchange_norm)
    all_residues = pd.concat(set_list, axis = 1)

    all_res_redo = pd.concat(set_list_res, axis = 1)
    all_res_redo = all_res_redo.fillna('NaN')
    all_residues.update(all_res_redo)
    order = ['Res '+str(x) for x in range(1, 307)]
    all_residues = all_residues[order]
    all_residues = all_residues.applymap(lambda x: x\
            if not isinstance(x, str) else np.nan)
    return(all_residues, mean_stop)

def replicate(rep, replicate_folder, cond_suffix, samples,
        sets, res_redo, set21):
    '''
    Tranform the raw foldchanges fro single biological replicates.
    __________
    rep: int denoting replicate number
    replicate folder: column in sample spreadsheet containing
        replicate info
    cond_suffix: file suffix for replicate files
    samples: sample spreadsheet
    sets: all complete sequencing sets
    res_redo: all individually resequenced all_residues
    set21: set 21 to be treated specially for the C terminal portion

    Output: dataframe with score at each amino acid at each residue
    '''
    mean_stop = {}
    len_set = {}
    rep1_set = []
    set_list_res = []
    # replicate 1
    for file in sets + res_redo + set21:
        replicate_dir = list(samples[samples['Set'] == \
            str(file)][replicate_folder])[0]
        if file in sets:
            fchange = pd.read_csv(replicate_dir + str(file)
                    + '_replicate'+str(rep)+cond_suffix, index_col = [0])
            start = list(samples[samples['Set'] == str(file)]['Start range'])[0]
            end = list(samples[samples['Set'] == str(file)]['End range'])[0]
            fchange.columns = ['Res '+str(x) for x in list(range(start, end))]
            wt_subseq = wt_full[start:end]
            cols = fchange.columns
            wt_vals = []
            for row, col in zip(wt_subseq, cols):
                wt_vals.append(fchange.loc[row, col])

            wt_mean = np.mean(wt_vals)
            fchange = fchange - wt_mean
            mean_stop[str(file)] = np.mean(fchange.loc['*'])
            len_set[str(file)] = len(fchange.columns)

            stop_mean = np.mean(fchange.loc['*'])
            scale_factor = -1/stop_mean
            fchange_norm = fchange*scale_factor
            rep1_set.append(fchange_norm)
        elif file in set21:
            fchange = pd.read_csv(replicate_dir + str(file)
                    + '_replicate'+str(rep)+cond_suffix, index_col = [0])
            start = list(samples[samples['Set'] == str(file)]['Start range'])[0]
            end = list(samples[samples['Set'] == str(file)]['End range'])[0]
            fchange.columns = ['Res '+str(x) for x in list(range(start, end))]
            # name the columns
            wt_subseq = wt_full[start:end] #find WT residues for the set

            #set average wt to 0
            cols = fchange.columns[:2]
            wt_vals = []
            for row, col in zip(wt_subseq, cols):
                wt_vals.append(fchange.loc[row, col])

            wt_mean = np.mean(wt_vals)
            fchange = fchange - wt_mean
            # add to dict for mean stop
            mean_stop[str(file)] = np.mean(fchange.loc['*'][:2])
            len_set[str(file)] = 2

            stop_mean = np.mean(fchange.loc['*'][:2])
            scale_factor = -1/stop_mean
            fchange_norm = (fchange - wt_mean)*scale_factor
            rep1_set.append(fchange_norm)

        elif file in res_redo:
            set_ind = file.find('R') #identify the R notation for the repeated set
            set_redo = file[:set_ind]
            start = list(samples[samples['Set'] == str(file)]['Start range'])[0]
            end = list(samples[samples['Set'] == str(file)]['End range'])[0]
            sites_ = list(samples[samples['Set'] == str(file)]['Sites'])[0]
            sites = ['Res '+ str(x) for x in sites_.split(',')]
            fchange = pd.read_csv(replicate_dir + str(file) + '_replicate'\
                    + str(rep) + cond_suffix, index_col = [0])
            fchange.columns = ['Res '+str(x) for x in list(range(start, end))]
            fchange = fchange[sites]
            wt_subseq = [wt_full[int(ind)] for ind in sites_.split(',')]
            cols = fchange.columns
            wt_vals = []
            for row, col in zip(wt_subseq, cols):
                wt_vals.append(fchange.loc[row, col])

            # Calculate scaling values for slotting in individual residues
            wt_mean = np.mean(wt_vals)
            stop_mean = np.mean(fchange.loc['*'])
            scale_factor = -1/stop_mean
            fchange_norm = fchange *scale_factor
            set_list_res.append(fchange_norm)

    all_residues = pd.concat(rep1_set, axis = 1)
    all_res_redo = pd.concat(set_list_res, axis = 1)
    all_res_redo = all_res_redo.fillna('NaN')
    all_residues.update(all_res_redo)
    order = ['Res '+ str(x) for x in range(1, 307)]
    all_residues = all_residues[order]
    all_residues = all_residues.applymap(lambda x: x if not \
            isinstance(x, str) else np.nan)
    return(all_residues)

def replicate_sigma(rep, replicate_folder, cond_suffix, samples,
        sets, res_redo):
    '''
    Tranform the raw foldchanges fro single biological replicates.
    __________
    rep: int denoting replicate number
    replicate folder: column in sample spreadsheet containing
        replicate info
    cond_suffix: file suffix for replicate files
    samples: sample spreadsheet
    sets: all complete sequencing sets
    res_redo: all individually resequenced all_residues
    set21: set 21 to be treated specially for the C terminal portion

    Output: dataframe with score at each amino acid at each residue
    '''
    mean_stop = {}
    len_set = {}
    rep1_set = []
    set_list_res = []
    # replicate 1
    for file in sets + res_redo:
        replicate_dir = list(samples[samples['Set']
                ==str(file)][replicate_folder])[0]
        if file in sets:
            fchange = pd.read_csv(replicate_dir + str(file)\
                    + '_replicate'+str(rep)+cond_suffix, index_col = [0])
            start = list(samples[samples['Set'] == str(file)]['Start range'])[0]
            end = list(samples[samples['Set'] == str(file)]['End range'])[0]
            fchange.columns = ['Res '+str(x) for x in list(range(start, end))]
            wt_subseq = wt_full[start:end] #find WT residues for the set

            flat_list = np.array([item for sublist \
                    in fchange.values for item in sublist])
            mean = flat_list[~np.isnan(flat_list)].mean()
            var = flat_list[~np.isnan(flat_list)].var()
            # set the variance of all set to 1
            var_norm = (flat_list-mean)/np.sqrt(var)+mean

            fchange_norm = (fchange-mean)/np.sqrt(var) + mean
            # add to dict for mean stop
            mean_stop[str(file)] = np.mean(fchange_norm.loc['*'])
            len_set[str(file)] = len(fchange_norm.columns)
            #new label for columns
            fchange_norm.columns = ['Res '+str(x) for x in \
                    list(range(start, end))] # name the columns

            #figure out average wt values for set and linear transform
            cols = fchange_norm.columns
            wt_vals = []
            for row, col in zip(wt_subseq, cols):
                wt_vals.append(fchange_norm.loc[row, col])
            fchange_norm = fchange_norm - np.mean(wt_vals)
            rep1_set.append(fchange_norm)

        elif file in res_redo:
            set_ind = file.find('R') #R denotes the repeated set
            set_redo = file[:set_ind]
            start = list(samples[samples['Set'] == str(file)]['Start range'])[0]
            end = list(samples[samples['Set'] == str(file)]['End range'])[0]
            sites_ = list(samples[samples['Set'] == str(file)]['Sites'])[0]
            sites = ['Res '+ str(x) for x in sites_.split(',')]
            fchange = pd.read_csv(replicate_dir + str(file) \
                    + '_replicate'+str(rep)+cond_suffix, index_col = [0])
            fchange.columns = ['Res '+str(x) for x in list(range(start, end))]
            fchange = fchange[sites]
            #find WT residues for the set
            wt_subseq = [wt_full[int(ind)] for ind in sites_.split(',')]

            flat_list = np.array([item for sublist in fchange.values\
                    for item in sublist])
            cols = fchange.columns
            wt_vals = []
            for row, col in zip(wt_subseq, cols):
                wt_vals.append(fchange.loc[row, col])

            # Calculate scaling values for slotting in individual residues
            wt_mean = np.mean(wt_vals)
            stop_mean = np.mean(fchange.loc['*'])
            scale_factor = mean_stop[set_redo]/(stop_mean-wt_mean)
            fchange_norm = (fchange - wt_mean)*scale_factor
            set_list_res.append(fchange_norm)

    all_residues = pd.concat(rep1_set, axis = 1)
    all_res_redo = pd.concat(set_list_res, axis = 1)
    all_res_redo = all_res_redo.fillna('NaN')
    all_residues.update(all_res_redo)
    order = ['Res '+ str(x) for x in range(1, 307)]
    all_residues = all_residues[order]
    all_residues = all_residues.applymap(lambda x: x if not \
            isinstance(x, str) else np.nan)
    return(all_residues)

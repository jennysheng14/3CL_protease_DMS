import pandas as pd
import function_bio_rep
import replicates_no_syn
from Bio.Seq import Seq
import numpy as np
from pathlib import Path
import glob
import re

seq_3CL = ("TACAAAATG"
"TACAAAATGAGTGGTTTTAGAAAAATGGCATTCCCATCTGGTAAAGTTGAGGGTTGTATGGT"
"ACAAGTAACTTGTGGTACAACTACACTTAACGGTCTTTGGCTTGATGACGTAGTTTACTGTCCAAGACATGT"
"GATCTGCACCTCTGAAGACATGCTTAACCCTAATTATGAAGATTTACTCATTCGTAAGTCTAATCATAATTTC"
"TTGGTACAGGCTGGTAATGTTCAACTCAGGGTTATTGGACATTCTATGCAAAATTGTGTACTTAAGCTTAAGG"
"TTGATACAGCCAATCCTAAGACACCTAAGTATAAGTTTGTTCGCATTCAACCAGGACAGACTTTTTCAGTGT"
"TAGCTTGTTACAATGGTTCACCATCTGGTGTTTACCAATGTGCTATGAGGCCCAATTTCACTATTAAGGGTTC"
"ATTCCTTAATGGTTCATGTGGTAGTGTTGGTTTTAACATAGATTATGACTGTGTCTCTTTTTGTTACATGCAC"
"CATATGGAATTACCAACTGGAGTTCATGCTGGCACAGACTTAGAAGGTAACTTTTATGGACCTTTTGTTGACA"
"GGCAAACAGCACAAGCAGCTGGTACGGACACAACTATTACAGTTAATGTTTTAGCTTGGTTGTACGCTGCTGT"
"TATAAATGGAGACAGGTGGTTTCTCAATCGATTTACCACAACTCTTAATGACTTTAACCTTGTGGCTATGAAGT"
"ACAATTATGAACCTCTAACACAAGACCATGTTGACATACTAGGACCTCTTTCTGCTCAAACTGGAATTGCCGT"
"TTTAGATATGTGTGCTTCATTAAAAGAATTACTGCAAAATGGTATGAATGGACGTACCATATTGGGTAGTGCTT"
"TATTAGAAGATGAATTTACACCTTTTGATGTTGTTAGACAATGCTCAGGTGTTACTTTCCAATAA")
wt_full = ('MSGFRKMAFPSGKVEGCMVQVTCGTTTLNGLWLDDVVYCPRHVICT'
           'SEDMLNPNYEDLLIRKSNHNFLVQAGNVQLRVIGHSMQNCVLKLKV'
           'DTANPKTPKYKFVRIQPGQTFSVLACYNGSPSGVYQCAMRPNFTIK'
           'GSFLNGSCGSVGFNIDYDCVSFCYMHHMELPTGVHAGTDLEGNFYG'
           'PFVDRQTAQAAGTDTTITVNVLAWLYAAVINGDRWFLNRFTTTLND'
           'FNLVAMKYNYEPLTQDHVDILGPLSAQTGIAVLDMCASLKELLQNG'
           'MNGRTILGSALLEDEFTPFDVVRQCSGVTFQ')

amino_acid_list = ['*', 'A', 'C', 'D', 'E', 'F', 'G', 'H',
                   'I', 'K', 'L', 'M', 'N', 'P', 'Q', 'R',
                   'S', 'T', 'V', 'W', 'Y']
amino_acid_list.reverse()
COMPARISONS = ['Glu_Gal', 'Glu_Gc', 'Glu_Grl', 'Gal_Gc', 'Gal_Grl']
spreadsheet = "sample_spreadsheet_021521.csv"
samples = pd.read_csv(spreadsheet, comment = '#')
sets = [1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18, \
       19,20,21,'R1', '8R', '13R1', '14R', '13R2', '16R',\
       '9R1','9R2', '10R1', '10R2']
sets = set(list(samples['Set']))

def get_fastq_files(sample_spreadsheet):
    '''
    Get the names of all the fastq files required.
    '''
    samples = pd.read_csv(sample_spreadsheet, comment = '#')
    fastq_files = []
    for s in sets:
        for condition in ['Glu', 'Gal', 'Gc', 'Grl']:
            x = str(s)
            start = list(samples[samples['Set'] == x]['Start range'])[0]
            end = list(samples[samples['Set'] == x]['End range'])[0]
            # test to see if I can directly get this data from google bucket
            # files = 'jenny_yeast/'+ list(sample[sample['Set'] == x]['Folder'] + \
            #     sample[sample['Set'] == x][condition] + '_R1.fastq.gz')
            files = list(samples[samples['Set'] == x]['Folder'] + \
                    samples[samples['Set'] == x][condition] + '_R1.fastq.gz')
            fastq_files.append(files)
    return(fastq_files)

def replicate_mean(rep1, rep2, thresh1, thresh2):
    '''
    Return the mean value of coding between two biological replicates. If
    only a single replicate exists take that as the mean.
    '''
    rep1 = pd.read_csv(rep1, index_col = 0)
    rep2 = pd.read_csv(rep2, index_col = 0)
    rep1_copy = rep1[(rep1['count_x']>thresh1)&(rep1['count_x']>thresh2)]
    rep2_copy = rep2[(rep2['count_x']>thresh1)&(rep2['count_x']>thresh2)]
    merged = rep1_copy.merge(rep2_copy, on = \
        ['site_1', 'site_2', 'site_3'], how = 'outer')
    merged['mean'] = merged[['ratio_x', 'ratio_y']].mean(axis=1)
    merged['std'] = merged[['ratio_x', 'ratio_y']].std(axis = 1)
    merged = merged.set_index(merged['site_1'] + \
                     merged['site_2'] + \
                     merged['site_3'])
    keep = ['N' not in x for x in merged.index]
    merged = merged[keep]
    return(merged)

def amino_acids_vals(cond, wt_site):
    '''
    Return a dataframe of mean values for the comparison being made--means
    are averaged between the two replicates.
    '''
    cond_df = pd.read_csv(cond, index_col = 0)
    cond_df.drop(wt_site)
    wt_first_aa = Seq(wt_site[0:3]).translate()[0]
    wt_last_aa = Seq(wt_site[6:9]).translate()[0]
    ind = list(cond_df.index)
    ind = [x for x in ind if len(x)==9 and '\n' not in x]
    ind.sort(key = lambda x:(x[3:6], x[0:3], x[6:9]))#sort display order
    # keep only synonymous on flanking
    keep_ind = [x for x in ind if Seq(x[0:3]).translate()[0]==\
            wt_first_aa and Seq(x[6:9]).translate()[0]== wt_last_aa]
    sorted_diff = cond_df.reindex(keep_ind)
    codons = [str(Seq(x).translate()) for x in sorted_diff.index]
    sorted_diff['Translation'] = codons
    # Dataframe for amino acid data
    g = sorted_diff.groupby('Translation')
    aa_df = pd.DataFrame(g['mean'].apply(list)) #mean of all codings
    # Propogate errors
    # add in standard error of mean here in addn to propogated errors
#         aa_df['std'] = g.apply(lambda x: np.sqrt(sum(x**2)))['std']
    return(aa_df)

def amino_acid_means(cond, wt_site):
    '''
    Given dataframe of mean values generated by amino_acid_vals,
    return mean values.
    '''
    cond_df = pd.read_csv(cond, index_col = 0)
    cond_df.drop(wt_site)
    wt_first_aa = Seq(wt_site[0:3]).translate()[0]
    wt_last_aa = Seq(wt_site[6:9]).translate()[0]
    ind = list(cond_df.index)
    ind = [x for x in ind if len(x)==9 and '\n' not in x]
    ind.sort(key = lambda x:(x[3:6], x[0:3], x[6:9]))#sort display order
    # keep only synonymous on flanking
    keep_ind = [x for x in ind if Seq(x[0:3]).translate()[0]==\
    wt_first_aa and Seq(x[6:9]).translate()[0]== wt_last_aa]
    sorted_diff = cond_df.reindex(keep_ind)
    codons = [str(Seq(x).translate()) for x in sorted_diff.index]
    sorted_diff['Translation'] = codons
    # Dataframe for amino acid data
    g = sorted_diff.groupby('Translation')
    aa_df = pd.DataFrame(g.mean()['mean']) #mean of all codings
    aa_df['std'] = pd.DataFrame(g.std()['mean'])
    aa_df['len'] = pd.DataFrame(g.size())
    return(aa_df)

def properties(spreadsheet):
    '''Define properties wildcard.'''
    samples = pd.read_csv(spreadsheet, comment = '#')
    props = []
    for s in sets:
        x = str(s)
        start = list(samples[samples['Set'] == x]['Start range'])[0]
        end = list(samples[samples['Set'] == x]['End range'])[0]
        for ind in range(start, end):
            property = str(x)+'_residue'+str(ind)
            props.append(property)
    return(props)

def sets_and_residues(spreadsheet):
    '''
    Define which residues to take from which sets, especially for repeated
    residues.
    '''
    samples = pd.read_csv(spreadsheet, comment = '#')
    set_ = []
    res = []
    for s in sets:
        x = str(s)
        if 'R' in str(s) and str(s)!= 'R1':
            sites = list(samples[samples['Set'] == str(x)]['Sites'])[0]
            sites = [str(x) for x in sites.split(',')]
            for site in sites:
                set_.append(x)
                res.append(site)
    for s in sets:
        x = str(s)
        if 'R' not in str(s) or str(s) == 'R1':
            start = list(samples[samples['Set'] == x]['Start range'])[0]
            end = list(samples[samples['Set'] == x]['End range'])[0]
            for site in range(start, end):
                if str(site) not in res:
                    set_.append(x)
                    res.append(str(site))
    return(list(zip(set_, res)))

<<<<<<< HEAD
def transform_matrix(spreadsheet, raw_matrix, std_matrix):
    '''
    Transform each set so that WT fixed at 0 and stop codon is normalized to
    -1 in each set.

    Also transform the raw standard deviations by the same transformation
=======
def transform_matrix(spreadsheet, raw_matrix):
    '''
    Transform each set so that WT fixed at 0 and stop codon is normalized to
    -1 in each set.
>>>>>>> ba44b4aecae034e218ff16b7fea508545c9a38a5
    __________
    Input:
    raw_matrix: matrix of untransformed values
    spreadsheet: spreadsheet indexing sets and residues
    set21: set21--treated separately because of the C terminus
    '''
    raw_matrix = pd.read_csv(raw_matrix, index_col = 0)
<<<<<<< HEAD
    std_matrix = pd.read_csv(std_matrix, index_col = 0)
=======
>>>>>>> ba44b4aecae034e218ff16b7fea508545c9a38a5
    set_res = sets_and_residues(spreadsheet)
    set_res = pd.DataFrame(set_res, columns = ['set', 'residue'])
    sets = list(set(pd.DataFrame(set_res, columns = ['set', 'residue'])['set']))
    mean_stop = {}
    len_set = {}
    set_list = []
<<<<<<< HEAD
    std_list = []
=======
>>>>>>> ba44b4aecae034e218ff16b7fea508545c9a38a5
    for set_ in sets:
        residues = [str(x) for x in list(set_res[set_res['set']==\
                set_]['residue'])]
        if set_ != '21':
            fchange = raw_matrix[residues]
<<<<<<< HEAD
            fchange_std = std_matrix[residues]
=======
>>>>>>> ba44b4aecae034e218ff16b7fea508545c9a38a5
            wt_subseq = [wt_full[int(i)] for i in residues] #find WT residues for the set
            flat_list = np.array([item for sublist in fchange.values\
                for item in sublist])
            mean = flat_list[~np.isnan(flat_list)].mean() # mean of the set
            var = flat_list[~np.isnan(flat_list)].var() # variance of the set

            wt_vals = []
            for row, col in zip(wt_subseq, residues):
                wt_vals.append(fchange.loc[row, col])
            wt_mean = np.mean(wt_vals)
            fchange = fchange - wt_mean
            mean_stop[str(set_)] = np.mean(fchange.loc['*'])
            len_set[str(set_)] = len(fchange.columns)

            stop_mean = np.mean(fchange.loc['*'])
            scale_factor = -1/stop_mean
            fchange_norm = fchange*scale_factor
<<<<<<< HEAD
            norm_std = abs(fchange_std*scale_factor)
            set_list.append(fchange_norm)
            std_list.append(norm_std)
        elif set_ == '21':
            fchange = raw_matrix[residues]
            fchange_std = std_matrix[residues]
=======
            set_list.append(fchange_norm)
        elif set_ == '21':
            fchange = raw_matrix[residues]
>>>>>>> ba44b4aecae034e218ff16b7fea508545c9a38a5
            wt_subseq = [wt_full[int(i)] for i in residues]
            cols = fchange.columns[:2]
            wt_vals = []
            for row, col in zip(wt_subseq, cols):
                wt_vals.append(fchange.loc[row, col])
            wt_mean = np.mean(wt_vals)
            fchange = fchange - wt_mean
            # add to dict for mean stop
            mean_stop[str(set_)] = np.mean(fchange.loc['*'][:2])
            len_set[str(set_)] = 2
            stop_mean = np.mean(fchange.loc['*'][:2])
            scale_factor = -1/stop_mean
            fchange_norm = (fchange - wt_mean)*scale_factor
<<<<<<< HEAD
            norm_std = abs(fchange_std*scale_factor)
            std_list.append(norm_std)
            set_list.append(fchange_norm)
    all_residues = pd.concat(set_list, axis = 1)
    all_std = pd.concat(std_list, axis = 1)
    order = [str(x) for x in range(1, 307)]
    all_residues = all_residues[order]
    all_std = all_std[order]
    all_residues = all_residues.applymap(lambda x: x if not \
        isinstance(x, str) else np.nan)
    all_std = all_std.applymap(lambda x: x if not \
        isinstance(x, str) else np.nan)
    return(all_residues, all_std)
=======
            set_list.append(fchange_norm)
    all_residues = pd.concat(set_list, axis = 1)
    order = [str(x) for x in range(1, 307)]
    all_residues = all_residues[order]
    all_residues = all_residues.applymap(lambda x: x if not \
        isinstance(x, str) else np.nan)
    return(all_residues)
>>>>>>> ba44b4aecae034e218ff16b7fea508545c9a38a5

def sum_counts_nosyn(file, wt_site):
    '''
    file--file of count_matrix
    residue--residue number of file being processed
    '''
    counts = pd.read_csv(file, index_col = 0)
    wt_first_aa = Seq(wt_site[0:3]).translate()[0]
    wt_last_aa = Seq(wt_site[6:9]).translate()[0]
    ind = counts['site_1']+counts['site_2']+counts['site_3']
    counts['index']=ind
    ind = [x for x in ind if len(x)==9 and '\n' not in x]
    counts.set_index('index', inplace = True)
    # keep only synonymous on flanking
    keep_ind = [x for x in ind if Seq(x[0:3]).translate()[0]==\
            wt_first_aa and Seq(x[6:9]).translate()[0]== wt_last_aa]
    sorted_diff = counts.reindex(keep_ind)
    no_syn = sorted_diff.groupby('site_2').sum()['count']
    return no_syn

def amino_acid_nosyn(df1, df2):
    '''
    Takes dataframes of comparison data (both replicates) and reports
    ratio(foldchanges) of individual variants across both replicates.
    '''
    df1['Translation'] = [str(Seq(x).translate()) for x in df1.index]
    df2['Translation'] = [str(Seq(x).translate()) for x in df2.index]
    g1 = df1.groupby('Translation')
    g2 = df2.groupby('Translation')
    aa_df1 = pd.DataFrame(g1['ratio'].apply(list)) #mean of all codings
    aa_df2 = pd.DataFrame(g2['ratio'].apply(list))
    merged = aa_df1.merge(aa_df2, on = 'Translation')
    merged['all_ratios'] = merged['ratio_x']+ merged['ratio_y']
    merged['len'] = merged['all_ratios'].apply(lambda x: len(x))
    merged['mean'] = merged['all_ratios'].apply(lambda x: np.mean(x))
    return merged
<<<<<<< HEAD

# amalgamate_count_matrix():
#     '''
#     Amalgamates all of the individual count_matrices into a single one
#     for feeding into deseq2.
#     '''
=======
>>>>>>> ba44b4aecae034e218ff16b7fea508545c9a38a5

import pandas as pd
import function_bio_rep
from Bio.Seq import Seq
import numpy as np
from pathlib import Path
# from snakemake.remote.GS import RemoteProvider as GSRemoteProvider
# GS = GSRemoteProvider()

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

amino_acid_list = ['*', 'A', 'C', 'D', 'E', 'F', 'G', 'H',
                   'I', 'K', 'L', 'M', 'N', 'P', 'Q', 'R',
                   'S', 'T', 'V', 'W', 'Y']
COMPARISONS = ['Glu_Gal', 'Glu_Gc', 'Glu_Grl', 'Gal_Gc', 'Gal_Grl']
spreadsheet = "sample_spreadsheet_021521.csv"

def get_fastq_files(sample_spreadsheet):
    sample = pd.read_csv(sample_spreadsheet, comment = '#')
    fastq_files = []
    for s in sets:
        for condition in ['Glu', 'Gal', 'Gc', 'Grl']:
            x = str(s)
            start = list(sample[sample['Set'] == x]['Start range'])[0]
            end = list(sample[sample['Set'] == x]['End range'])[0]
            # test to see if I can directly get this data from google bucket
            # files = 'jenny_yeast/'+ list(sample[sample['Set'] == x]['Folder'] + \
            #     sample[sample['Set'] == x][condition] + '_R1.fastq.gz')
            files = list(sample[sample['Set'] == x]['Folder'] + \
                    sample[sample['Set'] == x][condition] + '_R1.fastq.gz')
            fastq_files.append(files)
    return(fastq_files)

def replicate_mean(rep1, rep2, thresh1, thresh2):
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
    cond_df = pd.read_csv(cond, index_col = 0)
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

amino_acid_list.reverse()
sets = [1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18, \
       19,20,21,'R1', '8R', '13R1', '14R', '13R2', '16R',\
       '9R1','9R2', '10R1', '10R2']

rule all:
    # IDS, glob_wildcards('amino_acid_dir_{id}.csv')
    input:
        # "comparison_matrices_created.csv",
        # dynamic('{cond1}_{cond2}_comparison/set{set}_rep{rep}_residue{res}.csv')
#         "amino_acid_matrices_created.csv",
        # expand('amino_acid_dir_{id}.csv', id = IDS),
        dynamic('amino_acid_dir_{stuff2}.csv'),
#         # dynamic('{cond1}_{cond2}_replicate/set{set}_residue{res}.csv'),
#         dynamic("{cond1}_{cond2}_amino_acid/set{set}_residue{res}.csv")

rule amino_acid_vals:
    input:
        'replicatedir_{stuff2}.csv'
    output:
        "amino_acid_dir_{stuff2}.csv"
    run:
        sample = pd.read_csv(spreadsheet, comment = '#')
        for comparison in [['Glu', 'Gal'], ['Glu', 'Gc'], ['Glu', 'Grl'],
                           ['Gal', 'Gc'], ['Gal', 'Grl']]:
            Path('amino_acid_dir_' + comparison[0] + '_' + comparison[1]).mkdir(parents=True, exist_ok=True)
            for s in sets:
                x = str(s)
                start = list(sample[sample['Set'] == x]['Start range'])[0]
                end = list(sample[sample['Set'] == x]['End range'])[0]
                for ind in range(start, end):
                    file = 'replicatedir_' + comparison[0] + '_' + comparison[1]+'/set' +\
                        str(x)+'_residue'+str(ind)+'.csv'
                    site = function_bio_rep.mutations(list(range(ind, ind+1)), list(range(ind, ind+1)))
                    wt = site.wildtype()
                    aa_vals = amino_acids_vals(file, wt[0])
                    name = 'amino_acid_dir_' + comparison[0] + '_' + comparison[1]+'/set' +\
                            str(x)+'_residue'+str(ind)+'.csv'
                    aa_vals.to_csv(name)

rule replicate:
    input:
        dynamic('comparisondir_{stuff}.csv')
    output:
        dynamic('replicatedir_{stuff2}.csv')
    run:
        sample = pd.read_csv(spreadsheet, comment = '#')
        output_files = []
        thresh_dict = {'Glu': 100, 'Gal': 30, 'Grl': 30, 'Gc': 30}
        for comparison in [['Glu', 'Gal'], ['Glu', 'Gc'], ['Glu', 'Grl'],
                           ['Gal', 'Gc'], ['Gal', 'Grl']]:
            Path('replicatedir_' + comparison[0] + '_' + comparison[1]).mkdir(parents=True, exist_ok=True)
            for s in sets:
                x = str(s)
                start = list(sample[sample['Set'] == x]['Start range'])[0]
                end = list(sample[sample['Set'] == x]['End range'])[0]
                for ind in range(start, end):
                    file1 = 'comparisondir_'+comparison[0]+'_'+comparison[1]+'/set'+ x +\
                        '_rep0_residue'+str(ind)+'.csv'
                    file2 = 'comparisondir_'+comparison[0]+'_'+comparison[1]+'/set'+ x +\
                        '_rep1_residue'+str(ind)+'.csv'
                    replicates = replicate_mean(file1, file2,
                        thresh_dict[comparison[0]],thresh_dict[comparison[1]])
                    name = 'replicatedir_'+comparison[0] + '_' + comparison[1]+'/set' +\
                        str(x)+'_residue'+str(ind)+'.csv'
                    replicates.to_csv(name)

rule comp:
    input:
        spreadsheet = "sample_spreadsheet_021521.csv",
        count_matrix_list = "count_matrices_created.csv"
        # count_matrix_files = list(pd.read_csv("count_matrices_created.csv", index_col = 0)['0'])
    output:
        dynamic('comparisondir_{stuff}.csv')
    run:
        sample = pd.read_csv(input.spreadsheet, comment = '#')
        output_files = []
        for s in sets:
            x = str(s)
            start = list(sample[sample['Set'] == x]['Start range'])[0]
            end = list(sample[sample['Set'] == x]['End range'])[0]
            for comparison in [['Glu', 'Gal'], ['Glu', 'Gc'], ['Glu', 'Grl'],
                               ['Gal', 'Gc'], ['Gal', 'Grl']]:
                Path('comparisondir_'+ comparison[0] + '_' + comparison[1]).mkdir(parents=True, exist_ok=True)
                for rep in [0, 1]:
                    for ind in range(start, end):
                        file1 = comparison[0]+'_count_matrices/set'+ x +\
                            '_rep_'+str(rep)+'residue'+str(ind)+'.csv'
                        file2 = comparison[1]+'_count_matrices/set'+ x +\
                            '_rep_'+str(rep)+'residue'+str(ind)+'.csv'
                        cond1_dfs = pd.read_csv(file1, index_col = 0)
                        cond2_dfs = pd.read_csv(file2, index_col = 0)
                        #dataframes summarizing each cond2 site
                        cond1_norm = cond1_dfs.copy()
                        cond2_norm = cond2_dfs.copy()
                        cond1_norm['proportion'] = cond1_dfs['count']/\
                            cond1_dfs['count'].iloc[-1]
                        cond2_norm['proportion'] = cond2_dfs['count']/\
                            cond2_dfs['count'].iloc[-1]

                        merged = cond1_norm.merge(cond2_norm, on = \
                        ['site_1', 'site_2', 'site_3'])
                            # take ratio between conditions
                        merged['ratio'] = merged['proportion_x'].apply(lambda x: np.log2(x))-\
                                merged['proportion_y'].apply(lambda x: np.log2(x))
                        name = 'comparisondir_' + comparison[0] + '_' + comparison[1]+'/set' +\
                            str(x)+'_rep'+str(rep)+'_residue'+str(ind)+'.csv'
                        merged.to_csv(name)

rule count_matrix:
    input:
        fastq_files = get_fastq_files("sample_spreadsheet_021521.csv"),
        spreadsheet = "sample_spreadsheet_021521.csv"
    output:
        summary = "count_matrices_created.csv"
    run:
        sample = pd.read_csv(input.spreadsheet, comment = '#')
        output_files = []
        threshold = 1
        for s in sets:
            for condition in ['Glu', 'Gal', 'Gc', 'Grl']:
                x = str(s)
                start = list(sample[sample['Set'] == x]['Start range'])[0]
                end = list(sample[sample['Set'] == x]['End range'])[0]
                # test to see if I can directly get this data from google bucket
                # files = 'jenny_yeast/'+ list(sample[sample['Set'] == x]['Folder'] + \
                #     sample[sample['Set'] == x][condition] + '_R1.fastq.gz')
                files = list(sample[sample['Set'] == x]['Folder'] + \
                        sample[sample['Set'] == x][condition] + '_R1.fastq.gz')
                sequence = list(sample[sample['Set'] == x]['Sequence'])[0]
                position = list(sample[sample['Set'] == x]['Position'])[0]
                sites = function_bio_rep.mutations(list(range(start, end)),\
                        list(range(start, end)))
                Path(condition+'_count_matrices').mkdir(parents=True, exist_ok=True)
                print(sites.sites, sites.all_muts, position, threshold, sequence)
                for rep in [0, 1]:
                    count_mat = sites.count_matrix(files[rep], \
                                                sequence, position, threshold)
                    for ind, y in enumerate(count_mat):
                        name = condition+'_count_matrices/set'+ str(x)+\
                        '_rep_'+str(rep)+'residue'+str(start+ind)+'.csv'
                        y.to_csv(name)
                    output_files.append(name)
        output_df = pd.DataFrame(output_files)
        output_df.to_csv(output.summary)

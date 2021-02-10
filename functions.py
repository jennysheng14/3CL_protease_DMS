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
from Bio.Seq import Seq

class mutations:

    '''Characterizes mutation frequency at specific sites of interest within
    the 3CL protease. If working with a different protein, sequence in the
    seq_3CL variable should be changed accordingly. '''

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

    def __init__(self, sites, all_muts):
        self.sites = sites #codon site where SDM is expected. Enter as list
        self.wt = []
        check = all_muts
        check.insert(0, all_muts[0]-1)
        self.all_muts = check #should include one upstream and downstream of

        flank = check.copy()
        flank.insert(len(flank), all_muts[-1]+1)
        self.flank_muts = flank
        # actual mutated sites
        # position/index of site of interest in amongst all mutations
        self.pos = [i for i, x in enumerate(all_muts) if x in sites]

    def read_file_set(self, file, sequence, seq_position):

        '''
        Reads the lines in a single file and returns
        the amino acid in the desired positions.

        _______
        Input:
        file--fastq.gz file
        sequence--reference sequence for which the position is known. Must be
        string
        seq_position--the codon number encoded by the first codon in
        sequence + 2

        Output: Per read list of codons for all mutated positions
        '''

        N_list = []

        with gzip.open(file, 'r') as f:
            for line in itertools.islice(f, 1, None, 4):
                dna_codings = []
                line_dc = line.decode("utf-8")
                index = line_dc.find(sequence)
                for codon in self.all_muts:
                    dna_codings.append(line_dc[index+(codon-seq_position)*3:\
                    index+(codon-seq_position)*3+3])
                    #grab the three bases before and three after by default
                    #because of possible synonymous codings
                N_list.append(dna_codings)
        return(N_list)


    def wildtype(self):
        '''
        Takes the given codon position and finds wildtype sequence at that codon

        Input:
        codon_positions: list of sites where mutations are expected.
        optional argument: list of sites where additional wt sequences are
        desired.
        Returns:
        Wildtype codons at the codon position of interest
        '''
        wt_all = []
        wt_optional = []
        for x in self.sites:
            wt = self.seq_3CL[3*(x+4): 3*(x+7)]
            wt_all.append(wt)
        return(wt_all)

    def wt_codon(self):
        '''
        Takes the given codon position and finds wildtype sequence at that codon

        Input:
        codon_positions: list of sites where mutations are expected.
        optional argument: list of sites where additional wt sequences are
        desired.
        Returns:
        Wildtype codons at the codon position of interest
        '''
        wt_all = []
        wt_optional = []
        for x in self.flank_muts:
            wt = self.seq_3CL[3*(x+5): 3*(x+6)]
            wt = Seq(wt).translate()
            wt_all.append(wt)
        return(wt_all)

    def wildtype5(self):
        '''
        Takes the given codon position and finds wildtype sequence at that codon
        wildtype5 finds the wildtype sequence of 5 codon region. Two upstream
        and two downstream of site of interest. Useful for aligning border
        mutations to correct library.

        Input:
        codon_positions: list of sites where mutations are expected.
        optional argument: list of sites where additional wt sequences are
        desired.
        Returns:
        Wildtype codon at the codon position of interest
        '''
        wt_all = []
        wt_optional = []
        for x in self.sites:
            wt = self.seq_3CL[3*(x+3): 3*(x+8)]
            wt_all.append(wt)
        return(wt_all)

    def efficiency(self, file, sequence, seq_position):
        '''
        Takes the array of lists where the items in the lists are the -1 0 +1
        codons around the list in codon_positions. Note that dimensions must
        match. This does not take into account other mutations outside sites
        of interest.

        Returns the percentage of codons read at that site that are
        wildtype.
        '''
        N_lst = self.read_file_set(file, sequence, seq_position)
        wt_seq = self.wildtype()
        efficiency = []

        for ind, seq in list(enumerate(wt_seq)):
            single = [x[ind] for x in N_lst]
            num = Counter(single)[seq]
            denom = sum(Counter(single).values())
            frac = num/denom
            efficiency.append(frac)
        return(efficiency)

    def single_mutation(self, file, sequence, seq_position):

        '''
        Takes the array of lists where the items in the lists are the -1 0 +1
        codons around the position indicated in the list in codon_positions.
        Note that dimensions must match.

        all_muts are the positions of the mutations altered in the pool.

        Returns a matrix of counts with all detected unique combinations around
        codon of interest--including detected synonyous mutations.
        '''
        N_lst = self.read_file_set(file, sequence, seq_position)
        tot = len(N_lst)
        # N_join = np.asarray([''.join(x) for x in N_lst])
        codons = []
        wt = self.wildtype()
        wt5 = self.wildtype5()
        codon_dfs = []
        for codon in self.sites:
            codons.append([codon-1, codon, codon+1])
        #
        # #find all the wildtypes first for all of the positions

        N_lst = [x for x in N_lst if ''.join(x)!=\
        self.seq_3CL[3*(min(self.all_muts)+5) : 3*(max(self.all_muts)+5)+3]]
        wt_seq = tot - len(N_lst)

                # N_lst.remove(x) #remove the wildtypes for downstream. Find faster way
                #vectorize this part
        for ind, site in list(enumerate(self.sites)): #for each site of interest
            #full list will be sequences of 9bp with -1,0,+1
            sequence = []
            #if not first, second last, or last site
            if ind!=0 and ind!=len(self.sites)-1 and ind!=len(self.sites)-2:
                for x in N_lst: #for each read
                    # #for each site of interest check 5 amino acids
                    # to include one upstream and one downstream to
                    # make sure that no mutations outside of frame
                    flank = [x[self.pos[ind]-2], x[self.pos[ind]+2]]
                    coding = [x[self.pos[ind]-1], x[self.pos[ind]], \
                    x[self.pos[ind]+1]]

                # check if three codons are wiltype. If not add to list
                # Also check that flanking codons aren't mutated so we're not
                # looking at residuals of another library
                    if ''.join(coding) != wt[ind] and flank[0]==wt5[ind][0:3] and \
                    flank[1]==wt5[ind][-3:] and coding[1][2]!='A' and \
                    coding[1][2]!='C':
                    #add condition that we throw out non-NNK
                        sequence.append(coding)
            elif ind==0:
                for x in N_lst:
                    coding = [x[self.pos[ind]-1], x[self.pos[ind]], \
                    x[self.pos[ind]+1]]
                    flank = x[self.pos[ind]+2]
                    if ''.join(coding) != wt[ind] and len(''.join(coding))==9:
                        if coding[1][2]!='A' and \
                        coding[1][2]!='C':
                    #add condition that we throw out non-NNK
                            sequence.append(coding)
            elif ind==len(codons)-1:
                for x in N_lst:
                    coding = [x[self.pos[ind]-1], x[self.pos[ind]]]
                    flank = [x[self.pos[ind]-2]]
                    if ''.join(coding) != wt[ind][:6] and \
                    len(''.join(coding))==6:
                        if coding[1][2]!='A' and \
                        coding[1][2]!='C' and flank[0]==wt5[ind][0:3]:
                    #add condition that we throw out non-NNK
                            sequence.append(coding)
            elif ind==len(self.sites)-2:
                for x in N_lst:
                    coding = [x[self.pos[ind]-1], x[self.pos[ind]], \
                    x[self.pos[ind]+1]]
                    flank = [x[self.pos[ind]-2]]
                    if ''.join(coding) != wt[ind] and \
                    len(''.join(coding))==9:
                        if coding[1][2]!='A' and \
                        coding[1][2]!='C' and flank[0]==wt5[ind][0:3]:
                    #add condition that we throw out non-NNK
                            sequence.append(coding)

            codon_df = pd.DataFrame(sequence) #dataframe of three positions.
            codon_dfs.append(codon_df)
        return(codon_dfs, wt_seq)


    def count_matrix(self, file, sequence, seq_position):

        '''
        Builds the count matrices at the indicated positions for unique
        combinations of mutation at SDM site and synonymous sites.
        '''

        codon_df, wt_count = self.single_mutation(file, sequence, seq_position)

        count = [] # Count the occurences of each codon combo. List of
        # dataframes
        threshold = 1 # min number of counts needed to be included
        for ind, df in enumerate(codon_df):
            count_list = map(tuple, df.to_numpy().tolist())
            count_dict = Counter(count_list) #count of unique combos
            dict_keys = list(count_dict.keys())
            dict_keys.sort(key = lambda x:x[1]) #sort list based on second
            # codon: SDM codon
            index_map = {v: i for i, v in enumerate(dict_keys)}
            sorted_count = sorted(count_dict.items(), key=lambda pair: \
            index_map[pair[0]])

            if ind!=len(codon_df)-1:
                site = pd.DataFrame([x[0] for x in sorted_count], columns = \
                             ['site_1', 'site_2', 'site_3'])
                counts = pd.DataFrame([x[1] for x in sorted_count], columns = \
                ['count'])

                count_df = site.join(counts)
                #set threshold for min number of 50 counts
                count_df = count_df[count_df['count']>threshold] #set threshold
                count_df = count_df.reset_index(drop = True)
                count_df.loc[len(count_df)] =\
                [self.wildtype()[ind][0:3], self.wildtype()[ind][3:6],\
                 self.wildtype()[ind][6:], wt_count]

                count.append(count_df)
            else: #for last index of the set
                site = pd.DataFrame([x[0] for x in sorted_count], columns = \
                             ['site_1', 'site_2'])
                counts = pd.DataFrame([x[1] for x in sorted_count], columns = \
                ['count'])
                count_df = site.join(counts)
                #set threshold for min number of 50 counts
                count_df = count_df[count_df['count']>threshold] #set threshold
                count_df = count_df.reset_index(drop = True)
                count_df.loc[len(count_df)] = \
                [self.wildtype()[ind][0:3], self.wildtype()[ind][3:6], wt_count]
                count_df['site_3'] = pd.Series([self.wildtype()[ind][6:]]\
                *len(count_df))
                count.append(count_df)
        return(count)


    def replicate(self, files, sequence, seq_position):
        '''
        Reads multiple files which are all replicates and
        returns a mean and standard of the replicates and normalizes to the
        wildtype in the experiment.

        files--must be a list.
        '''
        count_amalgamate = [] #list of dataframes for each file with each
        # position in each sublist
        for ind, file in list(enumerate(files)):
            #ind indexes each replicate
            count_df = self.count_matrix(file, sequence, seq_position)
            count_per_file = [] # list of dataframes at each site for each file
            for x in range(len(count_df)): # x indexes each site of interest
                subcount = count_df[x].set_index(count_df[x]['site_1'] + \
                                 count_df[x]['site_2'] + \
                                 count_df[x]['site_3'])
                subcount = subcount.rename(columns={'count':'count' + str(ind)})
                subcount = subcount[~subcount.site_1.str.contains('\n')]
                subcount = subcount[~subcount.site_2.str.contains('\n')]
                subcount = subcount[~subcount.site_3.str.contains('\n')]
                subcount["site1_aa"] = subcount['site_1'].map(lambda name: \
                Seq(name).translate())
                subcount["site3_aa"] = subcount['site_3'].map(lambda name: \
                Seq(name).translate())

                subcount = subcount[subcount['site1_aa'] == self.wt_codon()[x]]
                subcount = subcount[subcount['site3_aa'] == self.wt_codon()[x+2]]
                subcount = subcount.drop(['site_1', 'site_2', 'site_3', \
                'site1_aa', 'site3_aa'], axis=1)
                norm_count = subcount.copy()
                norm_count = norm_count/norm_count.sum()
                # norm_count = norm_count[~norm_count.index.duplicated(keep='last')]
                # count number indexes the replicate
            # norm_count now needs to be normalized to the wildtype in that
            # file. Errors also need to be recorded.
# repeat of wildtype from sequencing error
                norm_wt = norm_count.div(norm_count.iloc[-1])
                count_per_file.append(norm_wt)
            count_amalgamate.append(count_per_file)

        # Find mean and sd.
        group_site = list(map(list, zip(*count_amalgamate)))
        # Invert list of lists so it's grouped by site
        summary = [] # each element in list is dataframe summarizing each site
        for site in group_site:
            grouped = pd.concat(site, axis = 1)
            summ = pd.DataFrame()
            summ['mean'] = grouped.mean(axis = 1)
            summ['std'] = grouped.std(axis = 1)
            summ['len'] = len(grouped.columns)
            summary.append(summ.dropna())
        return(summary) # each element in list is df summarizing each site

    def comparison(self, cond1, cond2, sequence, position):
        '''
        Compares expression values between two conditions.
        ___________
        Input:
        cond1-- list of files for cond1.
        cond2-- list of files for cond2.
        sequence-- reference sequence
        position-- codon position of reference sequence
        Output:
        Log(mut_cond1/wt_cond1)-Log(mut_cond2/wt_cond2)
        '''
        cond1_dfs = self.replicate(cond1, sequence, position)
        #dataframes summarizing each cond1 site
        cond2_dfs = self.replicate(cond2, sequence, position)
        #dataframes summarizing each cond2 site

        cond1_log = []
        cond1_error = []
        cond2_log = []
        cond2_error = []
        diff =[]
        err = []
        wt = self.wildtype()

        for ind, x in list(enumerate(cond1_dfs)):
            wt_seq = wt[ind]
            log = np.log2(x['mean']) #log2
            cond1_log.append(log)
            #variance of samples
            cond1_error.append(x['std']**2/(np.log(2)*x['mean']*x['len'])**2)
            #assume measurements are normal
            #assumption may not be entirely correct but easiest to deal w/rn

        for ind, x in list(enumerate(cond2_dfs)):
            wt_seq = wt[ind]
            log = np.log2(x['mean'])
            cond2_log.append(log)
            #variance--makes sum easier later
            cond2_error.append(x['std']**2/(np.log(2)*x['mean']*x['len'])**2)

        for index, (x, y, a, b) in enumerate(zip(cond1_log, cond2_log, \
        cond1_error, cond2_error)):
            difference = pd.DataFrame(x-y)
            #standard deviation of log ratio
            error =  pd.DataFrame(np.sqrt(a+b),columns = ['std'])
            wt_site = wt[index]
            # translate first and last amino acid to ensure synonymous flanking
            wt_first_aa = Seq(wt_site[0:3]).translate()[0]
            wt_last_aa = Seq(wt_site[6:9]).translate()[0]
            # sort dataframe by middle three base pairs, then first 3, last 3
            ind = list(difference.index)
            ind = [x for x in ind if len(x)==9 and '\n' not in x]
            ind.sort(key = lambda x:(x[3:6], x[0:3], x[6:9]))#sort display order
            keep_ind = [x for x in ind if Seq(x[0:3]).translate()[0]==\
            wt_first_aa and Seq(x[6:9]).translate()[0]== wt_last_aa]
            sorted_diff = difference.reindex(keep_ind)
            sorted_error = error.reindex(keep_ind)
            sorted_diff['std'] = sorted_error['std']
            # keep only synonymous codons in the pre- and post- positions
            diff.append(sorted_diff.dropna())
            # err.append(sorted_error)
        return (diff) # returns list of dataframes comparing sites

    def amino_acids(self, cond1, cond2, sequence, position):
        '''
        Compares expression values between two conditions based on amino acids.
        Processed from the codons from comparison function.
        ___________
        Input:
        cond1-- list of files for cond1.
        cond2-- list of files for cond2.
        sequence-- reference sequence
        position-- codon position of reference sequence
        Output:
        Log(mut_cond1/wt_cond1)-Log(mut_cond2/wt_cond2)
        '''

        seq_df = self.comparison(cond1, cond2, sequence, position)

        aa_amal = []
        for site in seq_df:
            codons = [str(Seq(x).translate()) for x in site.index]
            site['Translation'] = codons
            # Dataframe for amino acid data
            g = site.groupby('Translation')
            aa_df = pd.DataFrame(g.mean()['mean']) #mean of all codings
            aa_df['std_err_of_mean'] = pd.DataFrame(g.std()['mean'])**2
            # Propogate errors
            site['square'] = site['std']**2
            aa_df['sum_square'] = site.groupby('Translation').apply(lambda x: \
            x['square'].sum()/len(x)**2)
            # add in standard error of mean here in addn to propogated errors
            aa_df['std'] = np.sqrt(aa_df['sum_square']+aa_df['std_err_of_mean'])
    #         aa_df['std'] = g.apply(lambda x: np.sqrt(sum(x**2)))['std']
            aa_amal.append(aa_df)
        return(aa_amal)

    def mutation_matrix(self, cond1, cond2, sequence, position):
        dfs = amino_acids(cond1, cond2, sequence, position)
        aa_amalg = pd.DataFrame()
        for ind, df in enumerate(dfs):
            aa = [x[1] for x in df.index]
            mean = list(df['mean'])
            aa_df = pd.DataFrame({'Amino Acid': aa, 'Mean'+str(ind): mean})
            aa_df = aa_df.set_index('Amino Acid')
            aa_amalg = pd.concat([aa_amalg, aa_df], axis = 1, join = 'outer')
        return(aa_amalg)

def amino_acid_bar(df, show = True, save = False, *arg):
    '''
    Takes a single dataframe from list of dataframes
    from amino_acids attribute and plots it.
    ____________
    Input:
    df: dataframe
    show: whether to show the plot
    save: whethter to save the plot
    arg: filename to save the plot under
    '''

    fig = px.bar(df, y = 'mean', error_y='std')
    fig.update_layout(hovermode="x")
    if show == True:
        fig.show()
    if save == True:
        plotly.offline.plot(fig, filename = arg[0]+'.html')
def mutation_matrix(dfs):
    '''
    Takes the mutations from the amino_acids attribute
    (given in a dataframe with respective errors) and turns
    it into a panda dataframe (which can be turned into a matrix).
    '''
    aa_amalg = pd.DataFrame()
    for ind, df in enumerate(dfs):
        aa = [x[1] for x in df.index]
        mean = list(df['mean'])
        aa_df = pd.DataFrame({'Amino Acid': aa, 'Mean'+str(ind): mean})
        aa_df = aa_df.set_index('Amino Acid')
        aa_amalg = pd.concat([aa_amalg, aa_df], axis = 1, join = 'outer')
    return(aa_amalg)

def make_heatmap(df, x_, y_, wt_, show = True, save = False, **kwarg):
    '''
    Given a list of dataframe generated from
    amino_acid attribute, return an interactive plotly figure.
    _________________
    Input:
    df: dataframe of amino acids at each site.
    x_: x-axis labels
    y_: y-axis labels (list in inverted order)
    wt: a list of the wildtype residues at each of the amino acid positions

    Output: Plotly figure saved at location specified at filename
    '''
    fig = go.Figure(data = go.Heatmap(z = mutation_matrix(df).iloc[::-1],
                                 x = x_, y = y_))
    #Add marker to denote wildtype
    wt = wt_
    fig.add_scatter(y=wt, x = list(range(140, 150)), mode="markers",
                    marker=dict(size=4, color="White"),
                    name="wt")

    fig.update_layout({
        'plot_bgcolor': 'rgba(0, 0, 0, 0)',
        'paper_bgcolor': 'rgba(0, 0, 0, 0)',
        })

    if show == True:
        fig.show()
    if save == True:
        plotly.offline.plot(fig, filename = kwarg['name']+'.html')

def activity_km(dfs):
    '''
    With the comparison list of dataframes, predict
    protease activity categorically as on or off using
    K_means.

    Returns dataframe with a column that tells us whether ratio value is low
    or high.
    '''
    activity_df = []
    for ind, df in enumerate(dfs):
        reshaped = df['mean'].to_numpy().reshape(-1,1)
        kmeans = KMeans(init='k-means++', n_clusters=2, n_init=10)
        kmeans.fit(reshaped)
        labels = kmeans.labels_
        # the one closest to zero reflects the behavior of the wildtype
        centers = kmeans.cluster_centers_
        df2 = df.copy()
        if centers[0][0] < centers[1][0]:
            # 0 label is low relative change in abundance
            map_dict =  {0: 'low', 1: 'high'}
            df2['rel_proportion'] = [map_dict[i] for i in labels]
        else:
            map_dict =  {0: 'high', 1: 'low'}
            df2['rel_proportion'] = [map_dict[i] for i in labels]

        activity_df.append(df2)
    return(activity_df)

def drug_comparison(df_drug, df_gal):

    '''
    This function filters off of the conditions where the protease is still active using the
    activity function and compares it to the activity of the variant under drug (compared to gal)
    to report which residue mutations cause a drug escapee.
    ___________
    Input: Two activity dataframes. One for gal/glu and another for drug/gal.
    '''

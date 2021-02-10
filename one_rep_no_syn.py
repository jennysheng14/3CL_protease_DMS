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

                # N_lst.remove(x) #remove the wildtypes for downstream.
                # Find faster way to vectorize this part...
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
                    if ''.join(coding) != wt[ind] and flank[0]==wt5[ind][0:3] \
                    and \
                    flank[1]==wt5[ind][-3:] and coding[1][2]!='A' and \
                    coding[1][2]!='C':
                    #add condition that we throw out non-NNK
                        sequence.append(coding)
            elif ind==0:
                for x in N_lst:
                    coding = [x[self.pos[ind]-1], x[self.pos[ind]], \
                    x[self.pos[ind]+1]]
                    flank = x[self.pos[ind]+2]
                    if ''.join(coding) != wt[ind] and \
                    len(''.join(coding))==9:
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


    def count_matrix(self, file, sequence, seq_position, thresh):

        '''
        Builds the count matrices at the indicated positions for unique
        combinations of mutation at SDM site and synonymous sites.
        thresh--minimum number of counts of a barcode for considering it.
        '''

        codon_df, wt_count = self.single_mutation(file, sequence, seq_position)

        count = [] # Count the occurences of each codon combo. List of
        # dataframes
        # threshold = 1 # min number of counts needed to be included
        for ind, df in enumerate(codon_df):
            count_list = map(tuple, df.to_numpy().tolist())
            # dispose of information from flanking sites
            middle = [x[1] for x in count_list]
            count_dict = Counter(middle) #count of unique combos
            dict_keys = list(count_dict.keys())
            dict_keys.sort(key = lambda x:x[1]) #sort list based on second
            # codon: SDM codon
            index_map = {v: i for i, v in enumerate(dict_keys)}
            sorted_count = sorted(count_dict.items(), key=lambda pair: \
            index_map[pair[0]])

            if ind!=len(codon_df)-1:
                site = pd.DataFrame([x[0] for x in sorted_count], columns = \
                             ['site_2'])
                counts = pd.DataFrame([x[1] for x in sorted_count], columns = \
                ['count'])

                count_df = site.join(counts)
                #set threshold for min number of 50 counts
                count_df = count_df[count_df['count']>thresh] #set threshold
                count_df = count_df.reset_index(drop = True)
                count_df.loc[len(count_df)] = \
                    [self.wildtype()[ind][3:6],wt_count]

                count.append(count_df)
            else: #for last index of the set
                site = pd.DataFrame([x[0] for x in sorted_count], columns = \
                             ['site_2'])
                counts = pd.DataFrame([x[1] for x in sorted_count], columns = \
                ['count'])
                count_df = site.join(counts)
                #set threshold for min number of 50 counts
                count_df = count_df[count_df['count']>thresh] #set threshold
                count_df = count_df.reset_index(drop = True)
                count_df.loc[len(count_df)] = \
                [self.wildtype()[ind][3:6], wt_count]
                count.append(count_df)
        return(count)

    def comparison(self, cond1, cond2, sequence, position, thresh1, thresh2):
        '''
        Compares expression values between two conditions of the same biological
        replicate.
        ___________
        Input:
        cond1--file for cond1.
        cond2--file for cond2.
        sequence-- reference sequence
        position-- codon position of reference sequence
        thresh--threshold in count_matrix to count the barcode

        Output:
        Log2(mut_cond1/wt_cond1)-Log2(mut_cond2/wt_cond2)
        '''
        cond1_dfs = self.count_matrix(cond1, sequence, position, thresh1)
        # dataframes summarizing each cond1 site
        # normalize dataframes
        cond2_dfs = self.count_matrix(cond2, sequence, position, thresh2)
        # dataframes summarizing each cond2 site
        ratio = []
        for x in range(len(cond1_dfs)):
            # normalize to wildtype coding
            cond1_norm = cond1_dfs[x].copy()
            cond2_norm = cond2_dfs[x].copy()
            cond1_norm['count'] = cond1_dfs[x]['count']/\
                cond1_dfs[x]['count'].iloc[-1]
            cond2_norm['count'] = cond2_dfs[x]['count']/\
                cond2_dfs[x]['count'].iloc[-1]

            merged = cond1_norm.merge(cond2_norm, on = \
            ['site_2'])
            # take ratio between conditions
            merged['ratio'] = np.log2(merged['count_x'])-\
            np.log2(merged['count_y'])
            ratio.append(merged)
        return(ratio)

    def replicate(self, cond1, cond2, sequence, position, thresh1, thresh2):
        '''
        Given replicates of conditions, return dataframe for each codon in set
        with the foldchange of log2(cond1/cond2) each normalized to the wildtype
        coding.
        '''
        rep1 = self.comparison(cond1[0], cond2[0], sequence, position, thresh1,\
            thresh2)
        all_residues = []
        for x in range(len(rep1)): #for each residue
            rep1_copy = rep1[x].drop(['count_x','count_y'], axis = 1)
            rep1_copy = rep1_copy.set_index(rep1_copy['site_2'])
            keep = ['N' not in x for x in rep1_copy.index]
            rep1_copy = rep1_copy[keep]
            all_residues.append(rep1_copy)
        return(all_residues)

    def amino_acids(self, cond1, cond2, sequence, position, thresh1, thresh2):
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

        seq_df = self.replicate(cond1, cond2, sequence, position, thresh1, \
            thresh2)
        wt = self.wildtype()

        aa_amal = []
        for index, site in enumerate(seq_df):
            wt_site = wt[index]
            wt_first_aa = Seq(wt_site[0:3]).translate()[0]
            wt_last_aa = Seq(wt_site[6:9]).translate()[0]
            ind = list(site.index)
            ind = [x for x in ind if '\n' in x]
            sorted_diff = site.drop(ind)

            codons = [str(Seq(x).translate()) for x in sorted_diff.index]
            sorted_diff['Translation'] = codons
            # Dataframe for amino acid data
            g = sorted_diff.groupby('Translation')
            aa_df = pd.DataFrame(g.mean()['ratio']) #mean of all codings
            aa_df['std'] = pd.DataFrame(g.std()['ratio'])
            # Propogate errors
            # add in standard error of mean here in addn to propogated errors
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

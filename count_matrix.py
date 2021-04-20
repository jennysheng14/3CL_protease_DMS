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
                count_df = count_df[count_df['count']>thresh] #set threshold
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
                count_df = count_df[count_df['count']>thresh] #set threshold
                count_df = count_df.reset_index(drop = True)
                count_df.loc[len(count_df)] = \
                [self.wildtype()[ind][0:3], self.wildtype()[ind][3:6], wt_count]
                count_df['site_3'] = pd.Series([self.wildtype()[ind][6:]]\
                *len(count_df))
                count.append(count_df)
        return(count)

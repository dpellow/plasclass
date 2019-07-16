###
# Define the classifier class and provide a set of functions to enable classification
###

import numpy as np
import os
from sklearn.linear_model import LogisticRegression
from sklearn.preprocessing import StandardScaler
import itertools
from joblib import load

import sys
sys.path.insert(0, './lib')
import utils


import time

class classifier():
    def __init__(self):
        self._scales = [1000,10000,100000,500000]
        self._ks = [3,4,5,6,7]
        self._compute_kmer_inds()
        self._load_classifiers()


    def classify(self,seq):
        '''Classify the sequence(s), return the probability of the sequence(s) being a plasmid.
        Assumes seq is either an individual string or a list of strings
        Returns either an individual plasmid probability for seq or a list of
        plasmid probabilities for each sequence in seq
        '''
        if isinstance(seq, basestring): # single sequence
            print "Counting k-mers for sequence of length {}".format(len(seq))
            kmer_counts = utils.count_kmers(seq, self._ks)
            kmer_freqs = utils.counts2freqs(kmer_counts, self._ks)
            scale = self._get_scale(len(seq))
            standardized_freqs = self._standardize(kmer_freqs, scale)
            print "Classifying"
            return self.classifiers[scale]['clf'].predict_proba(standardized_freqs)[0,1]

        elif isinstance(seq, list): # list of sequences
            start_count = time.time()
            print "Counting k-mers for {} sequences".format(len(seq))
            all_kmer_counts = []
            scales = []
            for s in seq:
                kmer_counts = utils.count_kmers(s, self._ks, self._kmer_inds, self._kmer_count_lens)
                all_kmer_counts.append(kmer_counts)
                scales.append(self._get_scale(len(s)))
            end_count = time.time()

            print "{} s to count".format(end_count-start_count)
            kmer_freqs = utils.counts2freqs(all_kmer_counts,self._ks)
            end_freqs = time.time()
            print "{} s to convert freqs".format(end_freqs-end_count)
            # separate by length scale
            partitioned_classifications = {}
            print "Partitioning by length scale"
            for s in self._scales:
                # select the rows of kmer_freqs that are of this scale
                scale_freqs = kmer_freqs[[x == s for x in scales],:]
                if len(scale_freqs) <= 0: continue
                standardized_freqs = self._standardize(scale_freqs, s)
                print "Classifying sequences of length scale {}".format(s)
                partitioned_classifications[s] = self.classifiers[s]['clf'].predict_proba(standardized_freqs)[:,1]
            # recollate the results:
            results = []
            scale_inds = {s:0 for s in self._scales}
            for s in scales:
                results.append(partitioned_classifications[s][scale_inds[s]])
                scale_inds[s] += 1
            return np.array(results)

        else:
            raise TypeError('Can only classify strings or lists of strings')


    def _load_classifiers(self):
        ''' Load the multi-scale classifiers and scalers
        '''
        data_path = os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))),'data')

        self.classifiers = {}
        for i in self._scales:
            print "Loading classifier " + str(i)
            self.classifiers[i] = {'clf': load(os.path.join(data_path,'m'+str(i))), 'scaler': load(os.path.join(data_path,'s'+str(i)))}


    def _get_scale(self, length):
        ''' Choose which length scale to use for the sequence
        '''
        #TODO: binary search to make this more efficient
        if length <= self._scales[0]: return self._scales[0]
        for i,l in enumerate(self._scales[:-1]):
            if length <= float(l + self._scales[i+1])/2.0:
                return l
        return self._scales[-1]

    def _standardize(self, freqs, scale):
        ''' Use sklearn's standard scaler to standardize
        Choose the appropriate scaler based on sequence length
        '''
        return self.classifiers[scale]['scaler'].transform(freqs)

    def _compute_kmer_inds(self):
        ''' Compute the indeces of each canonical kmer in the kmer count vectors
        '''

        self._kmer_inds = {k: {} for k in self._ks}
        self._kmer_count_lens = {k: 0 for k in self._ks}

        alphabet = 'ACGT'
        for k in self._ks:
            all_kmers = [''.join(kmer) for kmer in itertools.product(alphabet,repeat=k)]
            all_kmers.sort()
            ind = 0
            for kmer in all_kmers:
                bit_mer = utils.mer2bits(kmer)
                rc_bit_mer = utils.mer2bits(utils.get_rc(kmer))
                if rc_bit_mer in self._kmer_inds[k]:
                    self._kmer_inds[k][bit_mer] = self._kmer_inds[k][rc_bit_mer]
                else:
                    self._kmer_inds[k][bit_mer] = ind
                    self._kmer_count_lens[k] += 1
                    ind += 1


if __name__ == '__main__':

    c = classifier()
    # some test sequences - first two are plasmid, third is chrom
    seqs = ['CGAGAAAGCTAAACCGTCGCCAGCCTGTTGGCCGCGATCTAAAGAACGAGGATTGCTCCGGCGGCTAGAACATTTTCGCTTGATGCCGTGCCGGCGGGCAAGCCAACATACCTGCGCGGCTTACGCCTACAGCTTCCAGCTGCTGCTCTCCTTCGCGGCGGCTCGTTTGGGGACCACTCCGTCCCGCATCGAGATCGAGCAGTTCGATGCACCGCTCCTTCTCGGCTTCCTGGAACGGCGATGTCACGTTGTCGGCAACTTAACGGGTGGCGGATAATTCGGTTGAGATGACCTCATCGACCATCAGCTACAAGAACCACCGCTTTCCACCGCAGATCATCGCCCGTGCGGTCTGGCTGTATTTCCGGTTCCCGTTGAGCCTGCGGCTGGTCGAGGAAATGTTGCTGGAGCGCGGCATCGTCGTTTCTTACGAAACGATACGGCGGTGGTGCCGCAAATTCGGGGCGGCCTACGCCAAGCAGTTGCGTAGAAAGAAACCGTCGCGGAAGGATATCTGGCATCTGGACGAGGTCGTGATTTCCATCGGCGGTCAAAAACATTGGCTTTGGCGCGCCGTTGACCAGGACGGCTACGTTCTCGATGAGATCGTCCAGACCCGCCGCGATACCAAGGCTGCCAAGCGATTGCTGGTCAGGCTGCTGAAGAAGCAAGGTCTGTCGCCGAAGCGCATCGTCACCGACAAACTGCGCTCATACGGCGCGGCAAAACGGGATGTAATGCCCGGCGTCGAACATCGATCGCACAAGGGCCTTAACAATCGCGCAGAGAATTCTCACGTGCCGCTTCGAAAACGAGAGCGGATGATGCAGGGCTTTCGATCCGTCGGAGGTTTGCAACGTTTCATCTCGGTCTTTTCAGCAGTCCGAAATCTTTTCGTCGCGCCGCACCAGAGACA']#, \
#            'GCGCGACGGGTGTCGCCGGGTGCGATGCCATAAAGAGAGCGGTACGCCCTGGAGAAATGCGCGCCGCTGGAAAAGCCGGTGGCGATGGCGATTTCGGAAATTGAAAGCGGGCTCTGCTCCAGGAGCCGCCGGGCATGTTGTAGCCGTATCAGCCGATATTGCTCCAGGAATGTCGCGCCGAGATGCGCGGCGAACAGCCGGTCGAGATGACGCGGCGTGACGCCGGCAATCCGCGCCATGGCGGCACGGTCGAGCGGCATCTCGATCGTCTCTTCCATTTTTTCGAGCACGCTCAGCAAACCCGGGTGGTGAACGCCATAACGCTCGGCGAGCGAACCGCGCTGCGGTGCCGTCGGCTCGCCGACTTCGGTGTGCAGATACCAGTCGCTGACGCGCCGGGCAAAATCGCCGCCCATGCGCTCTGATATCAAGACATGCATCATGTCGAGCGGCGCAACGCCGCCGCCGCAAGTAATTCGATTGCCGTCGATGACGAAACGTGCCTGGCGCGGCGAAAGCGCCGGAAAGGCTTCCAGCAAGGCCGGCGCATGCTCCCAATGAATGGTGAAATCGCGCTCGGAAAGAAGACCGGCGGCGGCAAGCAGATAGGGACCGCCGGAAATCCCGCCGATGCGGACGCCCTCGCGCGCCAGTTGGCGAAGGCAGGCGAGAACGGTGGGATAGTGCCAGTCCCGCGGCGATCCGCCGGCGCAGACAAAGACCGTGCCGAGGCCAGACCCGCGACCCGGAAGAGGCTCGCTCGGCACTGGCACGCCGCCCGATGATACCGCCGGAGTGCCGTCCGGCGAGAACGTCGACAATCGATAGATTTCACGCCCGGCCAGAAGATTGGCGGCGCGCAGCGGCTCGCTCACCGAAGCGTAGGACATCAGCGCAAATCCAGGGATCAATATGAAGCCGATGGTCTGGGTATTTCTTGCATCGATTGGGGCCATGTCCCCTTTCT', \
#            'CACGTCATTGATATGCCCCGGGTCCTGGTCCGACCCTAATTCGGTATTGGACAGCCATCAATGGGTGCTGACACCGTGGGGTCCGCGGATGCGGACGATGAGGAGAACCGTGGACCGGACGGCCGGCAAGGTCATGTCGATGAGGGAGGCCGTCGCCGCCTTCGTCCACGACGGGGACACCGTGAGCCTCGAAGGGTTCACCCATCTCGTCCCGACCGCCGCCGGTCACGAGATCATCCGGCAGGGCCGCCGGGATCTCACGGTCGTCCGGATGACGGCCGACATCGTGGTGGACCAGATGCTGGCGGCCGGGTGCGTGACGCGGCTGGTCTCCTCCTTCGTGGGGAACTCCTCCGCCGGTTCGCTCGGGGAGCTGCGGCGGCGGATCGAGCACGCGGACCCGGAGCCGCTGGCGTTCGAGGAGTACAGCCACTACGGGATGGTGTGCCGCTATCTCGCGGGGGCGCAGCGGCTGCCGTTCTACCCGCTGCGCTCGTACGCCGGCAGCGACCTGCCGGCCGTCAACCCCGGTATCCGCCAGGTCGTCTCCCCCTACCCGGCCGCCGACGGAGGGACCGAGCGGATCCACGTCGTGCCGCCCGTCAACCCCGACGTGACGATCGTCCACGCGCAGCGCGCCGACCGCCGGGGGAACACCCAGATCTGGGGGCTGACCGGCGTCCAGGCCGAGGCGGTGTACGCGGCGGGGAAGGCGGTCGTGGTGGTGGAGGAGCTGGTGGCGGACGAGGTGGTCCGCTCGGACCCCAATCGCACGCTGATTCCGGCGCACGCGGTCGACGCCGTGGTGGTCTGTCCGCGGGGGGCGCACCCCTCCTTCGCGCAGGGTTACTACGACCGGGACAACGCCTTCTACCGCTCCTGGTCCGCGATCAGCAAGGACCCGGCGCGGCTGCGGGCCTGGCTGGCGGAGTGGGTGTACGGGA']
    print c.classify(seqs)

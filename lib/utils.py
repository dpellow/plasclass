# Utility functions for the classifier module

complements = {'A':'T', 'C':'G', 'G':'C', 'T':'A'}

import itertools
import numpy as np
import sys


def readfq(fp): # this is a generator function
    ''' Adapted from https://github.com/lh3/readfq
    '''
    last = None # this is a buffer keeping the last unprocessed line
    while True: # mimic closure; is it a bad idea?
        if not last: # the first record or a record following a fastq
            for l in fp: # search for the start of the next record
                if l[0] in '>@': # fasta/q header line
                    last = l[:-1] # save this line
                    break
        if not last: break
        name, seqs, last = last[1:].partition(" ")[0], [], None
        for l in fp: # read the sequence
            if l[0] in '@+>':
                last = l[:-1]
                break
            seqs.append(l[:-1])
        if not last or last[0] != '+': # this is a fasta record
            yield name, ''.join(seqs), None # yield a fasta record
            if not last: break
        else: # this is a fastq record
            seq, leng, seqs = ''.join(seqs), 0, []
            for l in fp: # read the quality
                seqs.append(l[:-1])
                leng += len(l) - 1
                if leng >= len(seq): # have read enough quality
                    last = None
                    yield name, seq, ''.join(seqs); # yield a fastq record
                    break
            if last: # reach EOF before reading enough quality
                yield name, seq, None # yield a fasta record instead
                break


def get_rc(seq):
    ''' Return the reverse complement of seq
    '''
    rev = reversed(seq)
    return "".join([complements.get(i,i) for i in rev])


def count_kmers(seq, ks):
    ''' Count the k-mers in the sequence
        Return a dictionary of counts
        Assumes ks is sorted
    '''
    
    kmer_counts = {}
    for i in range(len(seq)-ks[-1] + 1):
        for k in ks:
            kmer = seq[i:i+k]
            try:
                kmer = min(kmer, get_rc(kmer))
                kmer_counts[kmer] = kmer_counts.get(kmer,0) + 1
            except:
                pass # out of alphabet

    # count the last few kmers
    end = len(ks)-1
    for i in range(len(seq)-ks[-1]+1,len(seq)-ks[0]+1):
        for k in ks[:end]:
            kmer = seq[i:i+k]
            try:
                kmer = min(kmer, get_rc(kmer))
                kmer_counts[kmer] = kmer_counts.get(kmer,0) + 1
            except:
                pass
        end -= 1

    return kmer_counts


def get_kmer_list(ks):
    ''' Return a sorted list of the lexicographically minimum of each kmer and
    its reverse complement for all the values of k in ks
    '''
    alphabet = 'ACGT'
    low_kmers = []
    for k in ks:
        all_kmers = [''.join(kmer) for kmer in itertools.product(alphabet,repeat=k)]
        for kmer in all_kmers:
            low_kmers.append(min(kmer, get_rc(kmer)))
    low_kmers.sort()
    return low_kmers


def counts2freqs(kmer_counts, ks):
    ''' Return a vector of the frequencies of k-mer occurences in a sequence
    Parameters: kmer_counts - dictionary of the kmer counts for the sequence,
                                or list of dictionaries for all the sequences
                ks - sorted list of the k-mer lengths
    Returns numpy array of frequencies
    '''
    kmers_list = get_kmer_list(ks)

    print "Converting counts to frequencies"

    kmer_freqs = []
    if isinstance(kmer_counts,dict): # single sequence
        count_list = [kmer_counts.get(kmer,0) for kmer in kmers_list]
        kmer_freqs.append(count_list)

    elif isinstance(kmer_counts,list): # set of sequences
        for k in kmer_counts:
            count_list = [k.get(kmer,0) for kmer in kmers_list]
            kmer_freqs.append(count_list)


    kmer_freqs = np.array(kmer_freqs, dtype=float)
    row_sums = np.sum(kmer_freqs,axis=1,keepdims=True)
    kmer_freqs = np.where(row_sums>0, np.divide(kmer_freqs,row_sums,dtype=float),0)
    return kmer_freqs

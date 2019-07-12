# Utility functions for the classifier module

complements = {'A':'T', 'C':'G', 'G':'C', 'T':'A'}
nt_bits = {'A':0,'C':1,'G':2,'T':3}

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


def mer2bits(kmer):
    ''' convert kmer to bit representation
    '''
    bit_mer=nt_bits[kmer[0]]
    for c in kmer[1:]:
        bit_mer = (bit_mer << 2) | nt_bits[c]
    return bit_mer


# def count_kmers(seq, ks):
#     ''' Count the k-mers in the sequence
#         Return a dictionary of counts
#         Assumes ks is sorted
#     '''
#
#     kmer_counts = {k:{} for k in ks}
#
#     k_masks = [2**(2*k)-1 for k in ks]
#     ind=0
#     bit_mers = [0 for k in ks]
#     rc_bit_mers = [0 for k in ks]
#
#     # get the first set of kmers
#     while True:
#         for i,k in enumerate(ks):
#             try:
#                 bit_mers[i] = mer2bits(seq[ind:ind+k])
#                 rc_bit_mers[i] = mer2bits(get_rc(seq[ind:ind+k]))
#                 canonical = min(bit_mers[i],rc_bit_mers[i])
#                 kmer_counts[k][canonical] = kmer_counts[k].get(canonical,0)+1
#             except:
#                 ind += 1
#                 break
#         if k == ks[-1]:
#             ind += 1
#             break
#     # count all other kmers
#     while ind<=len(seq)-ks[-1]: # iterate through sequence until last k-mer for largest k
#         for i,k in enumerate(ks):
#             try:
#                 c = nt_bits[seq[ind+k]]
#                 rc = 3^c
#                 bit_mers[i] = ((bit_mers[i]<<2)|c)&k_masks[i]
#                 rc_bit_mers[i] = ((rc_bit_mers[i]>>2)|rc<<(2*k-2))
#                 canonical = min(bit_mer,rc_bit_mer)
#                 kmer_counts[k][canonical] = kmer_counts[k].get(canonical,0)+1
#             except: # out of alphabet
#                 ind += 1
#                 # get the next set of legal kmers
#                 while ind<=len(seq)-ks[-1]:
#                     for i2,k2 in enumerate(ks):
#                         try:
#                             bit_mers[i2] = mer2bits(seq[ind:ind+k2])
#                             rc_bit_mers[i2] = mer2bits(get_rc(seq[ind:ind+k2]))
#                             canonical = min(bit_mers[i2],rc_bit_mers[i2])
#                             kmer_counts[k2][canonical] = kmer_counts[k2].get(canonical,0)+1
#                         except:
#                             ind += 1
#                             break
#                     if k2 == ks[-1]:
#                         ind += 1
#                         break
#         ind += 1 # move on to next letter in sequence
#
#     # count the last few kmers
#     end = len(ks)-1
#     for i in range(len(seq)-ks[-1]+1,len(seq)-ks[0]+1):
#         for k in ks[:end]:
#             kmer = seq[i:i+k]
#             try:
#                 canonical = min(mer2bits(kmer), mer2bits(get_rc(kmer)))
#                 kmer_counts[k][canonical] = kmer_counts[k].get(canonical,0) + 1
#             except:
#                 pass
#         end -= 1
#
#     return kmer_counts


def get_kmer_lists(ks):
    ''' Return a list of sorted lists of the lexicographically minimum of each kmer and
    its reverse complement for all the values of k in ks
    Assumes ks sorted increasing
    '''
    alphabet = 'ACGT'
    all_low_kmers = []
    for k in ks:
        low_kmers = []
        all_kmers = [''.join(kmer) for kmer in itertools.product(alphabet,repeat=k)]
        for kmer in all_kmers:
            low_kmers.append(min(mer2bits(kmer), mer2bits(get_rc(kmer))))
        low_kmers = list(set(low_kmers)) # deduplicate
        low_kmers.sort()
        all_low_kmers.append(low_kmers)
    return all_low_kmers
#
def count_kmers(seq, ks, kmer_lists, kmer_inds):
    ''' Count the k-mers in the sequence
        Return a dictionary of counts
        Assumes ks is sorted
    '''
    kmer_counts = {k:np.zeros(len(kmer_lists[i])) for i,k in enumerate(ks)}

    k_masks = [2**(2*k)-1 for k in ks]
    ind=0
    bit_mers = [0 for k in ks]
    rc_bit_mers = [0 for k in ks]

    # get the first set of kmers
    while True:
        for i,k in enumerate(ks):
            try:
                bit_mers[i] = mer2bits(seq[ind:ind+k])
                rc_bit_mers[i] = mer2bits(get_rc(seq[ind:ind+k]))
                canonical = min(bit_mers[i],rc_bit_mers[i])
                kmer_counts[k][kmer_inds[k][canonical]] += 1.
            except:
                ind += 1
                break
        if k == ks[-1]:
            break
    # count all other kmers
    while ind<len(seq)-ks[-1]: # iterate through sequence until last k-mer for largest k
        for i,k in enumerate(ks):
            try:
                c = nt_bits[seq[ind+k]]
                rc = 3^c
                bit_mers[i] = ((bit_mers[i]<<2)|c)&k_masks[i]
                rc_bit_mers[i] = ((rc_bit_mers[i]>>2)|rc<<(2*k-2))
                canonical = min(bit_mers[i],rc_bit_mers[i])
                kmer_counts[k][kmer_inds[k][canonical]] += 1.
            except: # out of alphabet
                ind += 1
                # get the next set of legal kmers
                while ind<=len(seq)-ks[-1]:
                    for i2,k2 in enumerate(ks):
                        try:
                            bit_mers[i2] = mer2bits(seq[ind:ind+k2])
                            rc_bit_mers[i2] = mer2bits(get_rc(seq[ind:ind+k2]))
                            canonical = min(bit_mers[i2],rc_bit_mers[i2])
                            kmer_counts[k2][kmer_inds[k2][canonical]] += 1.
                        except:
                            ind += 1
                            break
                    if k2 == ks[-1]:
                        ind += 1
                        break
        ind += 1 # move on to next letter in sequence

    # count the last few kmers
    end = len(ks)-1
    for i in range(len(seq)-ks[-1]+1,len(seq)-ks[0]+1):
        for k in ks[:end]:
            kmer = seq[i:i+k]
            try:
                canonical = min(mer2bits(kmer), mer2bits(get_rc(kmer)))
                kmer_counts[k][kmer_inds[k][canonical]] += 1.
            except:
                pass
        end -= 1

    return kmer_counts



def counts2freqs(kmer_counts, ks):
    ''' Return a vector of the frequencies of k-mer occurences in a sequence
    Parameters: kmer_counts - dictionary of the kmer counts for the sequence,
                              or list of dictionaries for all the sequences
                ks - sorted list of the k-mer lengths
    Returns numpy array of frequencies
    '''

    print "Converting counts to frequencies"

    if isinstance(kmer_counts,dict): # single sequence
        all_kmer_freqs = []
        for i,k in enumerate(ks):
            count_list = kmer_counts[k]
            count_sum = np.sum(count_list)
            if count_sum != 0:
                count_list = count_list/float(count_sum)
            all_kmer_freqs = np.hstack([all_kmer_freqs,count_list])

    elif isinstance(kmer_counts,list): # set of sequences
        kmer_freqs_mat = np.zeros((len(kmer_counts),sum([len(kmer_counts[0][k]) for k in ks])))
        for i,c in enumerate(kmer_counts):
            start_ind = 0
            for j,k in enumerate(ks):
                count_list = c[k]
                count_sum = np.sum(count_list)
                if count_sum != 0:
                    count_list = count_list/float(count_sum)
                kmer_freqs_mat[i,start_ind:start_ind+len(count_list)] = count_list
                start_ind += len(count_list)
    return kmer_freqs_mat


    # def counts2freqs(kmer_counts, ks):
    #     ''' Return a vector of the frequencies of k-mer occurences in a sequence
    #     Parameters: kmer_counts - dictionary of the kmer counts for the sequence,
    #                               or list of dictionaries for all the sequences
    #                 ks - sorted list of the k-mer lengths
    #     Returns numpy array of frequencies
    #     '''
    #     kmers_lists = get_kmer_lists(ks)
    #
    #     print "Converting counts to frequencies"
    #
    #     if isinstance(kmer_counts,dict): # single sequence
    #         all_kmer_freqs = []
    #         for i,k in enumerate(ks):
    #             count_list = [kmer_counts[k].get(kmer,0) for kmer in kmers_lists[i]]
    #             count_sum = sum(count_list)
    #             if count_sum != 0:
    #                 kmer_freqs = np.array(count_list,dtype=float)/float(count_sum)
    #             else:
    #                 kmer_freqs = np.zeros(len(count_list))
    #             all_kmer_freqs = np.hstack([all_kmer_freqs,kmer_freqs])
    #
    #     elif isinstance(kmer_counts,list): # set of sequences
    #         kmer_freqs_mat = np.zeros((len(kmer_counts),sum([len(klist) for klist in kmers_lists])))
    #         for i,c in enumerate(kmer_counts):
    #             all_kmer_freqs = np.array([])
    #             for j,k in enumerate(ks):
    #                 count_list = [c[k].get(kmer,0) for kmer in kmers_lists[j]]
    #                 count_sum = sum(count_list)
    #                 if count_sum != 0:
    #                     kmer_freqs = np.array(count_list,dtype=float)/float(count_sum)
    #                 else:
    #                     kmer_freqs = np.zeros(len(count_list))
    #                 all_kmer_freqs = np.hstack([all_kmer_freqs,kmer_freqs])
    #             kmer_freqs_mat[i,:] = all_kmer_freqs
    #     return kmer_freqs_mat

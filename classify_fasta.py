###
# Provide a command line script to classify sequences in a fasta file
###

import sys
sys.path.insert(0, './lib')
sys.path.insert(0, './bin')
import utils
import classifier


import argparse

def parse_user_input():

    parser = argparse.ArgumentParser(
        description=
        'classify_fasta classifies the sequences in a fasta file as plasmid origin or not'
        )
    parser.add_argument('-f','--fasta',
     help='fasta file of the sequences to be classified',
     required=True, type=str
     )
    parser.add_argument('-o','--outfile',
     help='output file prefix',
     required=False, type=str
     )

    return parser.parse_args()


def main(args):
    ''' Main function that classifies the sequences
    '''
    infile = args.fasta
    if args.outfile: outfile = args.outfile
    else: outfile = infile + '.probs.out'

    c = classifier.classifier()
    seq_names = []
    seqs = []
    print "Classifying {} in batches of 100k sequences".format(infile)
    print "Reading in the input fasta"
    i = 0
    # TODO: buffer this and run classification on batches (say ~10k sequences at a time)
    fp = open(infile)
    with open(outfile,'w') as o:
        for name, seq, _ in utils.readfq(fp):
            seq_names.append(name)
            seqs.append(seq)
            i += 1
            if i % 100000 == 0:
                print "Read {} sequences".format(i)
                probs = c.classify(seqs)
                for j,p in enumerate(probs):
                    o.write('>'+seq_names[j]+ '\n')
                    o.write(str(p) + '\n')
                seq_names = []
                seqs = []


        # last bunch of sequences:
        print "Read {} sequences".format(i)
        probs = c.classify(seqs)
        for j,p in enumerate(probs):
            o.write('>'+seq_names[j]+ '\n')
            o.write(str(p) + '\n')
    fp.close()
    print "Finished classifying"
    print "Class scores written in: {}".format(outfile)


if __name__=='__main__':
    args = parse_user_input()
    main(args)

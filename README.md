# classification
This module allows for easy classification of sequences as either plasmid or chromosomal.
For example, it can be used to classify the contigs in a (metagenomic) assembly.

## Installation

To install, download and run setup.py:

    git clone https://github.com/dpellow/classification.git
    cd classification
    python setup.py install

<!--- `classification` can also be installed using `pip`. Just do `pip install classification` --->

We recommend doing this in a virtual environment. For example, before running setup.py:
```
python -m virtualenv classification-env
source classification-env/bin/activate
```

`classification` is written in Python2.7 and requires NumPy and scikit-learn (versions compatible with Python2.7) and their dependencies. These will be installed by the setup.py script.

## Usage

The script `classify_fasta.py` can be used to classify the sequences in a fasta file:
```
python classify_fasta.py -f <fasta file> [-o <output file> default: <fasta file>.probs.out] [-p <num processes> default: 8]
```
The command line options for this script are:

`-f/--fasta`: The fasta file to be classified

`-o/--outfile`: The name of the output file. If not specified, <input filename>.probs.out

`-p/--num_processes`: The number of processes to use. Default=8

The format of the output file will be similar to the fasta format:\
Each header line starts with '>' and contains the sequence name matching that in the input file.\
Each header line is followed by a line with the plasmid score between 0 (not plasmid) and 1 (plasmid). \
The sequences will be classified in the same order as they appear in the fasta.

For example:
> \>AE015451.2\
0.13815347569215672

The classifier can also be imported and used directly in your own python code. For example, once the `classification` module has been installed you can use the following lines in your own code:
```
from classification import classifier
my_classifier = classifier()
my_classifier.classify(seqs)
```
The `classifier()` constructor takes an optional parameter of the number of processes to use for classification. The default is 1.

The sequence(s) to classify, `seqs`, can be either a single string or a list of strings. The strings must be uppercase.

The function `classifier.classify(seqs)` returns a list of plasmid scores, one per input sequence, in the same order as the input.


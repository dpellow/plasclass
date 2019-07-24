# classification
This module allows for easy classification of sequences as either plasmid or chromosomal.
For examlpe, it can be used to classify the contigs in a (metagenomic) assembly.

## Installation

To install, download and run setup.py:

    git clone https://github.com/dpellow/classification.git
    cd classification
    python setup.py install --user

We recommend doing this in a virtual environment. For example, before running setup.py:
```
    python -m virtualenv classification-env
    source classification-env/bin/activate
```

`classification` is written in Python2.7 and requires NumPy and scikit-learn (versions compatible with Python2.7) and their dependencies.

## Usage



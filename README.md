# Bioitools

A suite of python libraries for interacting with bioinformatics files.

## Installation

As a local install
```bash
$ python setup.py install --user
```
As a system install
```bash
$ python setup.py install
```

## Testing

Make sure to test the install to see if all the dependencies were correctly installed.
```bash
$ python setup.py test
```

## Usage

Bioitools installs an executable on the path consisting of two main subprograms:

### smooth

usage: `bioitools smooth [-h] [-b N] [-p F] method BEDGRAPH`

Takes in a bedgraph file, applies either hann or haar smoothing and prints the
transformed bedgraph.

#### Positional Arguments

| Argument | Description |
|----------|-------------|
| method   | smoothing method \[haar\|hann\] |
| BEDGRAPH | Bedgraph file to be smoothed |

#### Optional Arguments

| Argument | Description |
|----------|-------------|
| -b N     | Use N bins in hann smoothing \(Default 20\) |
| -p F     | Remove lower F percent of variation using after haar wavelet transform \(Default 80\) |

### fastaRemove

usage: bioitools fastaRemove [-h] REGEX FASTA

Prints fasta records from a specified fasta that don't match the given regular
expression.

positional arguments:
  REGEX       Regular expression matching fasta records to remove
  FASTA       Fasta file for input

optional arguments:
  -h, --help  show this help message and exit


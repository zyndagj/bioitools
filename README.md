# Bioitools

A suite of python libraries for interacting with bioinformatics files.

## Dependencies

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

Bioitools installs an executable on the path consisting of three main subprograms:

### normalize
```
usage: bioitools normalize [-h] [-m STR] BEDGRAPH

Normalizes a bedgraph file

positional arguments:
  BEDGRAPH              Bedgraph file to be smoothed

optional arguments:
  -h, --help            show this help message and exit
  -m STR, --method STR  Normalization method: - [RPGC] - 1x depth (reads pergenome coverage)
```

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

usage: `bioitools fastaRemove [-h] REGEX FASTA`

Prints fasta records from a specified fasta that don't match the given regular
expression.

#### Positional Arguments

| Argument | Description |
|----------|-------------|
| REGEX | Regular expression matching fasta records to remove |
| FASTA | Fasta file for input |

### repliCorr

usage: `bioitools repliCorr [-h] [-o PNG] [-s BOOL] [-r BOOL] BG [BG ...]`

repliCorr calculates the correlation of 2 re more replication bedgraph files using the [Phi Correlation coefficient](http://en.wikipedia.org/wiki/Phi_coefficient) and plots the result. 
                                                                        
#### Positional Arguments

| Argument | Description |
|----------|-------------|
| BG | Bedgraph files |
                                                                        
#### Optional Arguments

| Argument | Description |
|----------|-------------|
| -o PNG | File to write figure to (Default: figure.png) |
| -s BOOL | Save figure (Default: True) |
| -r BOOL | Render figure (Default: False) |

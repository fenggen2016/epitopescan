# epitopescan

This is a simple Python script for scanning sequences for epitopes.

# Usage

```
python epitopescan.py <fasta file> <csv file> <output file>
```

## Inputs

- `<fasta file>`: in-frame nucleotide alignment
- `<csv file>`: csv file with header:
  - the first column contains (linear) amino acid epitopes
  - the second column contains a text description of the epitope

## Outputs

A tabulated file with the name of the sequence in the first column, and one column per epitope (0 for absent, 1 for present).

# Example data

An example dataset is included of Ebola virus glycoprotein, along with a CSV file extracted from an Excel spreadsheet downloaded from the very useful [LANL HFV database](http://hfv.lanl.gov).

# To do

- Add option for regular expression searches
- Add option for position-specific (with respect to a reference sequence) epitopes

# Authors

- Selene Zarate
- Simon Frost (@sdwfrost)

This code uses a collection of utilities in `sequence.py` written by Selene Zarate; these have been modified slightly e.g. to accept different line termination characters.

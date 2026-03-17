# Extra_project

# FASTA Sequence Identity Analysis

This project provides a small Python 3 script, `analyze_fasta.py`, that reads a FASTA file and computes a pairwise k-mer identity matrix for all sequences in the file.

## Files

- `analyze_fasta.py` - main script
- `sequences.fasta` - input FASTA file for the assignment
- `environment.yml` - Conda environment file
- `LICENSE` - MIT license

## Installation

Create the Conda environment and activate it:

```bash
conda env create -f environment.yml
conda activate fasta_analysis
```

## Usage

Run with k = 1:

```bash
python analyze_fasta.py --input sequences.fasta --kmer 1
```

Run with k = 2:

```bash
python analyze_fasta.py --input sequences.fasta --kmer 2
```

Run with k = 3:

```bash
python analyze_fasta.py --input sequences.fasta --kmer 3
```

Write output to a file:

```bash
python analyze_fasta.py --input sequences.fasta --kmer 1 --output matrix.tsv
```

## Output format

The script writes a tab-separated table.

- The first row contains the column headers.
- The first column contains the sequence IDs.
- Each cell contains the pairwise identity value formatted to four decimal places.
- The diagonal is always `1.0000`.

Example layout:

```text
ID	MK278857.1	MK278840.1	MK278831.1	MK278830.1
MK278857.1	1.0000	...	...	...
MK278840.1	...	1.0000	...	...
MK278831.1	...	...	1.0000	...
MK278830.1	...	...	...	1.0000
```

## Notes

- FASTA records are parsed with `re`.
- Unexpected sequence characters trigger warnings to `stderr`.
- Pairwise identities are computed with NumPy.
- Shared k-mers are counted as a multiset intersection.

## References

- Pearson, W. R., & Lipman, D. J. (1988). Improved tools for biological sequence comparison.
- Compeau, P. E. C., Pevzner, P. A., & Tesler, G. (2011). How to apply de Bruijn graphs to genome assembly.
- Tarvin, R. D., et al. (2017). Interacting amino acid replacements allow poison frogs to evolve epibatidine resistance.

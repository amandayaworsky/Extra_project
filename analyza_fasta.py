#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import re
import sys
import argparse
import numpy as np


class FastaAnalyzer:
    """Analyze pairwise k-mer identity from a FASTA file."""

    def __init__(self, input_path: str, kmer: int, output_path: str | None = None) -> None:
        """Store arguments and validate k-mer size.

        Parameters:
            input_path: Path to the FASTA file.
            kmer: Length of k-mers.
            output_path: Optional output file path.

        Returns:
            None.
        """
        if kmer < 1:
            raise ValueError("k-mer size must be >= 1")
        self.input_path: str = input_path
        self.kmer: int = kmer
        self.output_path: str | None = output_path
        self.sequences: dict[str, str] = {}

    def parse_fasta(self) -> dict[str, str]:
        """Read and parse sequences from a FASTA file.

        Parameters:
            None.

        Returns:
            A dictionary mapping sequence IDs to sequences.
        """
        try:
            with open(self.input_path, "r", encoding="utf-8") as handle:
                text = handle.read()
        except OSError as exc:
            raise OSError(f"cannot read input file: {self.input_path}") from exc
        records = re.findall(r">([^\n]+)\n([^>]*)", text, flags=re.MULTILINE)
        for header, body in records:
            match = re.match(r"(\S+)(?:\s+(.*))?", header.strip())
            if not match:
                continue
            seq_id = match.group(1)
            sequence = re.sub(r"\s+", "", body).upper()
            bad = sorted(set(re.findall(r"[^ACGTN]", sequence, flags=re.IGNORECASE)))
            if bad:
                sys.stderr.write(f"Warning: {seq_id} contains unexpected characters: {''.join(bad)}\n")
            self.sequences[seq_id] = sequence
        if len(self.sequences) < 2:
            raise ValueError("FASTA file must contain at least two sequences")
        if any(len(seq) < self.kmer for seq in self.sequences.values()):
            raise ValueError("all sequences must be at least as long as k-mer size")
        return self.sequences

    def extract_kmers(self, sequence: str) -> np.ndarray:
        """Extract overlapping k-mers from one sequence.

        Parameters:
            sequence: Input nucleotide sequence.

        Returns:
            A NumPy array of k-mer strings.
        """
        return np.array([sequence[i:i + self.kmer] for i in range(len(sequence) - self.kmer + 1)])

    def compute_identity_matrix(self) -> tuple[list[str], np.ndarray]:
        """Compute the symmetric pairwise k-mer identity matrix.

        Parameters:
            None.

        Returns:
            A list of IDs and the identity matrix.
        """
        ids = list(self.sequences)
        n = len(ids)
        matrix = np.eye(n, dtype=float)
        counts = []
        for seq_id in ids:
            kmers = self.extract_kmers(self.sequences[seq_id])
            values, freq = np.unique(kmers, return_counts=True)
            counts.append((len(kmers), values, freq))
        for i in range(n):
            for j in range(i + 1, n):
                len_i, val_i, cnt_i = counts[i]
                len_j, val_j, cnt_j = counts[j]
                _, idx_i, idx_j = np.intersect1d(val_i, val_j, return_indices=True)
                shared = np.minimum(cnt_i[idx_i], cnt_j[idx_j]).sum()
                matrix[i, j] = matrix[j, i] = shared / min(len_i, len_j)
        return ids, matrix

    def write_output(self, ids: list[str], matrix: np.ndarray) -> None:
        """Write the identity matrix as a tab-separated table.

        Parameters:
            ids: Sequence identifiers.
            matrix: Pairwise identity matrix.

        Returns:
            None.
        """
        lines = ["ID\t" + "\t".join(ids)]
        for i, seq_id in enumerate(ids):
            row = "\t".join(f"{value:.4f}" for value in matrix[i])
            lines.append(f"{seq_id}\t{row}")
        text = "\n".join(lines) + "\n"
        if self.output_path is None:
            sys.stdout.write(text)
            return
        try:
            with open(self.output_path, "w", encoding="utf-8") as handle:
                handle.write(text)
        except OSError as exc:
            raise OSError(f"cannot write output file: {self.output_path}") from exc


def build_parser() -> argparse.ArgumentParser:
    """Create the command-line argument parser.

    Parameters:
        None.

    Returns:
        A configured ArgumentParser instance.
    """
    parser = argparse.ArgumentParser()
    parser.add_argument("-i", "--input", required=True, type=str)
    parser.add_argument("-k", "--kmer", default=1, type=int)
    parser.add_argument("-o", "--output", type=str, default=None)
    return parser


def main() -> int:
    """Run the FASTA identity analysis workflow.

    Parameters:
        None.

    Returns:
        Process exit status code.
    """
    args = build_parser().parse_args()
    try:
        analyzer = FastaAnalyzer(args.input, args.kmer, args.output)
        analyzer.parse_fasta()
        ids, matrix = analyzer.compute_identity_matrix()
        analyzer.write_output(ids, matrix)
        return 0
    except (OSError, ValueError) as exc:
        sys.stderr.write(f"Error: {exc}\n")
        return 1


if __name__ == "__main__":
    sys.exit(main())

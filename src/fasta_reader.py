from __future__ import annotations

import re
from collections import Counter
from pathlib import Path
from typing import Iterator


BV_BRC_GENOME_ID_PATTERN = re.compile(r"\|\s*(\d+\.\d+)\s*\]")


def extract_bvbrc_genome_id(header: str) -> str | None:
    """
    Extract the BV-BRC genome ID from a FASTA header.

    Example:
    >accn|JAFLRC010000001 ... [Acinetobacter bereziniae IJ5 | 106648.145]
    """
    match = BV_BRC_GENOME_ID_PATTERN.search(header)
    return match.group(1) if match else None


def iter_fasta_headers(fasta_path: str | Path) -> Iterator[str]:
    """Yield FASTA headers without loading sequence data into memory."""
    fasta_path = Path(fasta_path)

    if not fasta_path.exists():
        raise FileNotFoundError(f"FASTA file not found: {fasta_path}")

    with fasta_path.open("rt", encoding="utf-8", errors="replace") as handle:
        for line in handle:
            if line.startswith(">"):
                yield line.rstrip("\n")


def count_genomes_and_contigs(
    fasta_path: str | Path,
) -> tuple[Counter[str], int, int]:
    """
    Count contigs per BV-BRC genome ID.

    Returns:
        genome_contig_counts
        total_contigs
        headers_missing_genome_id
    """
    genome_contig_counts: Counter[str] = Counter()
    total_contigs = 0
    missing_genome_ids = 0

    for header in iter_fasta_headers(fasta_path):
        total_contigs += 1
        genome_id = extract_bvbrc_genome_id(header)

        if genome_id is None:
            missing_genome_ids += 1
        else:
            genome_contig_counts[genome_id] += 1

    return genome_contig_counts, total_contigs, missing_genome_ids

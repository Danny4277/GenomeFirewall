from __future__ import annotations

import re
from collections import Counter
from pathlib import Path
from typing import Iterator


BV_BRC_GENOME_ID_PATTERN = re.compile(r"\|\s*(\d+\.\d+)\s*\]")

SPECIES_PATTERN = re.compile(
    r"\[([^|\]]+?)\s*\|\s*\d+\.\d+\s*\]"
)

def extract_bvbrc_genome_id(header: str) -> str | None:
    """
    Extract the BV-BRC genome ID from a FASTA header.

    Example:
    >accn|JAFLRC010000001 ... [Acinetobacter bereziniae IJ5 | 106648.145]
    """
    match = BV_BRC_GENOME_ID_PATTERN.search(header)
    return match.group(1) if match else None

def extract_species_name(header: str) -> str | None:
    """
    Extract the organism name from a BV-BRC FASTA header.

    Example:
    [Acinetobacter bereziniae IJ5 | 106648.145]

    Returns:
    Acinetobacter bereziniae IJ5
    """
    match = SPECIES_PATTERN.search(header)

    if match is None:
        return None

    return match.group(1).strip()

def normalize_species_name(organism_name: str) -> str | None:
    """
    Convert a full organism/strain name into genus + species.

    Example:
    Acinetobacter bereziniae IJ5
    becomes:
    Acinetobacter bereziniae
    """
    parts = organism_name.split()

    if len(parts) < 2:
        return None

    return " ".join(parts[:2])


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

def count_genomes_by_species(
    fasta_path: str | Path,
) -> tuple[dict[str, str], Counter[str], int]:
    """
    Determine the species associated with every genome.

    Returns:
        genome_to_species:
            Maps BV-BRC genome ID to normalized species name.

        species_genome_counts:
            Counts unique genomes per species.

        conflicting_genomes:
            Number of genome IDs associated with multiple species names.
    """
    genome_species_sets: dict[str, set[str]] = {}

    for header in iter_fasta_headers(fasta_path):
        genome_id = extract_bvbrc_genome_id(header)
        organism_name = extract_species_name(header)

        if genome_id is None or organism_name is None:
            continue

        species = normalize_species_name(organism_name)

        if species is None:
            continue

        genome_species_sets.setdefault(genome_id, set()).add(species)

    genome_to_species: dict[str, str] = {}
    conflicting_genomes = 0

    for genome_id, species_names in genome_species_sets.items():
        if len(species_names) == 1:
            genome_to_species[genome_id] = next(iter(species_names))
        else:
            conflicting_genomes += 1

    species_genome_counts = Counter(genome_to_species.values())

    return (
        genome_to_species,
        species_genome_counts,
        conflicting_genomes,
    )

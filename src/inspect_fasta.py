from pathlib import Path

from src.fasta_reader import (
    count_genomes_and_contigs,
    count_genomes_by_species,
)


def main() -> None:
    fasta_path = Path(
        "data/raw/genomes/BVBRC_genome_sequence (34).fasta"
    )

    genome_counts, total_contigs, missing_ids = (
        count_genomes_and_contigs(fasta_path)
    )

    print(f"Total FASTA records/contigs: {total_contigs:,}")
    print(f"Unique BV-BRC genome IDs: {len(genome_counts):,}")
    print(f"Headers missing genome IDs: {missing_ids:,}")

    print("\nFirst 10 genome IDs:")
    for genome_id, contig_count in list(genome_counts.items())[:10]:
        print(f"  {genome_id}: {contig_count:,} contigs")

    print("\nGenomes with the most contigs:")
    for genome_id, contig_count in genome_counts.most_common(10):
        print(f"  {genome_id}: {contig_count:,} contigs")

    print("\nCounting genomes by species...")

    (
        genome_to_species,
        species_counts,
        conflicting_genomes,
    ) = count_genomes_by_species(fasta_path)

    print(f"Unique normalized species: {len(species_counts):,}")
    print(f"Genomes assigned to a species: {len(genome_to_species):,}")
    print(f"Genome IDs with conflicting species: {conflicting_genomes:,}")

    print("\nSpecies with the most genomes:")

    for species, genome_count in species_counts.most_common(20):
        print(f"  {species}: {genome_count:,} genomes")


if __name__ == "__main__":
    main()

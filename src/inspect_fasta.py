from pathlib import Path

from src.fasta_reader import count_genomes_and_contigs


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


if __name__ == "__main__":
    main()

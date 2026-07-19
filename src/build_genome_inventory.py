from pathlib import Path

import pandas as pd

from src.fasta_reader import (
    count_genomes_and_contigs,
    count_genomes_by_species,
)


def main() -> None:
    project_root = Path(__file__).resolve().parents[1]

    fasta_path = (
        project_root
        / "data"
        / "raw"
        / "genomes"
        / "BVBRC_genome_sequence (34).fasta"
    )

    output_path = (
        project_root
        / "data"
        / "processed"
        / "genome_inventory.csv"
    )

    genome_contig_counts, _, _ = count_genomes_and_contigs(
        fasta_path
    )

    genome_to_species, _, conflicting_genomes = (
        count_genomes_by_species(fasta_path)
    )

    rows = []

    for genome_id, contig_count in genome_contig_counts.items():
        rows.append(
            {
                "genome_id": genome_id,
                "species": genome_to_species.get(genome_id),
                "num_contigs": contig_count,
            }
        )

    inventory_df = pd.DataFrame(rows)
    inventory_df = inventory_df.sort_values(
        ["species", "genome_id"],
        na_position="last",
    )

    output_path.parent.mkdir(parents=True, exist_ok=True)
    inventory_df.to_csv(output_path, index=False)

    print(f"Saved inventory to: {output_path}")
    print(f"Rows: {len(inventory_df):,}")
    print(
        "Missing species:",
        f"{inventory_df['species'].isna().sum():,}",
    )
    print(f"Conflicting genome IDs: {conflicting_genomes:,}")

    print("\nTop species:")
    print(
        inventory_df["species"]
        .value_counts(dropna=False)
        .head(20)
        .to_string()
    )


if __name__ == "__main__":
    main()

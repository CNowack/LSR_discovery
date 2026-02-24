"""
group_by_species.py

Reads the NCBI RefSeq assembly summary and groups genome accessions by species
(using the species_taxid column). Splits each species into batches and writes:

  1. data/species_batches.tsv  — one row per genome with batch assignment
  2. data/batches/             — one .txt file per batch listing accessions
  3. data/species_genome_lists/— one .txt file per species listing all genome
                                 FTP paths (used by FastANI --rl)

Usage (called by Snakemake):
    python scripts/group_by_species.py \
        --summary  data/assembly_summary.txt \
        --outdir   data \
        --batch-size 20 \
        --first-batch-size 50

Or standalone:
    python scripts/group_by_species.py --help
"""

import argparse
import os
import math
import pandas as pd


# ── Constants ─────────────────────────────────────────────────────────────────

# Columns we need from the assembly summary
REQUIRED_COLS = [
    "assembly_accession",
    "species_taxid",
    "organism_name",
    "assembly_level",
    "ftp_path",
]

# Only keep assemblies with an FTP path and a usable assembly level
ASSEMBLY_LEVELS = {"Complete Genome", "Chromosome", "Scaffold", "Contig"}


# ── Helpers ───────────────────────────────────────────────────────────────────

def load_assembly_summary(path: str) -> pd.DataFrame:
    """
    Load the NCBI assembly summary TSV.
    The file has a two-line header; the second line holds column names.
    Lines starting with '#' are comments.
    """
    df = pd.read_csv(
        path,
        sep="\t",
        comment="#",
        header=0,
        low_memory=False,
    )

    # NCBI prepends '# ' to the first header line; strip if present
    df.columns = [c.lstrip("# ") for c in df.columns]

    # Validate required columns exist
    missing = [c for c in REQUIRED_COLS if c not in df.columns]
    if missing:
        raise ValueError(
            f"Assembly summary is missing expected columns: {missing}\n"
            f"Found columns: {list(df.columns[:10])} ..."
        )

    return df[REQUIRED_COLS].copy()


def filter_assemblies(df: pd.DataFrame) -> pd.DataFrame:
    """Remove rows with missing FTP paths or unusable assembly levels."""
    before = len(df)

    df = df[df["ftp_path"].notna()]
    df = df[df["ftp_path"] != "na"]
    df = df[df["assembly_level"].isin(ASSEMBLY_LEVELS)]

    after = len(df)
    print(f"  Kept {after:,} / {before:,} assemblies after filtering.")
    return df.reset_index(drop=True)


def make_fna_url(ftp_path: str) -> str:
    """
    Build the full FTP URL to the genomic FASTA file.
    NCBI FTP paths look like:
        https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/005/845/GCF_000005845.2_ASM584v2
    The genomic fna.gz file is named after the last path component.
    """
    basename = os.path.basename(ftp_path.rstrip("/"))
    return f"{ftp_path}/{basename}_genomic.fna.gz"


def assign_batches(
    accessions: list,
    first_batch_size: int,
    batch_size: int,
    species_taxid: str,
) -> list:
    """
    Split a list of accessions into batches.
    The first batch for a species is larger (50 per paper), subsequent ones
    are smaller (20 per paper).
    Returns a list of (accession, batch_id) tuples.
    """
    results = []
    remaining = list(accessions)

    # First batch
    first = remaining[:first_batch_size]
    remaining = remaining[first_batch_size:]
    if first:
        batch_id = f"{species_taxid}_batch_00"
        results.extend((acc, batch_id) for acc in first)

    # Subsequent batches
    n_subsequent = math.ceil(len(remaining) / batch_size) if remaining else 0
    for i in range(n_subsequent):
        chunk = remaining[i * batch_size : (i + 1) * batch_size]
        batch_id = f"{species_taxid}_batch_{i + 1:02d}"
        results.extend((acc, batch_id) for acc in chunk)

    return results


# ── Main ──────────────────────────────────────────────────────────────────────

def main(args):
    print(f"Loading assembly summary from: {args.summary}")
    df = load_assembly_summary(args.summary)
    df = filter_assemblies(df)

    # Add FTP URL column
    df["fna_url"] = df["ftp_path"].apply(make_fna_url)

    # ── Assign batches ────────────────────────────────────────────────────────
    print("Assigning batches by species...")
    batch_rows = []

    for species_taxid, group in df.groupby("species_taxid"):
        accessions = group["assembly_accession"].tolist()
        assigned = assign_batches(
            accessions,
            first_batch_size=args.first_batch_size,
            batch_size=args.batch_size,
            species_taxid=str(species_taxid),
        )
        for acc, batch_id in assigned:
            row = group[group["assembly_accession"] == acc].iloc[0]
            batch_rows.append({
                "accession":      acc,
                "species_taxid":  species_taxid,
                "organism_name":  row["organism_name"],
                "assembly_level": row["assembly_level"],
                "fna_url":        row["fna_url"],
                "batch_id":       batch_id,
            })

    result_df = pd.DataFrame(batch_rows)
    n_species = result_df["species_taxid"].nunique()
    n_batches = result_df["batch_id"].nunique()
    print(f"  {len(result_df):,} genomes across {n_species:,} species in {n_batches:,} batches.")

    # ── Write species_batches.tsv ─────────────────────────────────────────────
    os.makedirs(args.outdir, exist_ok=True)
    batches_tsv = os.path.join(args.outdir, "species_batches.tsv")
    result_df.to_csv(batches_tsv, sep="\t", index=False)
    print(f"Wrote: {batches_tsv}")

    # ── Write per-batch accession lists ──────────────────────────────────────
    batches_dir = os.path.join(args.outdir, "batches")
    os.makedirs(batches_dir, exist_ok=True)

    for batch_id, group in result_df.groupby("batch_id"):
        batch_file = os.path.join(batches_dir, f"{batch_id}.txt")
        group["accession"].to_csv(batch_file, index=False, header=False)

    print(f"Wrote {n_batches:,} batch files to: {batches_dir}/")

    # ── Write per-species genome lists (for FastANI --rl) ────────────────────
    # Each line is the local path where download_genomes.py will save the FASTA
    species_dir = os.path.join(args.outdir, "species_genome_lists")
    os.makedirs(species_dir, exist_ok=True)

    for species_taxid, group in result_df.groupby("species_taxid"):
        species_file = os.path.join(species_dir, f"{species_taxid}.txt")
        # Local paths follow the convention used by download_genomes.py
        local_paths = group.apply(
            lambda r: os.path.join(
                args.genomes_dir, r["batch_id"], f"{r['accession']}.fna"
            ),
            axis=1,
        )
        local_paths.to_csv(species_file, index=False, header=False)

    print(f"Wrote {n_species:,} species genome lists to: {species_dir}/")
    print("Done.")


# ── CLI ───────────────────────────────────────────────────────────────────────

if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Group RefSeq genomes by species and assign to batches."
    )
    parser.add_argument(
        "--summary",
        required=True,
        help="Path to NCBI assembly_summary.txt",
    )
    parser.add_argument(
        "--outdir",
        default="data",
        help="Output directory (default: data/)",
    )
    parser.add_argument(
        "--genomes-dir",
        default="data/genomes",
        help="Directory where download_genomes.py saves FASTAs (default: data/genomes/)",
    )
    parser.add_argument(
        "--batch-size",
        type=int,
        default=20,
        help="Genomes per batch after the first (default: 20)",
    )
    parser.add_argument(
        "--first-batch-size",
        type=int,
        default=50,
        help="Genomes in the first batch per species (default: 50)",
    )

    args = parser.parse_args()
    main(args)
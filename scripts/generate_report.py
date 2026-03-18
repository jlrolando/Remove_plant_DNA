"""
Generate summary statistics for the plant contig filtering pipeline.

Produces a TSV report with counts at each filtering stage.
"""

import sys
from Bio import SeqIO


def count_lines(filepath):
    """Count non-empty lines in a file."""
    n = 0
    with open(filepath) as fh:
        for line in fh:
            if line.strip():
                n += 1
    return n


def count_fasta(filepath):
    """Count sequences and total bp in a FASTA file."""
    n_seqs = 0
    total_bp = 0
    for record in SeqIO.parse(filepath, "fasta"):
        n_seqs += 1
        total_bp += len(record.seq)
    return n_seqs, total_bp


def load_id_set(filepath):
    """Load contig IDs from file."""
    ids = set()
    with open(filepath) as fh:
        for line in fh:
            line = line.strip()
            if line:
                ids.add(line)
    return ids


def main():
    contigs_file = snakemake.input.contigs
    comparison_file = snakemake.input.comparison
    cat_plant_all_file = snakemake.input.cat_plant_all
    cat_plant_clean_file = snakemake.input.cat_plant_clean
    cat_plant_prok_file = snakemake.input.cat_plant_prokaryotic
    kraken2_plant_file = snakemake.input.kraken2_plant
    gtdb_hits_file = snakemake.input.gtdb_hits
    confirmed_ids_file = snakemake.input.confirmed_ids
    out_report = snakemake.output.report

    # Count input contigs
    n_input, bp_input = count_fasta(contigs_file)

    # Load ID sets
    cat_all = load_id_set(cat_plant_all_file)
    cat_clean = load_id_set(cat_plant_clean_file)
    cat_prok = load_id_set(cat_plant_prok_file)
    kraken2_plant = load_id_set(kraken2_plant_file)
    gtdb_hits = load_id_set(gtdb_hits_file)
    confirmed = load_id_set(confirmed_ids_file)

    # Cross-method statistics
    cat_and_kraken = cat_clean & kraken2_plant
    cat_not_kraken = cat_clean - kraken2_plant
    kraken_not_cat = kraken2_plant - cat_all
    removed_by_gtdb = cat_clean & gtdb_hits

    with open(out_report, "w") as fh:
        fh.write("metric\tvalue\n")
        fh.write(f"input_contigs\t{n_input}\n")
        fh.write(f"input_total_bp\t{bp_input}\n")
        fh.write(f"cat_plant_contigs\t{len(cat_all)}\n")
        fh.write(f"cat_plant_removed_prokaryotic_orfs\t{len(cat_prok)}\n")
        fh.write(f"cat_plant_clean\t{len(cat_clean)}\n")
        fh.write(f"kraken2_plant_contigs\t{len(kraken2_plant)}\n")
        fh.write(f"gtdb_marker_hit_contigs\t{len(gtdb_hits)}\n")
        fh.write(f"removed_by_gtdb_from_cat_clean\t{len(removed_by_gtdb)}\n")
        fh.write(f"confirmed_plant_contigs\t{len(confirmed)}\n")
        fh.write(f"concordant_cat_kraken2\t{len(cat_and_kraken)}\n")
        fh.write(f"cat_only_not_kraken2\t{len(cat_not_kraken)}\n")
        fh.write(f"kraken2_only_not_cat\t{len(kraken_not_cat)}\n")

    # Print summary to stderr
    print("=" * 60, file=sys.stderr)
    print("PIPELINE SUMMARY", file=sys.stderr)
    print("=" * 60, file=sys.stderr)
    print(f"Input contigs:                    {n_input} ({bp_input:,} bp)", file=sys.stderr)
    print(f"CAT plant (all):                  {len(cat_all)}", file=sys.stderr)
    print(f"  Removed (prokaryotic ORFs):     {len(cat_prok)}", file=sys.stderr)
    print(f"  CAT plant (clean):              {len(cat_clean)}", file=sys.stderr)
    print(f"Kraken2 plant:                    {len(kraken2_plant)}", file=sys.stderr)
    print(f"GTDB marker hits (all contigs):   {len(gtdb_hits)}", file=sys.stderr)
    print(f"  Removed from CAT-clean:         {len(removed_by_gtdb)}", file=sys.stderr)
    print(f"CONFIRMED PLANT:                  {len(confirmed)}", file=sys.stderr)
    print(f"Concordant (CAT + Kraken2):       {len(cat_and_kraken)}", file=sys.stderr)
    print(f"CAT-only (not in Kraken2):        {len(cat_not_kraken)}", file=sys.stderr)
    print(f"Kraken2-only (not in CAT):        {len(kraken_not_cat)}", file=sys.stderr)
    print("=" * 60, file=sys.stderr)


main()

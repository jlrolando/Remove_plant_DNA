"""
Cross-reference CAT, Kraken2, and BLASTn results to produce:
  - A per-contig comparison table
  - The final list of confirmed plant contig IDs

Logic:
  Final confirmed plant = CAT-clean plant MINUS contigs with GTDB marker hits.
  Kraken2 agreement is recorded but not used as a hard filter.
"""

import sys


def load_id_set(filepath):
    """Load a file of contig IDs (one per line) into a set."""
    ids = set()
    with open(filepath) as fh:
        for line in fh:
            line = line.strip()
            if line:
                ids.add(line)
    return ids


def main():
    cat_clean = load_id_set(snakemake.input.cat_plant_clean)
    cat_prok = load_id_set(snakemake.input.cat_plant_prokaryotic)
    kraken2_plant = load_id_set(snakemake.input.kraken2_plant)
    gtdb_hits = load_id_set(snakemake.input.gtdb_hits)

    out_comparison = snakemake.output.comparison
    out_confirmed = snakemake.output.confirmed_ids

    # Union of all contigs seen across methods
    all_plant_candidates = cat_clean | cat_prok

    confirmed = []

    with open(out_comparison, "w") as fh:
        fh.write(
            "contig_id\tcat_plant\tcat_has_prokaryotic_orfs\t"
            "kraken2_plant\tgtdb_marker_hit\tfinal_verdict\n"
        )

        for contig_id in sorted(all_plant_candidates):
            is_cat_plant = "yes"
            has_prok_orfs = "yes" if contig_id in cat_prok else "no"
            is_kraken2_plant = "yes" if contig_id in kraken2_plant else "no"
            has_gtdb_hit = "yes" if contig_id in gtdb_hits else "no"

            # Determine final verdict
            if contig_id in cat_prok:
                verdict = "removed_prokaryotic_orfs"
            elif contig_id in gtdb_hits:
                verdict = "removed_gtdb_marker_hit"
            else:
                verdict = "confirmed_plant"
                confirmed.append(contig_id)

            fh.write(
                f"{contig_id}\t{is_cat_plant}\t{has_prok_orfs}\t"
                f"{is_kraken2_plant}\t{has_gtdb_hit}\t{verdict}\n"
            )

    with open(out_confirmed, "w") as fh:
        for cid in confirmed:
            fh.write(cid + "\n")

    n_concordant = len(cat_clean & kraken2_plant - gtdb_hits)
    n_discordant = len(cat_clean - kraken2_plant - gtdb_hits)

    print(f"CAT-clean plant contigs: {len(cat_clean)}", file=sys.stderr)
    print(f"Removed (prokaryotic ORFs): {len(cat_prok)}", file=sys.stderr)
    print(f"Removed (GTDB marker hits): {len(cat_clean & gtdb_hits)}", file=sys.stderr)
    print(f"Confirmed plant: {len(confirmed)}", file=sys.stderr)
    print(f"Concordant with Kraken2: {n_concordant}", file=sys.stderr)
    print(f"Discordant (CAT only, not Kraken2): {n_discordant}", file=sys.stderr)


main()

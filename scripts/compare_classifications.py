"""
Cross-reference CAT and Kraken2 results to produce:
  - A per-contig comparison table
  - The final list of confirmed plant contig IDs

Logic:
  Final confirmed plant = CAT-clean plant contigs.
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

    out_comparison = snakemake.output.comparison
    out_confirmed = snakemake.output.confirmed_ids

    # Union of all contigs seen across methods
    all_plant_candidates = cat_clean | cat_prok

    confirmed = []

    with open(out_comparison, "w") as fh:
        fh.write(
            "contig_id\tcat_plant\tcat_has_prokaryotic_orfs\t"
            "kraken2_plant\tfinal_verdict\n"
        )

        for contig_id in sorted(all_plant_candidates):
            is_cat_plant = "yes"
            has_prok_orfs = "yes" if contig_id in cat_prok else "no"
            is_kraken2_plant = "yes" if contig_id in kraken2_plant else "no"

            # Determine final verdict
            if contig_id in cat_prok:
                verdict = "removed_prokaryotic_orfs"
            else:
                verdict = "confirmed_plant"
                confirmed.append(contig_id)

            fh.write(
                f"{contig_id}\t{is_cat_plant}\t{has_prok_orfs}\t"
                f"{is_kraken2_plant}\t{verdict}\n"
            )

    with open(out_confirmed, "w") as fh:
        for cid in confirmed:
            fh.write(cid + "\n")

    n_concordant = len(cat_clean & kraken2_plant)
    n_discordant = len(cat_clean - kraken2_plant)

    print(f"CAT-clean plant contigs: {len(cat_clean)}", file=sys.stderr)
    print(f"Removed (prokaryotic ORFs): {len(cat_prok)}", file=sys.stderr)
    print(f"Confirmed plant: {len(confirmed)}", file=sys.stderr)
    print(f"Concordant with Kraken2: {n_concordant}", file=sys.stderr)
    print(f"Discordant (CAT only, not Kraken2): {n_discordant}", file=sys.stderr)


main()

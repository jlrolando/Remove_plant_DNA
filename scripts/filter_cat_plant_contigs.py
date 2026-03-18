"""
Filter CAT-classified plant contigs by inspecting ORF-level taxonomy.

Reads:
  - contig2classification.names.txt  (contig-level lineages)
  - ORF2LCA.names.txt                (ORF-level lineages)

Produces:
  - cat_plant_contigs.txt        — all contigs classified as Viridiplantae
  - cat_plant_clean.txt          — plant contigs with NO prokaryotic ORFs
  - cat_plant_prokaryotic_orfs.txt — plant contigs with ≥1 prokaryotic ORF
"""

import sys
from collections import defaultdict


def parse_contig_names(filepath):
    """Parse CAT contig2classification.names.txt.

    Returns dict: contig_id -> full_lineage_string (tab-joined taxonomy cols).
    """
    contigs = {}
    with open(filepath) as fh:
        header = fh.readline()
        for line in fh:
            line = line.rstrip("\n")
            if not line:
                continue
            cols = line.split("\t")
            contig_id = cols[0]
            # cols[1] = classification status, rest = lineage fields
            lineage = "\t".join(cols[1:])
            contigs[contig_id] = lineage
    return contigs


def parse_orf_names(filepath):
    """Parse CAT ORF2LCA.names.txt.

    Returns dict: orf_id -> full_lineage_string.
    """
    orfs = {}
    with open(filepath) as fh:
        header = fh.readline()
        for line in fh:
            line = line.rstrip("\n")
            if not line:
                continue
            cols = line.split("\t")
            orf_id = cols[0]
            lineage = "\t".join(cols[1:])
            orfs[orf_id] = lineage
    return orfs


def orf_to_contig(orf_id):
    """Derive contig ID from ORF ID (prodigal format: contigID_ORFnum)."""
    # Prodigal appends _1, _2, ... to the contig name
    parts = orf_id.rsplit("_", 1)
    return parts[0] if len(parts) == 2 and parts[1].isdigit() else orf_id


def is_viridiplantae(lineage_str):
    """Check if lineage contains Viridiplantae."""
    return "Viridiplantae" in lineage_str


def is_prokaryotic(lineage_str, keywords):
    """Check if lineage matches any prokaryotic keyword."""
    for kw in keywords:
        if kw in lineage_str:
            return True
    return False


def main():
    contig_names_file = snakemake.input.contig_names
    orf_names_file = snakemake.input.orf_names
    out_plant_all = snakemake.output.plant_all
    out_plant_clean = snakemake.output.plant_clean
    out_plant_prok = snakemake.output.plant_prokaryotic
    prok_keywords = snakemake.params.prok_keywords

    # Parse contig-level classifications
    contigs = parse_contig_names(contig_names_file)

    # Identify plant contigs
    plant_contigs = set()
    for contig_id, lineage in contigs.items():
        if is_viridiplantae(lineage):
            plant_contigs.add(contig_id)

    # Parse ORF-level classifications
    orfs = parse_orf_names(orf_names_file)

    # Group ORFs by contig
    contig_orfs = defaultdict(list)
    for orf_id, lineage in orfs.items():
        cid = orf_to_contig(orf_id)
        contig_orfs[cid].append((orf_id, lineage))

    # Inspect plant contigs for prokaryotic ORFs
    clean = []
    prokaryotic = []
    for contig_id in sorted(plant_contigs):
        has_prok = False
        for orf_id, lineage in contig_orfs.get(contig_id, []):
            if is_prokaryotic(lineage, prok_keywords):
                has_prok = True
                break
        if has_prok:
            prokaryotic.append(contig_id)
        else:
            clean.append(contig_id)

    # Write outputs
    with open(out_plant_all, "w") as fh:
        for cid in sorted(plant_contigs):
            fh.write(cid + "\n")

    with open(out_plant_clean, "w") as fh:
        for cid in clean:
            fh.write(cid + "\n")

    with open(out_plant_prok, "w") as fh:
        for cid in prokaryotic:
            fh.write(cid + "\n")

    print(f"Plant contigs (total): {len(plant_contigs)}", file=sys.stderr)
    print(f"Plant contigs (clean): {len(clean)}", file=sys.stderr)
    print(f"Plant contigs (prokaryotic ORFs): {len(prokaryotic)}", file=sys.stderr)


main()

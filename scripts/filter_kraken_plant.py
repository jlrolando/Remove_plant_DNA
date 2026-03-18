"""
Extract contig IDs classified under Viridiplantae by Kraken2.

Uses NCBI taxonomy nodes.dmp to walk up the lineage from each assigned
taxID and check whether it descends from Viridiplantae (taxid 33090).
"""

import sys


def load_taxonomy(nodes_file):
    """Load NCBI nodes.dmp into a parent dict: taxid -> parent_taxid."""
    parents = {}
    with open(nodes_file) as fh:
        for line in fh:
            parts = line.split("\t|\t")
            child = int(parts[0].strip())
            parent = int(parts[1].strip())
            parents[child] = parent
    return parents


def is_descendant_of(taxid, ancestor_taxid, parents, cache):
    """Check if taxid is a descendant of ancestor_taxid (or equal)."""
    if taxid in cache:
        return cache[taxid]

    visited = set()
    current = taxid
    while current != 1 and current not in visited:
        if current == ancestor_taxid:
            # Cache all visited nodes as True
            for node in visited:
                cache[node] = True
            cache[taxid] = True
            return True
        if current in cache:
            result = cache[current]
            for node in visited:
                cache[node] = result
            cache[taxid] = result
            return result
        visited.add(current)
        current = parents.get(current, 1)

    # Reached root without finding ancestor
    for node in visited:
        cache[node] = False
    cache[taxid] = False
    return False


def main():
    classification_file = snakemake.input.classification
    out_plant = snakemake.output.plant_contigs
    nodes_file = snakemake.params.nodes
    target_taxid = int(snakemake.params.viridiplantae_taxid)

    parents = load_taxonomy(nodes_file)
    cache = {}
    plant_contigs = []

    with open(classification_file) as fh:
        for line in fh:
            cols = line.strip().split("\t")
            if len(cols) < 3:
                continue
            classified = cols[0]  # C or U
            contig_id = cols[1]
            taxid = int(cols[2])

            if classified == "U":
                continue

            if is_descendant_of(taxid, target_taxid, parents, cache):
                plant_contigs.append(contig_id)

    with open(out_plant, "w") as fh:
        for cid in plant_contigs:
            fh.write(cid + "\n")

    print(f"Kraken2 plant contigs: {len(plant_contigs)}", file=sys.stderr)


main()

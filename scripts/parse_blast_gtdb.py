"""
Parse BLASTn output (outfmt 6) and identify contigs with significant hits
to GTDB-tk marker genes (bacterial/archaeal markers).

Any plant contig with a credible hit here is likely prokaryotic contamination.
"""

import sys


def main():
    hits_file = snakemake.input.hits
    out_flagged = snakemake.output.flagged
    min_identity = float(snakemake.params.perc_identity)
    min_length = int(snakemake.params.min_length)

    flagged_contigs = set()

    with open(hits_file) as fh:
        for line in fh:
            line = line.strip()
            if not line or line.startswith("#"):
                continue
            cols = line.split("\t")
            if len(cols) < 12:
                continue

            qseqid = cols[0]
            pident = float(cols[2])
            length = int(cols[3])

            if pident >= min_identity and length >= min_length:
                flagged_contigs.add(qseqid)

    with open(out_flagged, "w") as fh:
        for cid in sorted(flagged_contigs):
            fh.write(cid + "\n")

    print(f"Contigs with GTDB marker hits: {len(flagged_contigs)}", file=sys.stderr)


main()

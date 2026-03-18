"""
Pipeline: Filter plant contigs from seagrass root metagenomes
=============================================================
Three independent classification branches:
  1. CAT  — protein-level taxonomy (primary filter, per-ORF prokaryotic check)
  2. Kraken2 — k-mer-based validation
  3. BLASTn vs GTDB-tk marker genes — hard filter for prokaryotic markers

Final confirmed plant contigs = CAT-clean plant MINUS GTDB-marker-hit contigs.

Designed for HiPerGator (UF HPC) with SLURM resource declarations.
"""

configfile: "config.yaml"

OUTDIR = config["outdir"]


rule all:
    input:
        f"{OUTDIR}/confirmed_plant_contigs.fasta",
        f"{OUTDIR}/summary_report.tsv",
        f"{OUTDIR}/classification_comparison.tsv",


# =============================================================================
# Branch 1: CAT (Contig Annotation Tool)
# =============================================================================

rule cat_contigs:
    """Run CAT on contigs (prodigal ORF prediction + DIAMOND vs NR)."""
    input:
        contigs=config["contigs"],
    output:
        classification=f"{OUTDIR}/cat/contigs.contig2classification.txt",
        orf2lca=f"{OUTDIR}/cat/contigs.ORF2LCA.txt",
        proteins=f"{OUTDIR}/cat/contigs.predicted_proteins.faa",
        gff=f"{OUTDIR}/cat/contigs.predicted_proteins.gff",
    params:
        db=config["cat_db"],
        tax=config["cat_taxonomy"],
        prefix=f"{OUTDIR}/cat/contigs",
        r=config["cat_r"],
        f=config["cat_f"],
    threads: config["threads"]
    resources:
        mem_mb=config.get("cat_mem_mb", 64000),
        time_min=config.get("cat_time_min", 2880),
        slurm_partition=config.get("partition", "hpg-default"),
    shell:
        """
        CAT contigs \
            -c {input.contigs} \
            -d {params.db} \
            -t {params.tax} \
            -o {params.prefix} \
            --nproc {threads} \
            -r {params.r} \
            -f {params.f}
        """


rule cat_add_names:
    """Add human-readable taxonomic names to contig-level classifications."""
    input:
        classification=f"{OUTDIR}/cat/contigs.contig2classification.txt",
    output:
        f"{OUTDIR}/cat/contigs.contig2classification.names.txt",
    params:
        tax=config["cat_taxonomy"],
    resources:
        mem_mb=4000,
        time_min=30,
        slurm_partition=config.get("partition", "hpg-default"),
    shell:
        """
        CAT add_names \
            -i {input.classification} \
            -o {output} \
            -t {params.tax} \
            --only_official
        """


rule cat_add_names_orf:
    """Add human-readable taxonomic names to ORF-level classifications."""
    input:
        orf2lca=f"{OUTDIR}/cat/contigs.ORF2LCA.txt",
    output:
        f"{OUTDIR}/cat/contigs.ORF2LCA.names.txt",
    params:
        tax=config["cat_taxonomy"],
    resources:
        mem_mb=4000,
        time_min=30,
        slurm_partition=config.get("partition", "hpg-default"),
    shell:
        """
        CAT add_names \
            -i {input.orf2lca} \
            -o {output} \
            -t {params.tax} \
            --only_official
        """


rule filter_cat_plant_contigs:
    """
    Parse CAT output to identify Viridiplantae contigs, then inspect each
    ORF. Remove any plant contig where at least one ORF is classified as
    Bacteria or Archaea (chimeric or misclassified).
    """
    input:
        contig_names=f"{OUTDIR}/cat/contigs.contig2classification.names.txt",
        orf_names=f"{OUTDIR}/cat/contigs.ORF2LCA.names.txt",
    output:
        plant_all=f"{OUTDIR}/cat/cat_plant_contigs.txt",
        plant_clean=f"{OUTDIR}/cat/cat_plant_clean.txt",
        plant_prokaryotic=f"{OUTDIR}/cat/cat_plant_prokaryotic_orfs.txt",
    params:
        prok_keywords=config["prokaryotic_keywords"],
    resources:
        mem_mb=8000,
        time_min=30,
    script:
        "scripts/filter_cat_plant_contigs.py"


# =============================================================================
# Branch 2: Kraken2
# =============================================================================

rule kraken2_classify:
    """Run Kraken2 on contigs for independent k-mer-based classification."""
    input:
        contigs=config["contigs"],
    output:
        classification=f"{OUTDIR}/kraken2/kraken2_output.txt",
        report=f"{OUTDIR}/kraken2/kraken2_report.txt",
    params:
        db=config["kraken2_db"],
        confidence=config["kraken2_confidence"],
    threads: config["threads"]
    resources:
        mem_mb=config.get("kraken2_mem_mb", 64000),
        time_min=config.get("kraken2_time_min", 120),
        slurm_partition=config.get("partition", "hpg-default"),
    shell:
        """
        kraken2 \
            --db {params.db} \
            --output {output.classification} \
            --report {output.report} \
            --confidence {params.confidence} \
            --threads {threads} \
            {input.contigs}
        """


rule filter_kraken_plant:
    """Extract contig IDs classified under Viridiplantae by Kraken2."""
    input:
        classification=f"{OUTDIR}/kraken2/kraken2_output.txt",
    output:
        plant_contigs=f"{OUTDIR}/kraken2/kraken2_plant_contigs.txt",
    params:
        nodes=config["taxonomy_nodes"],
        viridiplantae_taxid=config["viridiplantae_taxid"],
    resources:
        mem_mb=4000,
        time_min=30,
    script:
        "scripts/filter_kraken_plant.py"


# =============================================================================
# Branch 3: BLASTn vs GTDB-tk marker genes
# =============================================================================

rule makeblastdb_gtdb:
    """Format GTDB-tk marker gene sequences as a BLAST nucleotide database."""
    input:
        markers=config["gtdb_markers_db"],
    output:
        nhr=f"{OUTDIR}/blast/gtdb_markers.nhr",
        nin=f"{OUTDIR}/blast/gtdb_markers.nin",
        nsq=f"{OUTDIR}/blast/gtdb_markers.nsq",
    params:
        db_prefix=f"{OUTDIR}/blast/gtdb_markers",
    resources:
        mem_mb=4000,
        time_min=30,
    shell:
        """
        makeblastdb \
            -in {input.markers} \
            -dbtype nucl \
            -out {params.db_prefix}
        """


rule blastn_gtdb:
    """BLASTn contigs against GTDB-tk marker gene database."""
    input:
        contigs=config["contigs"],
        nhr=f"{OUTDIR}/blast/gtdb_markers.nhr",
    output:
        hits=f"{OUTDIR}/blast/contigs_vs_gtdb_markers.txt",
    params:
        db_prefix=f"{OUTDIR}/blast/gtdb_markers",
        evalue=config["blast_evalue"],
    threads: config["threads"]
    resources:
        mem_mb=config.get("blast_mem_mb", 16000),
        time_min=config.get("blast_time_min", 480),
        slurm_partition=config.get("partition", "hpg-default"),
    shell:
        """
        blastn \
            -query {input.contigs} \
            -db {params.db_prefix} \
            -out {output.hits} \
            -outfmt "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore" \
            -evalue {params.evalue} \
            -num_threads {threads} \
            -max_target_seqs 5
        """


rule parse_blast_gtdb:
    """Filter BLAST hits by identity and alignment length thresholds."""
    input:
        hits=f"{OUTDIR}/blast/contigs_vs_gtdb_markers.txt",
    output:
        flagged=f"{OUTDIR}/blast/contigs_with_gtdb_hits.txt",
    params:
        perc_identity=config["blast_perc_identity"],
        min_length=config["blast_min_length"],
    resources:
        mem_mb=4000,
        time_min=10,
    script:
        "scripts/parse_blast_gtdb.py"


# =============================================================================
# Integration: cross-reference all three methods
# =============================================================================

rule compare_classifications:
    """
    Cross-reference CAT, Kraken2, and BLASTn results.
    Final confirmed plant = CAT-clean plant MINUS GTDB-marker-hit contigs.
    """
    input:
        cat_plant_clean=f"{OUTDIR}/cat/cat_plant_clean.txt",
        cat_plant_prokaryotic=f"{OUTDIR}/cat/cat_plant_prokaryotic_orfs.txt",
        kraken2_plant=f"{OUTDIR}/kraken2/kraken2_plant_contigs.txt",
        gtdb_hits=f"{OUTDIR}/blast/contigs_with_gtdb_hits.txt",
    output:
        comparison=f"{OUTDIR}/classification_comparison.tsv",
        confirmed_ids=f"{OUTDIR}/confirmed_plant_contig_ids.txt",
    resources:
        mem_mb=4000,
        time_min=10,
    script:
        "scripts/compare_classifications.py"


rule extract_confirmed_plant_contigs:
    """Extract the final set of confirmed plant contigs as FASTA."""
    input:
        contigs=config["contigs"],
        ids=f"{OUTDIR}/confirmed_plant_contig_ids.txt",
    output:
        fasta=f"{OUTDIR}/confirmed_plant_contigs.fasta",
    resources:
        mem_mb=4000,
        time_min=10,
    shell:
        """
        seqtk subseq {input.contigs} {input.ids} > {output.fasta}
        """


rule generate_report:
    """Produce summary statistics for the filtering pipeline."""
    input:
        contigs=config["contigs"],
        comparison=f"{OUTDIR}/classification_comparison.tsv",
        cat_plant_all=f"{OUTDIR}/cat/cat_plant_contigs.txt",
        cat_plant_clean=f"{OUTDIR}/cat/cat_plant_clean.txt",
        cat_plant_prokaryotic=f"{OUTDIR}/cat/cat_plant_prokaryotic_orfs.txt",
        kraken2_plant=f"{OUTDIR}/kraken2/kraken2_plant_contigs.txt",
        gtdb_hits=f"{OUTDIR}/blast/contigs_with_gtdb_hits.txt",
        confirmed_ids=f"{OUTDIR}/confirmed_plant_contig_ids.txt",
    output:
        report=f"{OUTDIR}/summary_report.tsv",
    resources:
        mem_mb=8000,
        time_min=10,
    script:
        "scripts/generate_report.py"

"""
Pipeline: Filter plant contigs from seagrass root metagenomes
=============================================================
Two independent classification branches:
  1. CAT  — protein-level taxonomy (primary filter, per-ORF prokaryotic check)
  2. Kraken2 — k-mer-based validation

Final confirmed plant contigs = CAT-clean plant contigs (Kraken2 agreement
is recorded but not used as a hard filter).

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
# Integration: cross-reference CAT and Kraken2
# =============================================================================

rule compare_classifications:
    """
    Cross-reference CAT and Kraken2 results.
    Final confirmed plant = CAT-clean plant contigs.
    """
    input:
        cat_plant_clean=f"{OUTDIR}/cat/cat_plant_clean.txt",
        cat_plant_prokaryotic=f"{OUTDIR}/cat/cat_plant_prokaryotic_orfs.txt",
        kraken2_plant=f"{OUTDIR}/kraken2/kraken2_plant_contigs.txt",
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
        confirmed_ids=f"{OUTDIR}/confirmed_plant_contig_ids.txt",
    output:
        report=f"{OUTDIR}/summary_report.tsv",
    resources:
        mem_mb=8000,
        time_min=10,
    script:
        "scripts/generate_report.py"

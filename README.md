# Remove prokaryotic DNA from seagrasas metagenome contigs

Snakemake pipeline for identifying and filtering plant (Viridiplantae) contigs from seagrass root metagenome assemblies. Uses three independent classification methods to produce a high-confidence set of plant contigs that can be used to remove plant-derived reads from future metagenomes before assembly.

## Pipeline Overview

Three classification branches run in parallel and are cross-referenced:

1. **CAT** (Contig Annotation Tool) — protein-level taxonomy via Prodigal ORF prediction + DIAMOND against NCBI NR. Per-ORF inspection removes plant contigs harboring prokaryotic ORFs (chimeric/misclassified).
2. **Kraken2** — independent k-mer-based classification. Agreement/disagreement with CAT is reported but not used as a hard filter (limited seagrass representation in Kraken2 databases).
3. **BLASTn vs GTDB-tk markers** — contigs are searched against concatenated bac120 + ar53 prokaryotic marker genes. Any plant contig with a significant hit is removed.

**Final confirmed plant contigs** = CAT-clean plant contigs MINUS contigs with GTDB marker gene hits.

```
contigs.fasta
  ├── cat_contigs ─→ cat_add_names ─→ filter_cat_plant_contigs
  │                  cat_add_names_orf ─┘
  ├── kraken2_classify ─→ filter_kraken_plant
  └── makeblastdb_gtdb ─→ blastn_gtdb ─→ parse_blast_gtdb
                                  ↓
                     compare_classifications
                                  ↓
               extract_confirmed_plant_contigs + generate_report
```

## Installation

Create the conda environment:

```bash
conda env create -f environment.yaml
conda activate remove_plant_dna
```

### Dependencies

| Tool | Version | Purpose |
|------|---------|---------|
| Snakemake | >= 7.0 | Workflow engine |
| CAT | 5.3 | Contig annotation (Prodigal + DIAMOND) |
| Kraken2 | 2.1.3 | k-mer classification |
| BLAST+ | 2.14.0 | Nucleotide alignment |
| seqtk | 1.4 | FASTA extraction |
| Python | >= 3.9 | Script runtime |
| pandas | | Data manipulation |
| BioPython | | Sequence parsing |

## Required Databases

Before running the pipeline, prepare the following databases and update paths in `config.yaml`.

### CAT database

```bash
CAT prepare --fresh -d cat_db -t cat_tax
```

### Kraken2 database

Use the Standard, PlusPF, or a custom database:

```bash
kraken2-build --standard --db kraken2_db
```

### GTDB-tk marker genes

Concatenate bac120 and ar53 marker gene nucleotide sequences from the GTDB-tk reference data package:

```bash
cat gtdbtk_data/markers/tigrfam/*.fna gtdbtk_data/markers/pfam/*.fna > gtdb_markers.fasta
```

### NCBI taxonomy

`nodes.dmp` is needed for Kraken2 lineage resolution. It is typically included in the Kraken2 database directory or can be downloaded from the NCBI taxdump:

```bash
wget https://ftp.ncbi.nlm.nih.gov/pub/taxonomy/taxdump.tar.gz
tar xzf taxdump.tar.gz nodes.dmp
```

## Configuration

Edit `config.yaml` before running. Key fields:

| Field | Description |
|-------|-------------|
| `contigs` | Path to input metagenome assembly FASTA |
| `sample_name` | Sample identifier |
| `outdir` | Output directory |
| `cat_db` / `cat_taxonomy` | CAT database and taxonomy paths |
| `kraken2_db` | Kraken2 database path |
| `gtdb_markers_db` | GTDB-tk marker genes FASTA |
| `taxonomy_nodes` | Path to NCBI `nodes.dmp` |
| `threads` | Number of threads for parallel tools (default: 16) |

### Tool parameters

| Parameter | Default | Description |
|-----------|---------|-------------|
| `cat_r` | 10 | Top DIAMOND hits per ORF |
| `cat_f` | 0.5 | Fraction of hits supporting classification |
| `kraken2_confidence` | 0.1 | Kraken2 confidence threshold |
| `blast_evalue` | 1e-10 | BLASTn e-value threshold |
| `blast_perc_identity` | 70 | Minimum percent identity for BLAST hits |
| `blast_min_length` | 100 | Minimum alignment length (bp) |

## Usage

### Local execution

```bash
snakemake --configfile config.yaml --cores 16
```

Dry run to check the DAG:

```bash
snakemake --configfile config.yaml -n
```

### HiPerGator (UF HPC)

The pipeline includes SLURM support. Each rule declares `mem_mb` and `time_min` resources and is submitted as its own SLURM job via `snakemake-executor-plugin-slurm`.

1. Edit `config.yaml` — set database paths, `slurm_account`, `qos`, and resource overrides
2. Edit `profile/config.yaml` — set `slurm_account`
3. Edit `run_hipergator.sh` — set `--account` and `--qos` in SBATCH headers

```bash
mkdir -p logs
sbatch run_hipergator.sh
```

#### Default SLURM resources

| Rule | Memory | Time |
|------|--------|------|
| `cat_contigs` | 64 GB | 48 h |
| `kraken2_classify` | 64 GB | 2 h |
| `blastn_gtdb` | 16 GB | 8 h |
| Lightweight rules | 4 GB | 10-30 min |

These can be overridden in `config.yaml` via `cat_mem_mb`, `kraken2_mem_mb`, `blast_mem_mb`, etc.

## Outputs

All outputs are written to the directory specified by `outdir` (default: `results/`).

| File | Description |
|------|-------------|
| `confirmed_plant_contigs.fasta` | Final high-confidence plant contigs (main deliverable) |
| `classification_comparison.tsv` | Per-contig table with verdicts from all three methods |
| `summary_report.tsv` | Pipeline statistics: contig counts at each filtering stage, method concordance |

### Intermediate outputs

| Directory | Contents |
|-----------|----------|
| `cat/` | CAT classification files, named lineages, plant/prokaryotic contig lists |
| `kraken2/` | Kraken2 output, report, and plant contig list |
| `blast/` | BLAST database, hit table, and flagged contig list |

## Project Structure

```
Remove_plant_DNA/
├── Snakefile                        # Pipeline rules (3 branches + integration)
├── config.yaml                      # All configurable parameters
├── environment.yaml                 # Conda dependencies
├── run_hipergator.sh                # SLURM submission script for HiPerGator
├── profile/
│   └── config.yaml                  # Snakemake SLURM executor profile
└── scripts/
    ├── filter_cat_plant_contigs.py  # Parse CAT output, flag prokaryotic ORFs
    ├── filter_kraken_plant.py       # Extract Kraken2 plant contigs via taxonomy
    ├── parse_blast_gtdb.py          # Filter BLAST hits by identity/length
    ├── compare_classifications.py   # Cross-reference methods, produce final list
    └── generate_report.py           # Summary statistics
```

#!/bin/bash
#SBATCH --job-name=filter_plant
#SBATCH --output=logs/snakemake_%j.out
#SBATCH --error=logs/snakemake_%j.err
#SBATCH --mail-type=END,FAIL
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=4gb
#SBATCH --time=96:00:00
#SBATCH --qos=your_group_account
#SBATCH --account=your_group_account
# -------------------------------------------------------
# HiPerGator submission script for plant contig filtering
#
# This is the "master" Snakemake job that sits on a single
# core and submits each rule as its own SLURM job via the
# snakemake-executor-plugin-slurm.
#
# Usage:
#   1. Edit config.yaml with your paths, account, and QOS
#   2. Edit this script: replace your_group_account, --mail-user
#   3. mkdir -p logs
#   4. sbatch run_hipergator.sh
# -------------------------------------------------------

set -euo pipefail

# --- Load modules or activate conda ---
# Option A: use conda (recommended if CAT is installed via conda)
# module load conda
# conda activate remove_plant_dna

# Option B: load HiPerGator modules individually
# module load cat/5.3
# module load kraken2/2.1.3
# module load ncbi_blast/2.14.0
# module load seqtk/1.4

# Ensure the snakemake SLURM executor plugin is installed:
#   pip install snakemake-executor-plugin-slurm

mkdir -p logs

# --- Run Snakemake ---
# The --profile flag points to profile/ which configures SLURM submission.
# Each rule is submitted as its own SLURM job with resources from the Snakefile.
snakemake \
    --configfile config.yaml \
    --profile profile/ \
    --cores 16 \
    --set-resources "cat_contigs:slurm_account=${SLURM_JOB_ACCOUNT}" \
    --set-resources "kraken2_classify:slurm_account=${SLURM_JOB_ACCOUNT}" \
    --set-resources "blastn_gtdb:slurm_account=${SLURM_JOB_ACCOUNT}"

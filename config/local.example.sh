#!/bin/bash
# Site-specific settings for the PRS pipeline.
# Copy to config/local.sh and edit. config/local.sh is git-ignored and is sourced by env/env_prs.sh.

# ---- Software environment ----
# Space-separated environment modules to load (leave empty if not using Lmod/Environment Modules)
PRS_MODULES=""
# conda.sh to source; leave empty to locate it via 'conda info --base'
CONDA_SH=""
# Name of the conda environment created from env/prs_env.lock.yml
PRS_CONDA_ENV="PRS"

# ---- Reference resources for building panels (build_panel.sh, steps 01-02) ----
# Directory holding the CanRisk PRS files (<PRS_NAME>_PRS.prs, e.g. BCAC_313_PRS.prs)
CANRISK_DIR=/path/to/canrisk
# UCSC liftover chain
CHAIN=/path/to/hg19ToHg38.over.chain.gz
# GRCh38 reference FASTA (with .fai) used to check REF alleles
REF_FASTA=/path/to/GRCh38.fa
# PGS Catalog harmonised score file used only as an independent position cross-check ('none' to skip)
PGS_CHECK=/path/to/PGS000004_hmPOS_GRCh38.txt.gz
# gnomAD v4.1 genome sites VCFs, one per chromosome
GNOMAD_DIR=/path/to/gnomad_4.1/genomes
GNOMAD_PATTERN='gnomad.genomes.v4.1.sites.{chrom}.vcf.bgz'

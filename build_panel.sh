#!/bin/bash
# Build a frozen PRS panel (steps 01-02).
# Reference resource paths are defined once in config/local.sh (see config/local.example.sh).
# Usage: [PRS_NAME=BCAC_313|BRIDGES_306] bash build_panel.sh [all|01|02]

PKG=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
PRS_NAME=${PRS_NAME:-BCAC_313}
STEP=${1:-all}

# Load environment (and config/local.sh) before strict mode (Lmod/conda are not set -eu safe)
source "$PKG/env/env_prs.sh" || exit 1
set -euo pipefail

: "${CANRISK_DIR:?set CANRISK_DIR in config/local.sh}"
: "${CHAIN:?set CHAIN in config/local.sh}"
: "${REF_FASTA:?set REF_FASTA in config/local.sh}"
: "${GNOMAD_DIR:?set GNOMAD_DIR in config/local.sh}"
PGS_CHECK=${PGS_CHECK:-none}
GNOMAD_PATTERN=${GNOMAD_PATTERN:-'gnomad.genomes.v4.1.sites.{chrom}.vcf.bgz'}

CANRISK_FILE=$CANRISK_DIR/${PRS_NAME}_PRS.prs
CORRECTIONS=$PKG/resources/manual_corrections.tsv
PANEL_DIR=$PKG/panels/$PRS_NAME
VALID_DIR=$PANEL_DIR/validation

[[ "$(command -v python)" == "$CONDA_PREFIX/bin/python" ]] || { echo "ERROR: wrong python: $(command -v python)" >&2; exit 1; }
echo "PRS      : $PRS_NAME"
echo "panel    : $PANEL_DIR"

if [[ $STEP == all || $STEP == 01 ]]; then
    bash "$PKG/bin/01_make_sites.sh" "$CANRISK_FILE" "$CHAIN" "$REF_FASTA" "$CORRECTIONS" "$PGS_CHECK" "$PANEL_DIR"
fi

if [[ $STEP == all || $STEP == 02 ]]; then
    bash "$PKG/bin/02_validate_sites.sh" "$PANEL_DIR" "$GNOMAD_DIR" "$GNOMAD_PATTERN" "$VALID_DIR"
fi

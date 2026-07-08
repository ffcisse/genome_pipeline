#!/usr/bin/env bash
#
# download_from_jgi.sh
#
# STANDALONE, OPTIONAL helper -- this is NOT wired into the Snakemake DAG and the
# pipeline never calls it. It exists purely to document how to (re-)fetch input
# FASTA from JGI Mycocosm/Phycocosm when you need to.
#
# WHY THIS IS SEPARATE FROM THE PIPELINE
# ---------------------------------------
# JGI's DownloadService needs direct access to JGI's Lucene search indexes, which
# are only reachable from NERSC Perlmutter Data Transfer Nodes (DTNs) or a JGI web
# host -- NOT from DORI, which is where the rest of this pipeline runs. So this is
# a manual step you run occasionally on Perlmutter/DTN; its output then gets
# copied/rsynced over to DORI and pointed at by config/config.yaml's
# input.protein_dir / input.cds_dir.
#
# If you already have the FASTA files -- from a previous download, a collaborator,
# wherever -- skip this script entirely and just point config.yaml at that
# directory. Nothing in the pipeline requires this script to have been run.
#
# USAGE (run on Perlmutter/DTN, NOT on DORI)
# -------------------------------------------
#   ./download_from_jgi.sh <portal-id-file> <category> <output-dir>
#
# Example:
#   ./download_from_jgi.sh portal_ids.txt proteins /pscratch/.../jgi_download/proteins
#
# <portal-id-file>
#   Plain text, one JGI portal ID per line, e.g.:
#     Cyamer1
#     CyamerSoos_1_1
#     Cyanyang1
#     Galph1_1
#     Galsul1
#     Galyel1
#     Porcrue1
#     Porpu1328_1
#     Rhomari1
#   Tip: this is exactly the species_id column of config/species.tsv, so you can
#   generate it with:  cut -f1 config/species.tsv | tail -n +2 > portal_ids.txt
#
# <category>
#   What to download. Common values for this project:
#     proteins      -- protein FASTA
#     cds           -- coding-sequence FASTA
#     sigp6_info    -- SignalP 6 predictions
#     tmhmm         -- TMHMM transmembrane predictions
#
# KEY TOOL
# --------
# gov.doe.jgi.mycocosm.download.DownloadService, invoked through the
# portal-apps-bootstrap.sh wrapper JGI provides on Perlmutter/DTN (ask JGI/Mycocosm
# support for its location if you don't already have it -- it is not part of this
# repo).
#
# IMPORTANT: --use-non-public-portals true is REQUIRED for this project, because
# these 9 genomes live on PRIVATE Mycocosm/Phycocosm portals, not the public ones.
# Omitting this flag will silently return nothing (or an auth error) for our
# species.
#
# ---------------------------------------------------------------------------

set -euo pipefail

PORTAL_ID_FILE="${1:?Usage: $0 <portal-id-file> <category> <output-dir>}"
CATEGORY="${2:?Usage: $0 <portal-id-file> <category> <output-dir>}"
OUTPUT_DIR="${3:?Usage: $0 <portal-id-file> <category> <output-dir>}"

echo "This is a documented STUB -- it does not run a real download yet."
echo "Fill in the real invocation below once portal-apps-bootstrap.sh is"
echo "available on your Perlmutter/DTN session, then remove this notice."
exit 1

# Real invocation (uncomment and adjust once portal-apps-bootstrap.sh is available):
#
# ./portal-apps-bootstrap.sh gov.doe.jgi.mycocosm.download.DownloadService \
#     --portal-id-file "${PORTAL_ID_FILE}" \
#     --category "${CATEGORY}" \
#     --output "${OUTPUT_DIR}" \
#     --info true \
#     --use-non-public-portals true

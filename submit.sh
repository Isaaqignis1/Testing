#!/bin/bash
# submit.sh - thin wrapper around `sbatch` for this pipeline.
#
# It does three things, all of which you can do by hand if you prefer
# (see README for the raw sbatch form):
#   1. sources config.sh so $LOGS_DIR is defined
#   2. makes sure the logs directory exists
#   3. calls sbatch with --export=ALL,REPO_ROOT=$PWD and the right -o / -e
#
# Usage:
#   bash submit.sh jobs/00_download.job
#   bash submit.sh jobs/01_dataset.job --dependency=afterok:1234567
#
# Any extra args after the .job path are passed straight to sbatch.

set -euo pipefail

HERE="$(cd "$(dirname "$0")" && pwd)"
source "$HERE/config.sh"
mkdir -p "$LOGS_DIR"

if [ $# -lt 1 ]; then
  echo "usage: bash submit.sh <jobs/XX_name.job> [extra sbatch args]"
  exit 2
fi

JOB="$1"; shift
if [ ! -f "$HERE/$JOB" ]; then
  echo "no such job file: $HERE/$JOB"
  exit 1
fi

sbatch \
  --export="ALL,REPO_ROOT=$HERE" \
  --output="$LOGS_DIR/%x_%j.out" \
  --error="$LOGS_DIR/%x_%j.err" \
  "$@" \
  "$HERE/$JOB"

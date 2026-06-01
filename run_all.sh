#!/usr/bin/env bash
# run_all.sh - submit every stage to SLURM in order
#
# usage:
#   bash run_all.sh              # whole pipeline
#   bash run_all.sh --from 03    # resume from stage 03 onwards
#   bash run_all.sh --only 11    # submit only stage 11
#   bash run_all.sh --dry-run    # print sbatch commands, do not submit

set -euo pipefail

REPO_ROOT="$(cd "$(dirname "$0")" && pwd)"
source "$REPO_ROOT/config.sh"
mkdir -p "$LOGS_DIR"

FROM=""
ONLY=""
DRYRUN=0
while [ $# -gt 0 ]; do
  case "$1" in
    --from)    FROM="$2"; shift 2 ;;
    --only)    ONLY="$2"; shift 2 ;;
    --dry-run) DRYRUN=1; shift ;;
    -h|--help) sed -n '2,10p' "$0"; exit 0 ;;
    *) echo "unknown arg: $1"; exit 2 ;;
  esac
done

STAGES="00_download 01_prokka 02_dataset 03_stage_phasomeit 04_phasomeit 05_parse_phasomeit 06_eggnog_db 07_eggnog_merge 08_eggnog_archaea 09_eggnog_bacteria 10_post_eggnog 11_quantify 12_tract_variation 13_figures"

submit_stage() {
  STAGE="$1"
  DEP="$2"
  JOB="$REPO_ROOT/pipeline/${STAGE}.job"
  if [ ! -f "$JOB" ]; then
    echo "missing job file: $JOB" >&2
    exit 1
  fi
  if [ -n "$DEP" ]; then
    CMD="sbatch --parsable --dependency=afterok:${DEP} ${JOB}"
  else
    CMD="sbatch --parsable ${JOB}"
  fi
  if [ "$DRYRUN" -eq 1 ]; then
    echo "  + ${CMD}" >&2
    echo "DRY_${STAGE}"
  else
    eval "$CMD"
  fi
}

PREV=""
for STAGE in $STAGES; do
  NUM="${STAGE%%_*}"
  if [ -n "$ONLY" ] && [ "$NUM" != "$ONLY" ]; then
    continue
  fi
  if [ -n "$FROM" ] && [ "$NUM" \< "$FROM" ]; then
    continue
  fi
  JID="$(submit_stage "$STAGE" "$PREV")"
  if [ -n "$PREV" ]; then
    printf "  submitted %-22s  jobid=%s  after %s\n" "$STAGE" "$JID" "$PREV"
  else
    printf "  submitted %-22s  jobid=%s\n" "$STAGE" "$JID"
  fi
  PREV="$JID"
  if [ -n "$ONLY" ]; then
    break
  fi
done

echo
echo "tail logs with:  tail -f $LOGS_DIR/*.out"

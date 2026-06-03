#!/usr/bin/env bash
# run_all.sh - pipeline orchestrator
#
# usage:
#   bash run_all.sh                        submit all pending stages in order
#   bash run_all.sh --next                 submit just the next pending stage
#   bash run_all.sh --only NN              force-submit stage NN (even if done)
#   bash run_all.sh --from NN              submit from stage NN onward
#   bash run_all.sh --dry-run              print sbatch commands, do not submit
#   bash run_all.sh --status               progress dashboard
#   bash run_all.sh --check NN             verify prereqs for one stage
#   bash run_all.sh --info [NN]            stage details (header)
#   bash run_all.sh --logs NN [--err]      tail latest log
#   bash run_all.sh --reset NN [--force]   delete marker (and outputs with --force)

set -uo pipefail

REPO_ROOT="$(cd "$(dirname "$0")" && pwd)"

STAGES="00_download 01_prokka 02_dataset 03_stage_phasomeit 04_phasomeit 05_parse_phasomeit 06_eggnog_db 07_eggnog_merge 08_eggnog_archaea 09_eggnog_bacteria 10_post_eggnog 11_quantify 12_tract_variation 13_figures"

# stage dependencies
deps_of() {
  case "$1" in
    00_download)         echo "" ;;
    01_prokka)           echo "00_download" ;;
    02_dataset)          echo "01_prokka" ;;
    03_stage_phasomeit)  echo "02_dataset" ;;
    04_phasomeit)        echo "03_stage_phasomeit" ;;
    05_parse_phasomeit)  echo "04_phasomeit" ;;
    06_eggnog_db)        echo "" ;;
    07_eggnog_merge)     echo "01_prokka" ;;
    08_eggnog_archaea)   echo "06_eggnog_db 07_eggnog_merge" ;;
    09_eggnog_bacteria)  echo "06_eggnog_db 07_eggnog_merge" ;;
    10_post_eggnog)      echo "05_parse_phasomeit 08_eggnog_archaea 09_eggnog_bacteria" ;;
    11_quantify)         echo "05_parse_phasomeit" ;;
    12_tract_variation)  echo "05_parse_phasomeit" ;;
    13_figures)          echo "10_post_eggnog 11_quantify 12_tract_variation" ;;
  esac
}

# stage-specific prereq checks (returns 0 ok, prints diagnostics on stderr)
preflight() {
  case "$1" in
    00_download)
      [ -f "$ARCHAEA_GENERA_FILE" ] || { echo "  missing: ARCHAEA_GENERA_FILE = $ARCHAEA_GENERA_FILE" >&2; return 1; }
      ;;
    01_prokka)
      [ -f "$ARCHAEA_FASTA_LIST" ] || { echo "  missing: ARCHAEA_FASTA_LIST = $ARCHAEA_FASTA_LIST  (stage 00 not run?)" >&2; return 1; }
      [ -s "$ARCHAEA_FASTA_LIST" ] || { echo "  empty:   $ARCHAEA_FASTA_LIST" >&2; return 1; }
      ;;
    02_dataset)
      [ -f "$ATB_METADATA" ] || { echo "  missing: ATB_METADATA  (stage 00)" >&2; return 1; }
      [ -d "$PROKKA_OUT" ]   || { echo "  missing: PROKKA_OUT    (stage 01)" >&2; return 1; }
      ;;
    03_stage_phasomeit)
      [ -d "$GFF_LISTS_DIR" ] || { echo "  missing: GFF_LISTS_DIR  (stage 02)" >&2; return 1; }
      ;;
    04_phasomeit)
      [ -f "$PHASOME_ARCHAEA_BASE/_genus_dirs.txt" ] || { echo "  missing: $PHASOME_ARCHAEA_BASE/_genus_dirs.txt  (stage 03)" >&2; return 1; }
      [ -d "$PHASOMEIT_REPO/.git" ] || { echo "  missing PhasomeIt clone at $PHASOMEIT_REPO  (run setup.sh)" >&2; return 1; }
      ;;
    05_parse_phasomeit)
      [ -d "$PHASOME_ARCHAEA_BASE" ] || { echo "  missing: PHASOME_ARCHAEA_BASE  (stage 04)" >&2; return 1; }
      ;;
    06_eggnog_db) ;;
    07_eggnog_merge)
      [ -d "$PROKKA_OUT" ] || { echo "  missing: PROKKA_OUT  (stage 01)" >&2; return 1; }
      ;;
    08_eggnog_archaea|09_eggnog_bacteria)
      [ -d "$EGGNOG_DB_DIR" ] && [ -n "$(ls "$EGGNOG_DB_DIR" 2>/dev/null)" ] || { echo "  missing eggnog DB at $EGGNOG_DB_DIR  (stage 06)" >&2; return 1; }
      [ -f "$EGGNOG_DIR/input/archaea_merged.faa" ] || { echo "  missing merged proteome  (stage 07)" >&2; return 1; }
      ;;
    10_post_eggnog)
      [ -f "$EGGNOG_DIR/output/archaea.emapper.annotations" ] || { echo "  missing emapper output  (stage 08)" >&2; return 1; }
      ;;
    11_quantify|12_tract_variation)
      [ -d "$EXTRACTION_DIR" ] || { echo "  missing: EXTRACTION_DIR  (stage 05)" >&2; return 1; }
      ;;
    13_figures)
      [ -d "$QUANT_DIR" ] || { echo "  missing: QUANT_DIR  (stage 11)" >&2; return 1; }
      ;;
  esac
  return 0
}

# parse args
FROM="" ; ONLY="" ; NEXT=0 ; DRYRUN=0
CMD="" ; ARG="" ; FLAG=""
while [ $# -gt 0 ]; do
  case "$1" in
    --from)    FROM="$2"; shift 2 ;;
    --only)    ONLY="$2"; shift 2 ;;
    --next)    NEXT=1; shift ;;
    --dry-run) DRYRUN=1; shift ;;
    --status)  CMD="status"; shift ;;
    --check)   CMD="check"; ARG="$2"; shift 2 ;;
    --info)
      CMD="info"
      if [ "${2:-}" != "" ] && [ "${2:0:2}" != "--" ]; then ARG="$2"; shift 2; else shift; fi
      ;;
    --logs)    CMD="logs"; ARG="$2"; shift 2 ;;
    --err)     FLAG="err"; shift ;;
    --force)   FLAG="force"; shift ;;
    --reset)   CMD="reset"; ARG="$2"; shift 2 ;;
    -h|--help) sed -n '2,14p' "$0"; exit 0 ;;
    *) echo "unknown arg: $1"; exit 2 ;;
  esac
done

need_config() {
  if [ ! -f "$REPO_ROOT/config.sh" ]; then
    echo "config.sh missing - run: cp config.sh.example config.sh" >&2
    exit 2
  fi
  # shellcheck disable=SC1091
  source "$REPO_ROOT/config.sh"
  mkdir -p "$LOGS_DIR" "$MARKERS_DIR"
}

# helpers ------------------------------------------------------------------
find_stage() {
  # accept "03" or full name; return full name on stdout, 1 if not found
  local q="$1"
  for s in $STAGES; do
    if [ "${s%%_*}" = "$q" ] || [ "$s" = "$q" ]; then echo "$s"; return 0; fi
  done
  return 1
}

marker_path() { echo "$MARKERS_DIR/$1.done"; }
marker_exists() { [ -f "$(marker_path "$1")" ]; }
has_running_job() { [ -f "$MARKERS_DIR/$1.running" ]; }

latest_log() {
  local s="$1" suffix="${2:-out}"
  ls -t "$LOGS_DIR/${s}_"*.$suffix 2>/dev/null | head -1
}

get_header_block() {
  awk '
    NR==1 && /^#!/ {next}
    /^#/ { sub(/^# ?/, ""); print; next }
    NF==0 { exit }
    { exit }
  ' "$1"
}
get_field() {
  get_header_block "$1" \
    | awk -v k="$2" 'tolower($0) ~ "^"tolower(k)":" { sub(/^[^:]*: */, ""); print; exit }'
}

green(){ printf "\033[32m%s\033[0m" "$*"; }
yellow(){ printf "\033[33m%s\033[0m" "$*"; }
red(){ printf "\033[31m%s\033[0m" "$*"; }
grey(){ printf "\033[90m%s\033[0m" "$*"; }

# verb: info ---------------------------------------------------------------
cmd_info() {
  if [ -n "$ARG" ]; then
    local s
    if ! s=$(find_stage "$ARG"); then echo "no stage '$ARG'"; exit 2; fi
    echo "=== $s ==="
    get_header_block "$REPO_ROOT/pipeline/${s}.job"
    return
  fi
  printf "%-22s  %-10s  %s\n" "STAGE" "WALLTIME" "PURPOSE"
  printf "%-22s  %-10s  %s\n" "-----" "--------" "-------"
  for s in $STAGES; do
    f="$REPO_ROOT/pipeline/${s}.job"
    printf "%-22s  %-10s  %s\n" "$s" "$(get_field "$f" Walltime)" "$(get_field "$f" Purpose)"
  done
}

# verb: status -------------------------------------------------------------
cmd_status() {
  need_config
  printf "%-22s  %-9s  %s\n" "STAGE" "STATE" "MARKER / LATEST LOG"
  printf "%-22s  %-9s  %s\n" "-----" "-----" "-------------------"
  for s in $STAGES; do
    if marker_exists "$s"; then
      m=$(marker_path "$s")
      ts=$(stat -c '%y' "$m" 2>/dev/null | cut -d. -f1)
      printf "%-22s  $(green "%-9s")  %s\n" "$s" "done" "$ts"
    else
      # check if dep markers all present
      ok=1
      for d in $(deps_of "$s"); do
        if ! marker_exists "$d"; then ok=0; break; fi
      done
      log=$(latest_log "$s" out)
      if [ "$ok" -eq 1 ]; then
        printf "%-22s  $(yellow "%-9s")  %s\n" "$s" "ready" "${log:-—}"
      else
        printf "%-22s  $(grey "%-9s")  $(grey "needs %s")\n" "$s" "blocked" "$(deps_of "$s")"
      fi
    fi
  done
}

# verb: check --------------------------------------------------------------
cmd_check() {
  need_config
  local s
  if ! s=$(find_stage "$ARG"); then echo "no stage '$ARG'"; exit 2; fi
  echo "=== checking $s ==="
  local fail=0

  # deps
  local d
  for d in $(deps_of "$s"); do
    if marker_exists "$d"; then
      printf "  $(green "✓")  dep $d (marker present)\n"
    else
      printf "  $(red "✗")  dep $d not yet done\n"
      printf "       fix: bash run_all.sh --only ${d%%_*}\n"
      fail=1
    fi
  done

  # stage-specific input/tool checks
  if preflight "$s" 2> /tmp/.pf.$$; then
    printf "  $(green "✓")  preflight checks pass\n"
  else
    while IFS= read -r ln; do printf "  $(red "✗") $ln\n"; done < /tmp/.pf.$$
    fail=1
  fi
  rm -f /tmp/.pf.$$

  if marker_exists "$s"; then
    printf "  $(yellow "!")  marker already exists for $s — stage is done. Use --reset $ARG to re-run.\n"
  fi

  if [ "$fail" -eq 0 ]; then
    printf "\n  $(green "Ready to submit:")  bash run_all.sh --only ${s%%_*}\n"
  else
    printf "\n  $(red "Fix the above first.")\n"
    exit 1
  fi
}

# verb: logs ---------------------------------------------------------------
cmd_logs() {
  need_config
  local s
  if ! s=$(find_stage "$ARG"); then echo "no stage '$ARG'"; exit 2; fi
  local sfx=out; [ "$FLAG" = "err" ] && sfx=err
  local f; f=$(latest_log "$s" $sfx)
  if [ -z "$f" ]; then echo "no $sfx log found for $s in $LOGS_DIR"; exit 2; fi
  echo "=== tail $f ==="
  tail -100 "$f"
}

# verb: reset --------------------------------------------------------------
cmd_reset() {
  need_config
  local s
  if ! s=$(find_stage "$ARG"); then echo "no stage '$ARG'"; exit 2; fi
  local m; m=$(marker_path "$s")
  if [ ! -f "$m" ]; then echo "no marker to reset for $s"; exit 0; fi
  rm -f "$m"
  echo "removed marker: $m"
  if [ "$FLAG" = "force" ]; then
    echo "(--force: outputs are NOT auto-deleted - delete by hand if needed)"
  fi
}

# verb: submit -------------------------------------------------------------
sbatch_run() {
  if [ "$DRYRUN" -eq 1 ]; then echo "DRY  $*" >&2; echo "DRY_$1"
  else sbatch --parsable "$@"; fi
}

submit_one() {
  local stage="$1" dep_jid="$2"
  local job="$REPO_ROOT/pipeline/${stage}.job"
  local marker="$(marker_path "$stage")"
  local lout="$LOGS_DIR/${stage}_%A_%a.out"
  local lerr="$LOGS_DIR/${stage}_%A_%a.err"

  if ! preflight "$stage" 2>/dev/null; then
    echo "  $stage preflight failed - run 'bash run_all.sh --check ${stage%%_*}' for details" >&2
    return 1
  fi

  local args=(--output="$lout" --error="$lerr" "$job")
  [ -n "$dep_jid" ] && args=(--dependency=afterok:"$dep_jid" "${args[@]}")

  local jid; jid=$(sbatch_run "${args[@]}")

  # marker job: depend on jid, just touch the marker
  local mark_args=(--dependency=afterok:"$jid" --output="$LOGS_DIR/marker_${stage}_%j.out" --error="$LOGS_DIR/marker_${stage}_%j.err" --wrap="mkdir -p '$MARKERS_DIR' && touch '$marker' && echo done")
  local mjid; mjid=$(sbatch_run "${mark_args[@]}")

  printf "  submitted %-22s  jobid=%-12s  marker=%s\n" "$stage" "$jid" "$mjid"
  echo "$mjid"   # downstream stages depend on the MARKER job, not the work job
}

cmd_submit() {
  need_config
  local prev_jid=""
  for stage in $STAGES; do
    local num="${stage%%_*}"

    # filter by --only / --from / --next
    if [ -n "$ONLY" ] && [ "$num" != "$ONLY" ]; then continue; fi
    if [ -n "$FROM" ] && [ "$num" \< "$FROM" ]; then continue; fi

    # skip already-done unless --only forces it
    if marker_exists "$stage" && [ -z "$ONLY" ]; then
      echo "  $stage already done (marker present); skipping"
      prev_jid=""
      continue
    fi

    if [ -n "$ONLY" ] && marker_exists "$stage"; then
      echo "  forcing $stage despite existing marker (use --reset to clear first)"
    fi

    if ! mjid=$(submit_one "$stage" "$prev_jid"); then
      echo "  ABORT at $stage"; exit 1
    fi
    prev_jid="$mjid"

    if [ -n "$ONLY" ] || [ "$NEXT" -eq 1 ]; then break; fi
  done

  echo
  echo "logs   : $LOGS_DIR"
  echo "markers: $MARKERS_DIR"
  echo "watch  : bash run_all.sh --status"
}

# dispatch -----------------------------------------------------------------
case "$CMD" in
  info)   cmd_info ;;
  status) cmd_status ;;
  check)  cmd_check ;;
  logs)   cmd_logs ;;
  reset)  cmd_reset ;;
  *)      cmd_submit ;;
esac

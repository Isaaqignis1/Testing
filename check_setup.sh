#!/usr/bin/env bash
# check_setup.sh — single setup verifier.
#
# Confirms:
#   (1) config.sh parses without error
#   (2) standards files exist
#   (3) scratch paths are writable
#   (4) conda is on PATH
#   (5) each named conda env exists AND its key tool runs (--version)
#   (6) PhasomeIt is cloned AND its top-level script runs (--help)
#   (7) the ATB index URL is reachable (HEAD request)
#
# usage:
#   bash check_setup.sh             # full run
#   bash check_setup.sh --quick     # skip network + tool invocations

set -u  # not -e: we want every check to run even after one fails

REPO_ROOT="$(cd "$(dirname "$0")" && pwd)"

# ---- argument parsing -------------------------------------------------------
QUICK=0
for a in "$@"; do
  case "$a" in
    --quick) QUICK=1 ;;
    -h|--help) sed -n '2,12p' "$0"; exit 0 ;;
    *) echo "unknown arg: $a"; exit 2 ;;
  esac
done

# ---- output helpers ---------------------------------------------------------
N_PASS=0; N_WARN=0; N_FAIL=0
green(){  printf "  \033[32mPASS\033[0m  %s\n" "$*"; N_PASS=$((N_PASS+1)); }
yellow(){ printf "  \033[33mWARN\033[0m  %s\n" "$*"; N_WARN=$((N_WARN+1)); }
red(){    printf "  \033[31mFAIL\033[0m  %s\n" "$*"; N_FAIL=$((N_FAIL+1)); }
section(){ printf "\n== %s ==\n" "$*"; }

# ---- 1. config.sh -----------------------------------------------------------
section "1. config.sh"
if [ ! -f "$REPO_ROOT/config.sh" ]; then
  red "config.sh missing — run: cp config.sh.example config.sh"
  echo
  echo "summary: $N_PASS pass, $N_WARN warn, $N_FAIL fail"
  exit 1
fi
# shellcheck disable=SC1091
if source "$REPO_ROOT/config.sh" 2>/dev/null; then
  green "config.sh sources cleanly"
else
  red "config.sh has a shell syntax error"
  exit 1
fi

# ---- 2. standards files -----------------------------------------------------
section "2. inputs the user can change (standards)"
for f in "$ARCHAEA_GENERA_FILE" "$CONTROLS_SPECIES_TSV" "$REPO_ROOT/standards/contaminants.txt" "$REPO_ROOT/standards/colors.py"; do
  if [ -f "$f" ]; then
    green "exists: $f"
  else
    red "missing: $f"
  fi
done

# ---- 3. scratch paths -------------------------------------------------------
section "3. scratch paths (created if missing, must be writable)"
for V in SCRATCH_ROOT INPUTS_DIR OUTPUTS_DIR LOGS_DIR ATB_TARBALL_CACHE; do
  P="${!V}"
  if mkdir -p "$P" 2>/dev/null && [ -w "$P" ]; then
    green "writable: $V = $P"
  else
    red "NOT writable: $V = $P"
  fi
done

# ---- 4. conda ---------------------------------------------------------------
section "4. conda"
if [ -n "$CONDA_MODULE" ]; then
  if command -v module >/dev/null 2>&1; then
    if module load "$CONDA_MODULE" 2>/dev/null; then
      green "module load $CONDA_MODULE"
    else
      yellow "module load $CONDA_MODULE failed (ok if conda is already on PATH)"
    fi
  else
    yellow "'module' command not available — relying on conda being on PATH"
  fi
fi

if command -v conda >/dev/null 2>&1; then
  CB="$(conda info --base 2>/dev/null)"
  green "conda found at $(command -v conda)"
  green "conda base: $CB"
  # shellcheck disable=SC1091
  source "$CB/etc/profile.d/conda.sh" 2>/dev/null || true
else
  red "conda not on PATH — fix CONDA_MODULE in config.sh, or load conda manually before running this check"
  echo; echo "summary: $N_PASS pass, $N_WARN warn, $N_FAIL fail"
  exit 1
fi

# ---- 5. conda envs + key tools ---------------------------------------------
section "5. conda envs + tool versions"

check_env() {
  ENV_NAME="$1"
  TOOL_CMD="$2"          # command to run inside the env
  EXPECT_RE="$3"         # extended regex that should match output
  LABEL="$4"

  if ! conda activate "$ENV_NAME" 2>/dev/null; then
    red "env '$ENV_NAME' does not exist (conda activate failed)"
    return
  fi
  green "env '$ENV_NAME' activates"
  if [ "$QUICK" -eq 1 ]; then
    conda deactivate
    return
  fi
  OUT="$(eval "$TOOL_CMD" 2>&1 | head -3 || true)"
  if echo "$OUT" | grep -Eq "$EXPECT_RE"; then
    green "$LABEL works: $(echo "$OUT" | head -1)"
  else
    red "$LABEL did not produce expected output. saw:"
    echo "$OUT" | sed 's/^/         /'
  fi
  conda deactivate
}

check_env "$ENV_PROKKA"    "prokka --version"               "[Pp]rokka"        "prokka"
check_env "$ENV_PHASOMEIT" "python -c 'import Bio,pandas,numpy,matplotlib; print(Bio.__version__)'"  "[0-9]+\." "biopython+pandas"
check_env "$ENV_EGGNOG"    "emapper.py --version"           "emapper-[0-9]"    "eggnog-mapper"

# ---- 6. PhasomeIt ----------------------------------------------------------
section "6. PhasomeIt"
if [ ! -d "$PHASOMEIT_REPO/.git" ]; then
  red "PhasomeIt not cloned at $PHASOMEIT_REPO (run setup.sh)"
elif [ ! -f "$PHASOMEIT_REPO/phasomeit.py" ]; then
  red "PhasomeIt repo present but phasomeit.py missing at $PHASOMEIT_REPO"
else
  green "PhasomeIt repo at $PHASOMEIT_REPO"
  if [ "$QUICK" -eq 0 ]; then
    conda activate "$ENV_PHASOMEIT" 2>/dev/null || true
    if python "$PHASOMEIT_REPO/phasomeit.py" -h 2>&1 | head -1 | grep -q -i "phasome\|usage"; then
      green "phasomeit.py --help runs"
    else
      yellow "phasomeit.py exists but --help did not produce expected output"
    fi
    conda deactivate 2>/dev/null || true
  fi
fi

# ---- 7. network -------------------------------------------------------------
section "7. network (ATB index reachable)"
if [ "$QUICK" -eq 1 ]; then
  yellow "skipped (--quick)"
else
  URL="https://osf.io/download/4yv85/"
  if command -v curl >/dev/null 2>&1; then
    if curl -sI -L --max-time 30 "$URL" | head -1 | grep -q "200\|302"; then
      green "ATB index reachable: $URL"
    else
      red "ATB index not reachable: $URL"
    fi
  else
    yellow "curl not installed, skipping network check"
  fi
fi

# ---- summary ----------------------------------------------------------------
echo
echo "==========================================================="
printf "  %d pass, %d warn, %d fail\n" "$N_PASS" "$N_WARN" "$N_FAIL"
echo "==========================================================="
if [ "$N_FAIL" -eq 0 ]; then
  echo "  Ready to submit. Next:"
  echo "    bash run_all.sh --dry-run     # preview"
  echo "    bash run_all.sh               # for real"
  exit 0
else
  echo "  Fix the FAILs above, then re-run this script."
  exit 1
fi

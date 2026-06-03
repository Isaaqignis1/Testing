#!/usr/bin/env bash
# check_setup.sh - one-shot setup verifier.
# Confirms config parses, scratch is writable, conda envs exist and their
# key tools actually run, PhasomeIt is cloned and runnable, ATB is reachable.
#
# usage:
#   bash check_setup.sh           # full
#   bash check_setup.sh --quick   # skip tool/network invocations

set -u

REPO_ROOT="$(cd "$(dirname "$0")" && pwd)"
QUICK=0
for a in "$@"; do
  case "$a" in
    --quick) QUICK=1 ;;
    -h|--help) sed -n '2,9p' "$0"; exit 0 ;;
    *) echo "unknown arg: $a"; exit 2 ;;
  esac
done

N_PASS=0; N_WARN=0; N_FAIL=0
green(){  printf "  \033[32mPASS\033[0m  %s\n" "$*"; N_PASS=$((N_PASS+1)); }
yellow(){ printf "  \033[33mWARN\033[0m  %s\n" "$*"; N_WARN=$((N_WARN+1)); }
red(){    printf "  \033[31mFAIL\033[0m  %s\n" "$*"; N_FAIL=$((N_FAIL+1)); }
section(){ printf "\n== %s ==\n" "$*"; }

# 1. config.sh
section "1. config.sh"
if [ ! -f "$REPO_ROOT/config.sh" ]; then
  red "config.sh missing - run: cp config.sh.example config.sh"
  exit 1
fi
if source "$REPO_ROOT/config.sh" 2>/dev/null; then
  green "config.sh sources cleanly"
else
  red "config.sh has a shell syntax error"; exit 1
fi

# 2. standards
section "2. inputs the user can change (standards)"
for f in "$ARCHAEA_GENERA_FILE" "$CONTROLS_SPECIES_TSV" \
         "$REPO_ROOT/standards/contaminants.txt" "$REPO_ROOT/standards/colors.py"; do
  if [ -f "$f" ]; then green "exists: $f"
  else                  red "missing: $f"; fi
done

# 3. scratch
section "3. scratch paths (created if missing, must be writable)"
for V in SCRATCH_ROOT INPUTS_DIR OUTPUTS_DIR LOGS_DIR ATB_TARBALL_CACHE; do
  P="${!V}"
  if mkdir -p "$P" 2>/dev/null && [ -w "$P" ]; then
    green "writable: $V = $P"
  else
    red "NOT writable: $V = $P"
  fi
done

# 4. conda
section "4. conda"
if [ -n "$CONDA_MODULE" ]; then
  if command -v module >/dev/null 2>&1; then
    if module load "$CONDA_MODULE" 2>/dev/null; then
      green "module load $CONDA_MODULE"
    else
      yellow "module load $CONDA_MODULE failed (ok if conda already on PATH)"
    fi
  fi
fi
if command -v conda >/dev/null 2>&1; then
  CB="$(conda info --base 2>/dev/null)"
  green "conda found at $(command -v conda)"
  green "conda base: $CB"
  source "$CB/etc/profile.d/conda.sh" 2>/dev/null || true
else
  red "conda not on PATH"; exit 1
fi

# 5. envs + tools
section "5. conda envs + tool versions"

check_env_activate() {
  local NAME="$1"
  if ! conda activate "$NAME" 2>/dev/null; then
    red "env '$NAME' does not exist or cannot be activated"
    return 1
  fi
  green "env '$NAME' activates"
  return 0
}

# Prokka
if check_env_activate "$ENV_PROKKA"; then
  if [ "$QUICK" -eq 0 ]; then
    OUT="$(prokka --version 2>&1 | head -1 || true)"
    if echo "$OUT" | grep -Eq "prokka [0-9]"; then
      green "prokka --version: $OUT"
    else
      red "prokka --version did not run cleanly. saw: $OUT"
    fi
  fi
  conda deactivate 2>/dev/null || true
fi

# Phasome / analysis env — must have biopython + pandas + numpy + matplotlib
if check_env_activate "$ENV_PHASOMEIT"; then
  if [ "$QUICK" -eq 0 ]; then
    OUT="$(python -c 'import Bio, pandas, numpy, matplotlib; print("Bio", Bio.__version__, "pd", pandas.__version__)' 2>&1)"
    if echo "$OUT" | grep -q "^Bio "; then
      green "biopython+pandas+numpy+matplotlib: $OUT"
    else
      red "env '$ENV_PHASOMEIT' missing python deps."
      echo "         saw: $OUT" | head -3
      echo "         FIX:  conda env update -n $ENV_PHASOMEIT -f envs/phasome.yml"
      echo "         or:   edit config.sh and point ENV_PHASOMEIT at a working env"
    fi
  fi
  conda deactivate 2>/dev/null || true
fi

# Eggnog - emapper.py without a db prints an error, but if the binary runs
# at all the env is fine (the db is downloaded by stage 06).
if check_env_activate "$ENV_EGGNOG"; then
  if [ "$QUICK" -eq 0 ]; then
    if command -v emapper.py >/dev/null 2>&1; then
      green "emapper.py on PATH at $(command -v emapper.py)"
      # version check only — db is fine to be missing pre-stage-06
      VER="$(emapper.py --version 2>&1 | grep -Eo 'emapper-[0-9.]+' | head -1)"
      if [ -n "$VER" ]; then
        green "emapper version: $VER  (DB download is stage 06, expected missing)"
      else
        yellow "emapper.py runs but no version string seen — proceed anyway"
      fi
    else
      red "emapper.py not on PATH inside env '$ENV_EGGNOG'"
    fi
  fi
  conda deactivate 2>/dev/null || true
fi

# 6. PhasomeIt
section "6. PhasomeIt"
if [ ! -d "$PHASOMEIT_REPO/.git" ]; then
  red "PhasomeIt not cloned at $PHASOMEIT_REPO (run setup.sh)"
elif [ ! -f "$PHASOMEIT_REPO/phasomeit.py" ]; then
  red "PhasomeIt repo present but phasomeit.py missing at $PHASOMEIT_REPO"
else
  green "PhasomeIt repo at $PHASOMEIT_REPO"
  if [ "$QUICK" -eq 0 ]; then
    conda activate "$ENV_PHASOMEIT" 2>/dev/null || true
    if python "$PHASOMEIT_REPO/phasomeit.py" -h 2>&1 | head -3 | grep -qi "phasome\|usage"; then
      green "phasomeit.py -h runs"
    else
      yellow "phasomeit.py exists but -h gave unexpected output"
    fi
    conda deactivate 2>/dev/null || true
  fi
fi

# 7. network
section "7. network (ATB index reachable)"
if [ "$QUICK" -eq 1 ]; then
  yellow "skipped (--quick)"
else
  URL="https://osf.io/download/argkq/"
  # OSF rejects HEAD on download endpoints; use a 1-byte range GET instead.
  if command -v curl >/dev/null 2>&1; then
    HTTP_CODE="$(curl -sSL -o /dev/null -w '%{http_code}' \
                       --max-time 60 -r 0-0 "$URL" 2>/dev/null || echo 000)"
    if [ "$HTTP_CODE" = "200" ] || [ "$HTTP_CODE" = "206" ]; then
      green "ATB index reachable (HTTP $HTTP_CODE)"
    else
      red "ATB index not reachable (HTTP $HTTP_CODE). Try: curl -L --max-time 60 -o /tmp/atb_test '$URL'"
    fi
  else
    yellow "curl not installed, skipping"
  fi
fi

echo
echo "==========================================================="
printf "  %d pass, %d warn, %d fail\n" "$N_PASS" "$N_WARN" "$N_FAIL"
echo "==========================================================="
if [ "$N_FAIL" -eq 0 ]; then
  echo "  Ready to submit. Next:"
  echo "    bash run_all.sh --dry-run"
  echo "    bash run_all.sh"
  exit 0
else
  echo "  Fix the FAILs above, then re-run."
  exit 1
fi

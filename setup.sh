#!/usr/bin/env bash
# setup.sh - one-off setup. Safe to re-run.
#   1. write config.sh from the example (only if missing)
#   2. clone PhasomeIt into PHASOMEIT_REPO (only if missing)
#   3. create the three conda envs from envs/*.yml, OR update them to match
#      the yaml if they already exist (yaml is the source of truth).
#
# Conda envs live wherever conda already keeps them, NOT inside the repo.
# Data dirs live on scratch (see config.sh), NOT inside the repo.
#
# usage:
#   bash setup.sh
#   bash setup.sh --skip-envs        # only clone PhasomeIt + write config
#   bash setup.sh --skip-clone       # only handle envs
#   bash setup.sh --no-update-envs   # create missing envs but DON'T touch existing ones

set -euo pipefail

REPO_ROOT="$(cd "$(dirname "$0")" && pwd)"
cd "$REPO_ROOT"

SKIP_ENVS=0
SKIP_CLONE=0
NO_UPDATE_ENVS=0
for a in "$@"; do
  case "$a" in
    --skip-envs)      SKIP_ENVS=1 ;;
    --skip-clone)     SKIP_CLONE=1 ;;
    --no-update-envs) NO_UPDATE_ENVS=1 ;;
    -h|--help)        sed -n '2,14p' "$0"; exit 0 ;;
    *) echo "unknown arg: $a"; exit 2 ;;
  esac
done

# 1. config.sh
if [ ! -f config.sh ]; then
  cat > config.sh <<'EOF'
#!/bin/bash
# config.sh - your local pipeline configuration.
# defaults come from config.sh.example (sourced first; never edit it).
# put YOUR overrides below the line. anything not overridden uses the default.
#
# examples:
#   SCRATCH_ROOT="${HOME}/somewhere_else"
#   MIN_PER_SPECIES=3
#   ENV_PHASOMEIT="ar"

source "$(dirname "${BASH_SOURCE[0]}")/config.sh.example"

# ===== YOUR OVERRIDES BELOW =====

EOF
  echo "[setup] created minimal config.sh that sources the example."
  echo "        add overrides below the marker line if you need any."
else
  echo "[setup] config.sh already exists"
fi
source config.sh

# 2. PhasomeIt clone
if [ "$SKIP_CLONE" -eq 0 ]; then
  if [ -d "$PHASOMEIT_REPO/.git" ]; then
    echo "[setup] PhasomeIt already cloned at $PHASOMEIT_REPO"
  else
    echo "[setup] cloning PhasomeIt -> $PHASOMEIT_REPO"
    mkdir -p "$(dirname "$PHASOMEIT_REPO")"
    git clone "$PHASOMEIT_GIT_URL" "$PHASOMEIT_REPO"
    chmod +x "$PHASOMEIT_REPO/bossref" 2>/dev/null || true
  fi
fi

# 3. conda envs
if [ "$SKIP_ENVS" -eq 0 ]; then
  if [ -n "$CONDA_MODULE" ]; then
    module purge 2>/dev/null || true
    module load "$CONDA_MODULE" 2>/dev/null || {
      echo "[setup] WARN: 'module load $CONDA_MODULE' failed."
      echo "             If conda is already on PATH this is fine."
    }
  fi
  if ! command -v conda >/dev/null 2>&1; then
    echo "[setup] ERROR: conda not on PATH."
    echo "        Fix CONDA_MODULE in config.sh, or load conda manually first."
    exit 1
  fi
  source "$(conda info --base)/etc/profile.d/conda.sh"

  ensure_env() {
    local name="$1"
    local yml="$2"
    if conda env list | awk '{print $1}' | grep -qx "$name"; then
      if [ "$NO_UPDATE_ENVS" -eq 1 ]; then
        echo "[setup] env '$name' exists - leaving alone (--no-update-envs)"
      else
        echo "[setup] env '$name' exists - syncing to $yml (conda env update --prune)"
        conda env update -n "$name" -f "$REPO_ROOT/$yml" --prune
      fi
    else
      echo "[setup] creating env '$name' from $yml"
      conda env create -n "$name" -f "$REPO_ROOT/$yml"
    fi
  }
  ensure_env "$ENV_PROKKA"    "envs/prokka.yml"
  ensure_env "$ENV_PHASOMEIT" "envs/phasome.yml"
  ensure_env "$ENV_EGGNOG"    "envs/eggnog.yml"
fi

# 4. ensure scratch dirs
mkdir -p "$SCRATCH_ROOT" "$INPUTS_DIR" "$OUTPUTS_DIR" "$LOGS_DIR" "$ATB_TARBALL_CACHE"

echo
echo "[setup] done."
echo "        next:  bash check_setup.sh"
echo "        then:  bash run_all.sh --dry-run"
echo "        then:  bash run_all.sh"

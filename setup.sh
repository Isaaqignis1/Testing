#!/usr/bin/env bash
# one-off setup: create config.sh, clone PhasomeIt, create conda envs
# safe to re-run — every step is idempotent.
#
# usage: bash setup.sh
#        bash setup.sh --skip-envs    # only clone PhasomeIt + write config
#        bash setup.sh --skip-clone   # only create envs

set -euo pipefail

REPO_ROOT="$(cd "$(dirname "$0")" && pwd)"
cd "${REPO_ROOT}"

SKIP_ENVS=0; SKIP_CLONE=0
for a in "$@"; do
  case "$a" in
    --skip-envs)  SKIP_ENVS=1 ;;
    --skip-clone) SKIP_CLONE=1 ;;
    -h|--help)    sed -n '2,8p' "$0"; exit 0 ;;
    *) echo "unknown arg: $a"; exit 2 ;;
  esac
done

# ---------- 1. config.sh -----------------------------------------------------
if [[ ! -f config.sh ]]; then
  cp config.sh.example config.sh
  echo "[setup] created config.sh from config.sh.example"
else
  echo "[setup] config.sh already exists, leaving as-is"
fi
# shellcheck disable=SC1091
source config.sh

# ---------- 2. PhasomeIt clone ----------------------------------------------
if [[ "${SKIP_CLONE}" -eq 0 ]]; then
  if [[ ! -d "${PHASOMEIT_REPO}/.git" ]]; then
    echo "[setup] cloning PhasomeIt -> ${PHASOMEIT_REPO}"
    mkdir -p "$(dirname "${PHASOMEIT_REPO}")"
    git clone "${PHASOMEIT_GIT_URL}" "${PHASOMEIT_REPO}"
    chmod +x "${PHASOMEIT_REPO}/bossref" 2>/dev/null || true
  else
    echo "[setup] PhasomeIt already cloned at ${PHASOMEIT_REPO}"
  fi
fi

# ---------- 3. conda envs ----------------------------------------------------
if [[ "${SKIP_ENVS}" -eq 0 ]]; then
  if [[ -n "${CONDA_MODULE}" ]]; then
    module purge 2>/dev/null || true
    module load "${CONDA_MODULE}" 2>/dev/null || {
      echo "[setup] WARNING: 'module load ${CONDA_MODULE}' failed. If you're not on CSF, set CONDA_MODULE='' in config.sh."
    }
  fi
  if ! command -v conda >/dev/null; then
    echo "[setup] ERROR: conda not on PATH. Fix CONDA_MODULE in config.sh or load conda before running setup."
    exit 1
  fi
  source "$(conda info --base)/etc/profile.d/conda.sh"

  for spec in \
      "${ENV_PROKKA}:envs/prokka.yml" \
      "${ENV_PHASOMEIT}:envs/phasome.yml" \
      "${ENV_EGGNOG}:envs/eggnog.yml"; do
    name="${spec%%:*}"; yml="${spec##*:}"
    if [[ -d "${CONDA_ENVS_DIR}/${name}" ]]; then
      echo "[setup] env ${name} already exists at ${CONDA_ENVS_DIR}/${name}"
    else
      echo "[setup] creating env ${name} from ${yml}"
      conda env create -n "${name}" -f "${REPO_ROOT}/${yml}"
    fi
  done
fi

# ---------- 4. ensure dirs ---------------------------------------------------
mkdir -p "${LOGS_DIR}" "${INPUTS_DIR}" "${OUTPUTS_DIR}"

echo
echo "[setup] done."
echo "        next: bash check_config.sh"
echo "              bash run_all.sh"

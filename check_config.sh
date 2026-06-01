#!/usr/bin/env bash
# check_config.sh — verify config.sh resolves to valid targets.
# Distinguishes:
#   INPUT  (must exist before pipeline runs, except those created by stage 00)
#   OUTPUT (parent dir must be writable; the dir itself is created at runtime)
#   TOOL   (env or external repo)
#
# usage: bash check_config.sh

set -u

REPO_ROOT="$(cd "$(dirname "$0")" && pwd)"
if [[ ! -f "${REPO_ROOT}/config.sh" ]]; then
  echo "config.sh not found. Run: cp config.sh.example config.sh"
  exit 2
fi
# shellcheck disable=SC1091
source "${REPO_ROOT}/config.sh"

red(){    printf "\033[31m%s\033[0m\n" "$*"; }
green(){  printf "\033[32m%s\033[0m\n" "$*"; }
yellow(){ printf "\033[33m%s\033[0m\n" "$*"; }

n_err=0; n_warn=0

req(){    # required path must exist
  local k="$1" v="$2"
  if [[ -e "$v" ]]; then green "  OK     $k -> $v"
  else                   red   "  MISS   $k -> $v"; ((n_err++)); fi
}
opt(){    # optional path, warn if missing
  local k="$1" v="$2"
  if [[ -e "$v" ]]; then green  "  OK     $k -> $v"
  else                   yellow "  later  $k -> $v"; fi
}
parent(){ # parent dir must be writable (output target)
  local k="$1" v="$2"
  local p; p="$(dirname "$v")"
  if mkdir -p "$p" 2>/dev/null && [[ -w "$p" ]]; then
    green "  OK     $k -> $v (parent writable)"
  else
    red "  UNWR   $k -> $v"; ((n_err++))
  fi
}

echo "=== TOOLS ============================================================="
if [[ -n "${CONDA_MODULE}" ]]; then
  if module load "${CONDA_MODULE}" 2>/dev/null; then
    green "  OK     CONDA_MODULE=${CONDA_MODULE}"
  else
    yellow "  warn   CONDA_MODULE=${CONDA_MODULE} (module load failed; OK if conda already on PATH)"
    ((n_warn++))
  fi
fi
for e in "${ENV_PROKKA}" "${ENV_PHASOMEIT}" "${ENV_EGGNOG}"; do
  req "env:${e}" "${CONDA_ENVS_DIR}/${e}"
done
req PHASOMEIT_REPO "${PHASOMEIT_REPO}"

echo
echo "=== STANDARDS / INPUTS that must exist before submission =============="
req ARCHAEA_GENERA_FILE  "${ARCHAEA_GENERA_FILE}"
if [[ -n "${CONTROLS_SPECIES_TSV}" ]]; then
  req CONTROLS_SPECIES_TSV "${CONTROLS_SPECIES_TSV}"
fi

echo
echo "=== INPUTS that stage 00_download creates ============================="
opt GENOMES_ARCHAEA_DIR  "${GENOMES_ARCHAEA_DIR}"
opt GENOMES_CONTROLS_DIR "${GENOMES_CONTROLS_DIR}"
opt ATB_METADATA         "${ATB_METADATA}"
opt ARCHAEA_FASTA_LIST   "${ARCHAEA_FASTA_LIST}"
opt CONTROLS_FASTA_LIST  "${CONTROLS_FASTA_LIST}"

echo
echo "=== OUTPUT dirs (created at runtime, parents must be writable) ========"
for v in PROKKA_OUT PROKKA_CTRL_OUT SPECIES_CALLS_DIR GFF_LISTS_DIR \
         PHASOME_ARCHAEA_BASE PHASOME_CONTROLS_BASE EXTRACTION_DIR \
         EGGNOG_DIR EGGNOG_DB_DIR CLEANED_DIR QUANT_DIR ANALYSIS_DIR \
         FIGURES_DIR LOGS_DIR; do
  parent "${v}" "${!v}"
done

echo
echo "======================================================================="
if [[ $n_err -gt 0 ]]; then
  red "  ${n_err} problem(s) found. Fix them before submitting."
  exit 1
else
  green "  all green ($n_warn warning(s)). Ready to submit."
fi

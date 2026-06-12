#!/bin/bash
# config.sh - paths and env names for the pipeline.
#
# Every .job script in jobs/ sources this file. Edit the four lines
# in the EDIT-ME section below. Leave the rest alone unless you know
# what you are doing.
#
# This file is meant to live in the repo. Commit your edits; do not
# put secrets here.

# =============================================================================
# EDIT-ME
# =============================================================================
# Where the pipeline reads inputs from and writes outputs to.
# Use a path that is NOT inside the repo. On Manchester CSF this should be
# your scratch area, because outputs run into tens of GB.
SCRATCH_ROOT="${HOME}/scratch/archaea-pv-analysis"

# Conda envs the pipeline uses. Create these once (see README). If you
# already have envs with the right tools under different names, change
# these to point at them.
ENV_PROKKA="prokka_env"
ENV_PHASOMEIT="phasome_env"
ENV_EGGNOG="eggnog_env"

# The `module load` string for conda. On Manchester CSF the line below is
# correct. On a workstation where conda is already on PATH, set this to "".
CONDA_MODULE="apps/binapps/conda/miniforge3/25.9.1"

# SLURM partition. Manchester CSF default is multicore. Adjust for your cluster.
SLURM_PARTITION="multicore"

# Minimum strains per species used in the dataset filter (stage 01).
# A genus survives if any of its species has at least this many strains.
MIN_PER_SPECIES=5

# =============================================================================
# DERIVED PATHS - DO NOT EDIT BELOW THIS LINE
# =============================================================================
REPO_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

# Inputs and outputs live on scratch, NOT in the repo.
INPUTS_DIR="${SCRATCH_ROOT}/INPUTS"
OUTPUTS_DIR="${SCRATCH_ROOT}/OUTPUTS"
LOGS_DIR="${OUTPUTS_DIR}/logs"

# PhasomeIt source lives next to the data on scratch, not in the repo.
PHASOMEIT_REPO="${SCRATCH_ROOT}/external/PhasomeIt"
PHASOMEIT_GIT_URL="https://github.com/JackAidley/PhasomeIt.git"

# Input data layout (written by stage 00, read by later stages).
GENOMES_ARCHAEA_DIR="${INPUTS_DIR}/archaea/fasta"
GENOMES_CONTROLS_DIR="${INPUTS_DIR}/controls/fasta"
ATB_METADATA="${INPUTS_DIR}/archaea/atb_archaea_metadata.tsv"
ATB_TARBALL_CACHE="${INPUTS_DIR}/_tar_cache"
ARCHAEA_FASTA_LIST="${INPUTS_DIR}/archaea/fasta_list.txt"
CONTROLS_FASTA_LIST="${INPUTS_DIR}/controls/fasta_list.txt"
CONTROLS_GENUS_MAP="${INPUTS_DIR}/controls/controls_species.tsv"

# Standards (shipped with the repo, user-editable templates).
ARCHAEA_GENERA_FILE="${REPO_ROOT}/standards/default_archaea_genera.txt"
CONTROLS_SPECIES_TSV="${REPO_ROOT}/standards/default_controls_species.tsv"

# Output layout (one folder per stage).
PROKKA_OUT="${OUTPUTS_DIR}/01_prokka/archaea"
PROKKA_CTRL_OUT="${OUTPUTS_DIR}/01_prokka/controls"
SPECIES_CALLS_DIR="${OUTPUTS_DIR}/02_dataset/species_calls"
GFF_LISTS_DIR="${OUTPUTS_DIR}/02_dataset/genus_gff_lists"
PHASOME_ARCHAEA_BASE="${OUTPUTS_DIR}/03_phasomeit/archaea"
PHASOME_CONTROLS_BASE="${OUTPUTS_DIR}/03_phasomeit/controls"
EXTRACTION_DIR="${OUTPUTS_DIR}/04_extraction"
EGGNOG_DIR="${OUTPUTS_DIR}/05_eggnog"
EGGNOG_DB_DIR="${EGGNOG_DIR}/db"
CLEANED_DIR="${OUTPUTS_DIR}/06_cleaned"
QUANT_DIR="${OUTPUTS_DIR}/07_quantification"
ANALYSIS_DIR="${OUTPUTS_DIR}/08_analysis"
FIGURES_DIR="${OUTPUTS_DIR}/09_figures"

# Export so subshells (the python scripts the .job files launch) see them.
export REPO_ROOT SCRATCH_ROOT \
       CONDA_MODULE ENV_PROKKA ENV_PHASOMEIT ENV_EGGNOG \
       PHASOMEIT_REPO PHASOMEIT_GIT_URL \
       ARCHAEA_GENERA_FILE CONTROLS_SPECIES_TSV MIN_PER_SPECIES \
       INPUTS_DIR OUTPUTS_DIR LOGS_DIR \
       GENOMES_ARCHAEA_DIR GENOMES_CONTROLS_DIR \
       ATB_METADATA ATB_TARBALL_CACHE \
       ARCHAEA_FASTA_LIST CONTROLS_FASTA_LIST CONTROLS_GENUS_MAP \
       PROKKA_OUT PROKKA_CTRL_OUT \
       SPECIES_CALLS_DIR GFF_LISTS_DIR \
       PHASOME_ARCHAEA_BASE PHASOME_CONTROLS_BASE \
       EXTRACTION_DIR EGGNOG_DIR EGGNOG_DB_DIR \
       CLEANED_DIR QUANT_DIR ANALYSIS_DIR FIGURES_DIR \
       SLURM_PARTITION

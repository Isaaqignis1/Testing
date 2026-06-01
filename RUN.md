# Running the pipeline

Every stage submits to SLURM. Nothing here is intended to run on a CSF login
node. All `.job` scripts source `config.sh` from the repo root, so edit
`config.sh` once and the same paths are used everywhere.

The order below reflects the data dependencies. If a step's output already
exists the corresponding job script will detect it and skip work where it can.

## 0. one-time setup

```bash
cp config.sh.example config.sh
# edit config.sh: env names, scratch paths, AllTheBacteria locations

# create conda envs (off the login node — use an interactive job)
# prokka_env       : prokka
# phasome_env      : biopython + python 3.10+ (for PhasomeIt and parsers)
# eggnog           : eggnog-mapper v2

mkdir -p logs

# eggnog reference db (~50 GB download, run once)
sbatch eggnog/download_eggnog_db.job
```

Required environment variables (set in `config.sh`) that the new steps need
on top of the originals:

```bash
export ATB_METADATA="${HOME}/scratch/genomes/atb.archaea/metadata.tsv"
export SPECIES_CALLS_DIR="${HOME}/project1/species_calls"
export GFF_LISTS_DIR="${HOME}/project1/genus_gff_lists"
export PHASOME_ARCHAEA_BASE="${HOME}/scratch/phasome/genus_archaea"
export PHASOME_CONTROLS_BASE="${HOME}/project1/controls/phasome/genus"
export EXTRACTION_DIR="${HOME}/project1/Extraction"
export CLEANED_DIR="${HOME}/project1/Cleaned"
export QUANT_DIR="${HOME}/project1/quantification"
export ANALYSIS_DIR="${HOME}/project1/analysis"
export CONTROLS_SPECIES_TSV="${HOME}/project1/repo/standards/SpeciesCallsBacteria.tsv"
```

## 1. annotate genomes with Prokka

Independent: archaea and controls can run in parallel.

```bash
# adjust the SBATCH array range to match the number of fastas
sbatch annotation/prokka_archaea.job
sbatch annotation/prokka_controls.job
```

Outputs: `${PROKKA_OUT}/<sample_id>/<sample_id>.{gff,gbk,faa}` and the
matching files under `${PROKKA_CTRL_OUT}` for controls.

## 2. build the filtered species-call table

This is the step where the threshold has changed. The original write-up
applied `>= 5 strains per genus`; this build applies `>= 5 strains per
species`, then keeps a genus only if at least one of its species survives.

`prepare_dataset.job` also produces the per-genus `<Genus>_gffs.txt` files
that the PhasomeIt stage consumes.

```bash
sbatch dataset/prepare_dataset.job
```

Outputs in `${SPECIES_CALLS_DIR}`:

- `SpeciesCallsArchaea.csv`     — final per-sample call (replaces standards/SpeciesCallsArchaea.csv as the working copy)
- `species_count_summary.csv`   — pre- and post-threshold counts per species
- `genus_count_summary.csv`     — per-genus rollup

Outputs in `${GFF_LISTS_DIR}`:

- `<Genus>_gffs.txt`            — one path per kept sample, consumed by step 3
- `missing_samples.tsv`         — samples whose Prokka output was not found

Inspect `species_count_summary.csv` before continuing — this is where you
verify which species and genera dropped out under the new threshold.

## 3. stage PhasomeIt input folders

Copies each genus's Prokka `.gbk` files into a per-genus directory under
`${PHASOME_ARCHAEA_BASE}` and `${PHASOME_CONTROLS_BASE}`. Also writes a
`_genus_dirs.txt` list file in each base that the array job reads.

```bash
sbatch phasomeit/stage_phasomeit_inputs.job
```

## 4. run PhasomeIt

Submit twice, once per dataset. The job's `-a` array range must cover the
number of genera produced by step 3 (check `wc -l _genus_dirs.txt`).

```bash
# update -a in phasomeit/phasomeit_run.job if needed
sbatch --export=ALL,LIST=${PHASOME_ARCHAEA_BASE}/_genus_dirs.txt  phasomeit/phasomeit_run.job
sbatch --export=ALL,LIST=${PHASOME_CONTROLS_BASE}/_genus_dirs.txt phasomeit/phasomeit_run.job
```

Each genus directory will contain `summary_tracts/` (the PhasomeIt output)
after the run completes.

## 5. parse PhasomeIt HTML

Produces the extracted CSV tables consumed by quantification, eggnog
linking and the figure notebooks. Also produces the Prokka-function ranking
that backs the supplementary figure.

```bash
sbatch parsing/parse_phasomeit.job
```

Outputs (under `${EXTRACTION_DIR}`):

- `archaea/phasomeit_group_summary.csv`
- `archaea/phasomeit_annotated_functions.csv`
- `archaea/phasomeit_members.csv`
- `archaea/phasomeit_pairwise.csv`
- `archaea/phasomeit_genus_summary.csv`
- `archaea/prokka_function_ranking/pv_function_counts_by_genus.csv` (and total/genus summaries)
- the same set under `controls/`
- `phasomeit_tract_data.csv`              — combined archaeal + bacterial tract rows
- `phasomeit_full_group_summary.csv`      — domain-labelled group rows for the notebook

## 6. eggNOG annotation

Independent: archaea and bacteria runs are submitted separately. The
input-merger job concatenates Prokka .faa files into a single FASTA per
dataset; `archaea` is restricted to taxon 2157 to recover arCOG mappings.

```bash
sbatch eggnog/eggnog_input_merger.job
sbatch eggnog/run_eggnog_archaea.job
sbatch eggnog/run_eggnog_bacteria.job
```

When both emapper jobs are done:

```bash
sbatch eggnog/post_eggnog.job
```

Outputs (under `${CLEANED_DIR}`):

- `all_annotations.tsv`
- `archaea_annotations_clean.tsv`
- `per_sample_COG_summary.tsv`
- `per_genus_broad_role_normalised.tsv`
- `Mapped_functions/phasomeit_eggnog_linked.tsv`             — archaea
- `Mapped_functions/phasomeit_eggnog_linked_controls.tsv`    — controls

## 7. quantification

Per-genome PV counts, genome sizes, and the PV/Mb normalisation used in
Figures 1–2 of the report.

```bash
sbatch quantification/quantify.job
```

Outputs (under `${QUANT_DIR}`):

- `genome_lengths_archaea.csv`, `genome_lengths_controls.csv`
- `pv_per_fasta.csv`
- `pv_per_mb_per_genome.csv`
- `pv_per_mb_per_genus.csv`

## 8. method follow-up — tract-length variation

New analysis described in Section 4.4 of the report. Identifies gene
groups whose SSR tract length varies across strains of the same species
(high-confidence PV candidates) versus groups whose tract length is fixed
(SSR-bearing but no detectable switching).

```bash
sbatch analysis/tract_length_variation.job
```

Outputs (under `${ANALYSIS_DIR}`):

- `tract_length_variation_per_group.csv`
- `tract_length_variation_per_genus.csv`

## 9. figures

The notebooks under `figures/` read `${EXTRACTION_DIR}` and `${CLEANED_DIR}`
directly. They are written so that `PROJECT_ROOT` can be set with an
environment variable, but most users will edit the path constants in the
first cell. These notebooks are light enough to run interactively on a
data-transfer or visualisation node — they do not require the SLURM
compute partitions.

- `build_master_tables.ipynb` — joins extracted tracts, group summary and
  eggnog annotations; applies the four positional filters described in
  Section 2.5 of the report.
- `plot_pv_burden.ipynb`      — Figures 1 and 2.
- `plot_cog_functions.ipynb`  — Figures 3 and 4.

## stage-to-output summary

| stage | submit                                  | depends on                       |
|-------|-----------------------------------------|----------------------------------|
| 1a    | `annotation/prokka_archaea.job`         | fastas + ARCHAEA_FASTA_LIST      |
| 1b    | `annotation/prokka_controls.job`        | fastas + CONTROLS_FASTA_LIST     |
| 2     | `dataset/prepare_dataset.job`           | 1a                               |
| 3     | `phasomeit/stage_phasomeit_inputs.job`  | 2 + 1b                           |
| 4a    | `phasomeit/phasomeit_run.job` (arch)    | 3                                |
| 4b    | `phasomeit/phasomeit_run.job` (ctrls)   | 3                                |
| 5     | `parsing/parse_phasomeit.job`           | 4a + 4b                          |
| 6a    | `eggnog/eggnog_input_merger.job`        | 1a + 1b                          |
| 6b    | `eggnog/run_eggnog_archaea.job`         | 6a                               |
| 6c    | `eggnog/run_eggnog_bacteria.job`        | 6a                               |
| 6d    | `eggnog/post_eggnog.job`                | 5 + 6b + 6c                      |
| 7     | `quantification/quantify.job`           | 5                                |
| 8     | `analysis/tract_length_variation.job`   | 5                                |
| 9     | `figures/*.ipynb`                       | 5 + 6d + 7 + 8                   |

If you want to enforce ordering via SLURM dependencies, capture the job id
from each submission and pass `--dependency=afterok:<jid>` to the next.

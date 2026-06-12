# archaea-pv-analysis

SLURM pipeline that quantifies putative phase-variable (PPV) loci across
archaeal genomes, with bacterial controls. It downloads about 800
archaeal assemblies from AllTheBacteria, filters to a per-species
strain threshold, annotates with Prokka, runs PhasomeIt to find SSR
tracts in coding regions, runs eggNOG-mapper for functional annotation,
and produces the figures from the report.

The pipeline is 14 SLURM jobs you submit one at a time. There is no
orchestrator, no DAG runner, no marker system. Every `.job` file in
`jobs/` is a normal `sbatch` script with a clear header explaining its
inputs, outputs, and how to submit it.

## Repo layout

```
.
├── README.md             this file
├── LICENSE
├── config.sh             paths and conda env names (edit before first run)
├── submit.sh             thin wrapper around sbatch (optional, 30 lines)
│
├── jobs/                 one numbered .job per stage, 00 to 13
│   ├── 00_download.job
│   ├── 01_dataset.job
│   ├── 02_prokka.job
│   ├── 03_stage_phasomeit.job
│   ├── 04_phasomeit.job
│   ├── 05_parse_phasomeit.job
│   ├── 06_eggnog_db.job
│   ├── 07_eggnog_merge.job
│   ├── 08_eggnog_archaea.job
│   ├── 09_eggnog_bacteria.job
│   ├── 10_post_eggnog.job
│   ├── 11_quantify.job
│   ├── 12_tract_variation.job
│   └── 13_figures.job
│
├── scripts/              python called by the jobs
├── envs/                 conda yamls (prokka, phasome, eggnog)
└── standards/            archaeal genus list, controls, contaminants, colours
```

Three things to know about the layout:

1. The repo only contains code and small text files. Data, derived
   outputs, conda envs, and the PhasomeIt source all live elsewhere
   (under `$HOME/scratch/archaea-pv-analysis` by default; set in
   `config.sh`).
2. `submit.sh` is optional. It is 30 lines and you can read it. It
   sources `config.sh`, makes the logs directory, and calls `sbatch`
   with `--export=ALL,REPO_ROOT=$PWD` and the right `-o`/`-e` paths.
   Every example in this README also shows the raw `sbatch` form.
3. Each `.job` file has a comment block at the top that tells you
   what it needs, what it produces, how long it takes, and how to
   submit it. Read the file before submitting it.

## Requirements

- A SLURM cluster.
- conda (miniforge or miniconda). On Manchester CSF this is loaded with
  `module load apps/binapps/conda/miniforge3/25.9.1`.
- Roughly 200 GB of scratch (the eggNOG reference db is ~50 GB, the
  Prokka outputs for ~400 archaeal samples are another big chunk).

## One-time setup

### 1. Clone the repo

```bash
git clone <this-repo> archaea-pv-analysis
cd archaea-pv-analysis
```

### 2. Edit config.sh

Open `config.sh`. Edit only the lines in the `EDIT-ME` block at the
top. On Manchester CSF the only line you need to change is usually
`SCRATCH_ROOT`. Everything below the `DO NOT EDIT BELOW THIS LINE`
marker is derived from those few values.

### 3. Create the three conda envs

This can take 5 to 20 minutes per env on a fresh install:

```bash
module load apps/binapps/conda/miniforge3/25.9.1   # CSF; skip if conda is on PATH already
conda env create -f envs/prokka.yml
conda env create -f envs/phasome.yml
conda env create -f envs/eggnog.yml
```

If you already have envs with the right tools under different names,
point `ENV_PROKKA`, `ENV_PHASOMEIT`, `ENV_EGGNOG` in `config.sh` at them
instead.

### 4. Clone PhasomeIt

```bash
source config.sh
git clone "$PHASOMEIT_GIT_URL" "$PHASOMEIT_REPO"
chmod +x "$PHASOMEIT_REPO/bossref" 2>/dev/null || true
```

### 5. (Optional) verify the envs work

```bash
conda run -n "$ENV_PROKKA"    prokka --version
conda run -n "$ENV_PHASOMEIT" python -c 'import Bio, pandas, numpy, matplotlib; print("ok")'
conda run -n "$ENV_EGGNOG"    emapper.py --version
```

### 6. (Optional) supply controls

If you want bacterial controls in your run, drop their FASTAs into
`$GENOMES_CONTROLS_DIR` and write a TSV mapping sample ID to species
at `$CONTROLS_GENUS_MAP`. There is a template at
`standards/default_controls_species.tsv` you can adapt. Stage 03 and
stage 09 detect controls automatically and skip them if absent.

## Running the pipeline

The recommended workflow is to submit one job, wait for it to finish,
look at its output, then submit the next. Both forms are shown below
for stage 00 and they do exactly the same thing. For later stages
only the wrapper form is shown.

### Stage 00: download archaeal genomes

```bash
bash submit.sh jobs/00_download.job
```

Or, raw sbatch:

```bash
source config.sh
mkdir -p "$LOGS_DIR"
sbatch --export=ALL,REPO_ROOT=$PWD \
       -o "$LOGS_DIR/%x_%j.out" -e "$LOGS_DIR/%x_%j.err" \
       jobs/00_download.job
```

While it runs:

```bash
squeue -u $USER
```

When it finishes, verify:

```bash
ls "$GENOMES_ARCHAEA_DIR" | wc -l    # should be about 815 fastas
wc -l "$ATB_METADATA"                # should be 819 lines (818 samples + header)
tail -30 "$LOGS_DIR/00_download_<jobid>.out"
```

### Stage 01: dataset filter

```bash
bash submit.sh jobs/01_dataset.job
```

When it finishes, look at the filter summary before paying for Prokka:

```bash
cat "$SPECIES_CALLS_DIR/species_count_summary.csv"
cat "$SPECIES_CALLS_DIR/genus_count_summary.csv"
wc -l "$INPUTS_DIR/archaea/filtered_fasta_list.txt"
```

To use a different threshold, change `MIN_PER_SPECIES` in `config.sh`
and resubmit stage 01.

### Stage 02: Prokka annotation (array job)

This is a SLURM array, one task per filtered sample. Submit it twice.
The first call is a dry run that prints the array range to use:

```bash
bash submit.sh jobs/02_prokka.job
# prints something like:
#   filtered list has 414 samples
#   re-submit with:
#     bash submit.sh jobs/02_prokka.job --array=1-414%50

bash submit.sh jobs/02_prokka.job --array=1-414%50
```

The `%50` caps how many tasks run at once. Tune for your cluster's
queue limits. Each task can take up to a day on CSF.

Verify:

```bash
ls "$PROKKA_OUT" | wc -l
sacct -X --format=JobID,JobName,State,ExitCode -j <array_jobid>
```

### Stage 03: stage PhasomeIt inputs

```bash
bash submit.sh jobs/03_stage_phasomeit.job
```

Verify a few `.gbk` files exist under
`$PHASOME_ARCHAEA_BASE/<some_genus>/`.

### Stage 04: PhasomeIt (array job)

Same pattern as stage 02. Dry-run first, then re-submit with the array
range:

```bash
bash submit.sh jobs/04_phasomeit.job
bash submit.sh jobs/04_phasomeit.job --array=1-N%9
```

If you also have controls staged at `$PHASOME_CONTROLS_BASE`, submit a
second array for them by overriding `LIST`:

```bash
LIST=$PHASOME_CONTROLS_BASE/_genus_dirs.txt
sbatch --export=ALL,REPO_ROOT=$PWD,LIST=$LIST \
       -o "$LOGS_DIR/%x_%j.out" -e "$LOGS_DIR/%x_%j.err" \
       --array=1-M%9 \
       jobs/04_phasomeit.job
```

### Stages 05 through 13

Each is a single submission. Submit, wait, verify with the paths listed
at the top of each `.job` file.

```bash
bash submit.sh jobs/05_parse_phasomeit.job
bash submit.sh jobs/06_eggnog_db.job          # can run any time after setup
bash submit.sh jobs/07_eggnog_merge.job
bash submit.sh jobs/08_eggnog_archaea.job
bash submit.sh jobs/09_eggnog_bacteria.job    # skip if not using controls
bash submit.sh jobs/10_post_eggnog.job
bash submit.sh jobs/11_quantify.job
bash submit.sh jobs/12_tract_variation.job
bash submit.sh jobs/13_figures.job
```

Stage 06 has no dependencies. You can submit it in parallel with the
early stages to overlap the ~8 hours of database download with Prokka.

## Stage reference

| Stage | Job file                  | Depends on  | Wall time      |
|-------|---------------------------|-------------|----------------|
| 00    | 00_download.job           | nothing     | ~30 min        |
| 01    | 01_dataset.job            | 00          | ~10 min        |
| 02    | 02_prokka.job             | 01          | up to 24h/task |
| 03    | 03_stage_phasomeit.job    | 02          | ~1h            |
| 04    | 04_phasomeit.job          | 03          | up to 6h/task  |
| 05    | 05_parse_phasomeit.job    | 04          | ~6h            |
| 06    | 06_eggnog_db.job          | nothing     | ~8h            |
| 07    | 07_eggnog_merge.job       | 02          | ~4h            |
| 08    | 08_eggnog_archaea.job     | 06, 07      | up to 72h      |
| 09    | 09_eggnog_bacteria.job    | 06, 07      | up to 36h      |
| 10    | 10_post_eggnog.job        | 05, 08, 09  | ~2h            |
| 11    | 11_quantify.job           | 05          | ~1h            |
| 12    | 12_tract_variation.job    | 05          | ~1h            |
| 13    | 13_figures.job            | 10, 11, 12  | ~2h            |

## Submitting everything chained (optional)

If you want a single command that queues every stage and lets SLURM
handle ordering, the pattern is `--dependency=afterok:<prev_jobid>`
between consecutive `sbatch` calls. Save the snippet below as
`submit_all.sh` (it is not committed because the array sizes are
specific to your filtered cohort):

```bash
#!/bin/bash
set -euo pipefail
source config.sh
mkdir -p "$LOGS_DIR"

LAST=""
submit_after () {
  local job="$1"; shift
  local args=(--parsable
              --export="ALL,REPO_ROOT=$PWD"
              -o "$LOGS_DIR/%x_%j.out"
              -e "$LOGS_DIR/%x_%j.err")
  [ -n "$LAST" ] && args+=(--dependency=afterok:"$LAST")
  args+=("$@" "$job")
  LAST=$(sbatch "${args[@]}")
  echo "  $job  jobid=$LAST"
}

submit_after jobs/00_download.job
submit_after jobs/01_dataset.job
submit_after jobs/02_prokka.job              --array=1-414%50
submit_after jobs/03_stage_phasomeit.job
submit_after jobs/04_phasomeit.job           --array=1-50%9
submit_after jobs/05_parse_phasomeit.job
submit_after jobs/06_eggnog_db.job
submit_after jobs/07_eggnog_merge.job
submit_after jobs/08_eggnog_archaea.job
submit_after jobs/09_eggnog_bacteria.job
submit_after jobs/10_post_eggnog.job
submit_after jobs/11_quantify.job
submit_after jobs/12_tract_variation.job
submit_after jobs/13_figures.job
```

Set the array sizes (`1-414%50` and `1-50%9`) to the values stage 02's
and stage 04's dry runs print for your filtered cohort.

## What you might want to edit

Three files change what gets analysed; everything else is derived:

| File                                       | Edit to change                                   |
|--------------------------------------------|--------------------------------------------------|
| `config.sh`                                | paths, conda env names, threshold, partition    |
| `standards/default_archaea_genera.txt`     | which archaeal genera to keep                    |
| `standards/default_controls_species.tsv`   | which bacterial controls and minimum strain count|

The defaults target the 18 archaeal genera from the original report,
use Brucella + E. coli + Campylobacter as controls (8 strains per
species), apply a "5 or more strains per species" retention rule, and
write data to `$HOME/scratch/archaea-pv-analysis/`.

## Threshold note

The original report used a per-genus "5 or more strains" threshold. This
pipeline applies the rule at the species level instead: a genus is kept
if at least one of its species has `MIN_PER_SPECIES` or more strains.
Change `MIN_PER_SPECIES` in `config.sh` to override.

After stage 01 finishes, look at
`$SPECIES_CALLS_DIR/species_count_summary.csv` to see which species and
genera dropped out before committing compute to Prokka and PhasomeIt.

## Outputs you care about

Once stage 13 has run, the things to look at are under `$FIGURES_DIR`:

```
master_strains.csv          one row per strain
master_groups.csv           one row per (genus, gene group)
fig1_*.{pdf,png}            PV per Mb scatter
fig2_*.{pdf,png}            PV per Mb boxplot
fig3_*.{pdf,png}            COG composition
fig4_*.{pdf,png}            COG heatmap
supp_top_prokka_functions.{pdf,png}
```

The CSVs are the underlying tables for every figure and are useful in
their own right.

## When something breaks

Every job's stdout and stderr go to
`$LOGS_DIR/<stage>_<jobid>.{out,err}`. Start there.

Common failures:

- `ModuleNotFoundError: No module named '...'`. The env named in
  `config.sh` is missing a package. Verify with
  `conda run -n $ENV_PHASOMEIT python -c 'import Bio'`. Either fix the
  env (`conda env update -n $ENV_PHASOMEIT -f envs/phasome.yml`) or
  point `ENV_PHASOMEIT` at a working env.
- Job dies in 1 second with
  `slurm_script: ... config.sh: No such file or directory`. You
  submitted without `--export=ALL,REPO_ROOT=$PWD`. Use `submit.sh` or
  add the flag.
- `02_prokka.job` or `04_phasomeit.job` finishes in seconds with a
  usage message printed. You did not pass `--array=...`. Re-submit
  with the range it printed.
- `download_eggnog_data.py: command not found`. You used the wrong
  env. The eggNOG tools live in `$ENV_EGGNOG`, not `$ENV_PHASOMEIT`.

## Licence

MIT. See `LICENSE`.

# RUN.md — manual sbatch reference

For the impatient: `bash run_all.sh` chains everything below with
`--dependency=afterok`. The reference here is for re-running a single stage
or debugging.

Every job sources `${REPO_ROOT}/config.sh`, writes logs to
`OUTPUTS/logs/<stage>_<jobid>.out|err`, and is safe to re-run (existing
outputs are detected and skipped where possible).

## Stage map

| #  | Job file                              | Depends on | Outputs (under OUTPUTS/)         |
|----|---------------------------------------|------------|----------------------------------|
| 00 | `pipeline/00_download.job`            | —          | INPUTS/archaea, INPUTS/controls  |
| 01 | `pipeline/01_prokka.job`              | 00         | 01_prokka/{archaea,controls}/    |
| 02 | `pipeline/02_dataset.job`             | 01         | 02_dataset/                      |
| 03 | `pipeline/03_stage_phasomeit.job`     | 02         | 03_phasomeit/                    |
| 04 | `pipeline/04_phasomeit.job`           | 03         | 03_phasomeit/<genus>/summary_tracts/ |
| 05 | `pipeline/05_parse_phasomeit.job`     | 04         | 04_extraction/                   |
| 06 | `pipeline/06_eggnog_db.job`           | —          | 05_eggnog/db/ (one-off)          |
| 07 | `pipeline/07_eggnog_merge.job`        | 01         | 05_eggnog/input/*.faa            |
| 08 | `pipeline/08_eggnog_archaea.job`      | 06, 07     | 05_eggnog/output/archaea.*       |
| 09 | `pipeline/09_eggnog_bacteria.job`     | 06, 07     | 05_eggnog/output/bacteria.*      |
| 10 | `pipeline/10_post_eggnog.job`         | 05, 08, 09 | 06_cleaned/                      |
| 11 | `pipeline/11_quantify.job`            | 05         | 07_quantification/               |
| 12 | `pipeline/12_tract_variation.job`     | 05         | 08_analysis/                     |
| 13 | `pipeline/13_figures.job`             | 10, 11, 12 | 09_figures/                      |

## Stage 01 (Prokka) and stage 04 (PhasomeIt) are array jobs

Both auto-detect their array size on first invocation. Submit once without
`-a` and the script prints the correct `-a 1-N` command for re-submission.

```bash
# Stage 01 (Prokka)
sbatch pipeline/01_prokka.job              # prints "submit with -a 1-N"
sbatch -a 1-N%50 pipeline/01_prokka.job    # actually runs

# Stage 04 (PhasomeIt) — once per dataset
sbatch pipeline/04_phasomeit.job                                            # archaea, prints -a 1-N
sbatch -a 1-N%9 pipeline/04_phasomeit.job
sbatch --export=ALL,LIST=${PHASOME_CONTROLS_BASE}/_genus_dirs.txt pipeline/04_phasomeit.job
sbatch -a 1-M%9 --export=ALL,LIST=${PHASOME_CONTROLS_BASE}/_genus_dirs.txt pipeline/04_phasomeit.job
```

`run_all.sh` handles the array sizing automatically (it queries the list
file length and inserts the right `-a` range).

## Resuming after a failure

```bash
bash run_all.sh --from 05    # skips 00-04, resumes at 05
```

You can also re-submit a single stage by hand. The downstream stages will
not re-run unless you call `run_all.sh --from <next>` again, since
dependencies are wired through `run_all.sh`.

## Where to look when something breaks

| Symptom                                  | Look at                                                                 |
|------------------------------------------|--------------------------------------------------------------------------|
| Stage 00 fails partway through           | `OUTPUTS/logs/00_download_*.err` — usually a transient HTTP error; rerun |
| Some species missing from species_calls  | `OUTPUTS/02_dataset/species_calls/species_count_summary.csv`            |
| PhasomeIt genus task fails / hangs       | `OUTPUTS/logs/04_phasomeit_*.err` + `OUTPUTS/03_phasomeit/<genus>/`     |
| Eggnog out of memory                     | bump `-c 16` to `-c 32` in stage 08/09; rerun                            |
| Figures look wrong                       | rerun `bash run_all.sh --only 13` after editing `standards/colors.py`   |

## Notebook iteration

The `.py` figure scripts under `scripts/figures/` are mirrored as notebooks
under `scripts/figures/notebooks/`. The notebooks read the same CSVs and
can be opened on a transfer node or via a Jupyter server. They are NOT
re-submitted by `run_all.sh` — the `.py` versions are the canonical
producers of the deliverable figures.

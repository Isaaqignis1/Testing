# RUN.md — manual sbatch reference

For most uses: `bash run_all.sh` chains everything with `--dependency=afterok`
and injects the SLURM `--output`/`--error` paths into `${LOGS_DIR}`. This
file is for re-running a single stage by hand.

## Manual sbatch — important

The `.job` files do NOT specify `#SBATCH -o`/`-e`. If you submit a job
directly (without `run_all.sh`), pass `--output` and `--error` yourself or
SLURM will dump logs into your current directory:

```bash
sbatch \
    --output=${LOGS_DIR}/01_prokka_%A_%a.out \
    --error=${LOGS_DIR}/01_prokka_%A_%a.err \
    pipeline/01_prokka.job
```

Or just `cd "$LOGS_DIR"` before each `sbatch`.

## Stage map

| #  | Job file                              | Depends on  | Outputs (under OUTPUTS_DIR/)         |
|----|---------------------------------------|-------------|--------------------------------------|
| 00 | `pipeline/00_download.job`            | —           | INPUTS/{archaea,controls}/           |
| 01 | `pipeline/01_prokka.job`              | 00          | 01_prokka/{archaea,controls}/        |
| 02 | `pipeline/02_dataset.job`             | 01          | 02_dataset/                          |
| 03 | `pipeline/03_stage_phasomeit.job`     | 02          | 03_phasomeit/                        |
| 04 | `pipeline/04_phasomeit.job`           | 03          | 03_phasomeit/<genus>/summary_tracts/ |
| 05 | `pipeline/05_parse_phasomeit.job`     | 04          | 04_extraction/                       |
| 06 | `pipeline/06_eggnog_db.job`           | —           | 05_eggnog/db/ (one-off, ~50 GB)      |
| 07 | `pipeline/07_eggnog_merge.job`        | 01          | 05_eggnog/input/*.faa                |
| 08 | `pipeline/08_eggnog_archaea.job`      | 06, 07      | 05_eggnog/output/archaea.*           |
| 09 | `pipeline/09_eggnog_bacteria.job`     | 06, 07      | 05_eggnog/output/bacteria.*          |
| 10 | `pipeline/10_post_eggnog.job`         | 05, 08, 09  | 06_cleaned/                          |
| 11 | `pipeline/11_quantify.job`            | 05          | 07_quantification/                   |
| 12 | `pipeline/12_tract_variation.job`     | 05          | 08_analysis/                         |
| 13 | `pipeline/13_figures.job`             | 10, 11, 12  | 09_figures/                          |

## Array jobs (01, 04)

Both auto-detect their array size on first invocation. Submit once without
`-a` and the script prints the correct `-a 1-N` to use. `run_all.sh` does
this for you.

```bash
# stage 01 (Prokka) — combined list of all fastas
sbatch pipeline/01_prokka.job                     # prints "submit with -a 1-N"
sbatch -a 1-N%50 pipeline/01_prokka.job

# stage 04 (PhasomeIt) — once per dataset
sbatch pipeline/04_phasomeit.job                  # archaea, prints -a 1-N
sbatch -a 1-N%9 pipeline/04_phasomeit.job
sbatch --export=ALL,LIST=${PHASOME_CONTROLS_BASE}/_genus_dirs.txt pipeline/04_phasomeit.job
sbatch -a 1-M%9 --export=ALL,LIST=${PHASOME_CONTROLS_BASE}/_genus_dirs.txt pipeline/04_phasomeit.job
```

## Common failure → fix

| Symptom                                  | Look at                                                       |
|------------------------------------------|----------------------------------------------------------------|
| Stage 00 dies partway                    | `${LOGS_DIR}/00_download_*.err` — usually transient HTTP, rerun |
| Lots of species dropped after stage 02   | `${SPECIES_CALLS_DIR}/species_count_summary.csv`              |
| PhasomeIt task hangs / OOMs              | `${LOGS_DIR}/04_phasomeit_*.err`, then bump `-t` or `-c`       |
| Eggnog OOM                               | bump `-c 16` to `-c 32` in stage 08/09                         |
| Figures look wrong                       | rerun `bash run_all.sh --only 13` after editing `standards/colors.py` |
| Conda env missing on first run           | `bash setup.sh` again, or `conda env create -n NAME -f envs/NAME.yml` |

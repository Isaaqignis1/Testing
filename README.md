# archaea-pv-analysis

End-to-end SLURM pipeline that quantifies putative phase-variable (PPV)
loci across archaeal genomes, with bacterial controls. Runs from a fresh
clone with no manual data preparation — it downloads everything it needs,
filters to a per-species threshold, annotates, runs PhasomeIt + eggNOG,
quantifies PPV per Mb and per strain, and produces the report figures.

## Three commands to start

```bash
git clone <this-repo> Testing && cd Testing
bash setup.sh                  # writes config.sh, clones PhasomeIt, builds 3 conda envs
bash check_setup.sh            # confirms every path, env and tool works
bash run_all.sh                # submits all 14 stages, chained with --dependency=afterok
```

Resume / partial reruns:

```bash
bash run_all.sh --from 05      # resume from stage 05
bash run_all.sh --only 11      # run only the quantification stage
bash run_all.sh --dry-run      # print sbatch commands, do not submit
```

## Where things live

The repo itself stays small (text only). Heavy stuff goes elsewhere:

| Where             | What's there                                            | Default location                                            |
|-------------------|----------------------------------------------------------|--------------------------------------------------------------|
| Repo              | scripts, jobs, standards, env yamls, README             | wherever you cloned it                                       |
| Scratch (data)    | ATB downloads, derived data, figures, SLURM logs        | `${HOME}/scratch/archaea-pv-analysis/{INPUTS,OUTPUTS}`       |
| Scratch (tools)   | the PhasomeIt clone                                     | `${HOME}/scratch/archaea-pv-analysis/external/PhasomeIt`     |
| Conda             | the three named environments (`prokka_env`, `phasome_env`, `eggnog_env`) | wherever conda keeps your envs                |

Everything is configurable in `config.sh` — most users edit nothing.

## What you might want to edit

Three things change what gets analysed:

| File                                       | Edit to...                                  |
|--------------------------------------------|----------------------------------------------|
| `standards/default_archaea_genera.txt`     | analyse a different archaeal genus set       |
| `standards/default_controls_species.tsv`   | use different bacterial controls / strain counts |
| `config.sh`                                | change paths, threshold, cluster module      |

The default config:
- targets the 18 archaeal genera from the original report,
- pulls Brucella + E. coli + Campylobacter controls (≥8 strains per species),
- applies a `>= 5 strains per species` retention rule (genus survives if any species clears),
- writes downloads to `${HOME}/scratch/archaea-pv-analysis/INPUTS`,
- writes outputs to `${HOME}/scratch/archaea-pv-analysis/OUTPUTS`,
- expects conda envs named `prokka_env`, `phasome_env`, `eggnog_env`.

## Repo layout

```
.
├── README.md                       this file
├── RUN.md                          per-stage sbatch reference
├── LICENSE
├── config.sh.example               heavily documented; copy to config.sh
├── setup.sh                        one-off env + PhasomeIt setup
├── check_setup.sh                  comprehensive setup verifier
├── run_all.sh                      chains all stages with afterok
│
├── pipeline/                       one numbered .job per stage (00-13)
├── scripts/                        python modules called by the jobs
├── standards/                      genus list, controls list, colours, contaminants
└── envs/                           conda yamls for the 3 envs
```

Pipeline stages:

| #  | Stage                  | What it does                                          |
|----|------------------------|-------------------------------------------------------|
| 00 | download               | ATB tarballs → INPUTS                                 |
| 01 | prokka                 | annotate every assembly                               |
| 02 | dataset                | per-species filter + per-genus gff lists              |
| 03 | stage_phasomeit        | copy .gbk files into per-genus folders                |
| 04 | phasomeit              | run PhasomeIt per genus                               |
| 05 | parse_phasomeit        | HTML → CSV; PPV groups, tracts, members, pairwise     |
| 06 | eggnog_db              | one-off ~50 GB reference download                     |
| 07 | eggnog_merge           | concat per-genome .faa into per-domain proteomes      |
| 08 | eggnog_archaea         | emapper, tax_scope Archaea (arCOG)                    |
| 09 | eggnog_bacteria        | emapper, default tax scope                            |
| 10 | post_eggnog            | parse + link PV members to COG annotations            |
| 11 | quantify               | PV per genome, per Mb, per genus                      |
| 12 | tract_variation        | within-species tract-length variance (Section 4.4)    |
| 13 | figures                | Figures 1-4 + supplementary as PDF + PNG              |

## Threshold note

The original report used a per-genus `≥5` strain threshold. This pipeline
applies `≥5 strains per species` (a genus is kept if any of its species
clears the threshold). Change `MIN_PER_SPECIES` in `config.sh` to override.

After stage 02 finishes, inspect
`${OUTPUTS_DIR}/02_dataset/species_calls/species_count_summary.csv`
to see which species and genera dropped out before committing compute to
stages 03+.

## When something breaks

`bash check_setup.sh` reports the most common problems with a pass/fail per
check. For runtime failures, every job's stdout/stderr is in
`${LOGS_DIR}/<stage>_<jobid>.{out,err}`.

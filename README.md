# archaea-pv-analysis

End-to-end SLURM pipeline that quantifies putative phase-variable (PPV) loci
across archaeal genomes, with bacterial controls for reference. The pipeline
runs from a fresh clone with no external inputs — it downloads everything it
needs, organises outputs neatly inside the repo, and produces the figures
from the report.

## Quick start

```bash
# on CSF (or any SLURM cluster)
git clone <this repo> Testing
cd Testing

bash setup.sh                  # creates config.sh, clones PhasomeIt, builds 3 conda envs
bash check_config.sh           # verifies every path/env resolves
bash run_all.sh                # submits all stages with --dependency=afterok
```

Resume / partial reruns:

```bash
bash run_all.sh --from 05      # resume from stage 05
bash run_all.sh --only 11      # run only the quantification stage
bash run_all.sh --dry-run      # print the sbatch commands, do not submit
```

## What you might want to edit

By default the pipeline analyses the 18 archaeal genera from the original
report against Brucella / E. coli / Campylobacter controls (8 strains per
species), with a `>=5 strains per species` retention threshold.

To change the targets, edit one of three files — nothing else:

| File                                           | Edit to...                              |
|------------------------------------------------|------------------------------------------|
| `standards/default_archaea_genera.txt`         | analyse different archaeal genera        |
| `standards/default_controls_species.tsv`       | use different bacterial controls          |
| `config.sh`                                    | change paths, threshold, cluster module   |

All paths in `config.sh` default to locations inside the repo — `INPUTS/`
for downloaded data and `OUTPUTS/` for everything derived. Override
`INPUTS_DIR` or `OUTPUTS_DIR` if you'd rather keep raw data on scratch.

## Layout

```
.
├── README.md                       ← this file
├── RUN.md                          step-by-step sbatch reference
├── LICENSE
├── config.sh.example               heavily documented, copy to config.sh
├── setup.sh                        one-off env setup + PhasomeIt clone
├── check_config.sh                 sanity checker
├── run_all.sh                      submit chain, with afterok dependencies
│
├── pipeline/                       one numbered .job per stage
│   ├── 00_download.job             ATB metadata + fastas
│   ├── 01_prokka.job               array, one task per fasta
│   ├── 02_dataset.job              per-species filtering + gff lists
│   ├── 03_stage_phasomeit.job      stage per-genus .gbk folders
│   ├── 04_phasomeit.job            array, one task per genus
│   ├── 05_parse_phasomeit.job      HTML -> CSV tables
│   ├── 06_eggnog_db.job            one-off ~50 GB download
│   ├── 07_eggnog_merge.job         concat per-genome .faa
│   ├── 08_eggnog_archaea.job       emapper, tax_scope Archaea
│   ├── 09_eggnog_bacteria.job      emapper, default
│   ├── 10_post_eggnog.job          parse + link to PV members
│   ├── 11_quantify.job             PV per genome, per Mb, per genus
│   ├── 12_tract_variation.job      within-species tract-length variance
│   └── 13_figures.job              Figures 1-4 + supplementary
│
├── scripts/                        python modules called by the jobs
│   ├── dataset/                    download_atb, build_species_calls, build_genus_gff_lists
│   ├── phasomeit/                  setup_archaea, setup_controls
│   ├── parsing/                    parse_groups, parse_tracts, summarise_runs, build_full_group_summary, summarise_prokka_functions
│   ├── eggnog/                     merge_proteomes, parse_eggnog, link_eggnog
│   ├── quantification/             genome_lengths, pv_per_fasta, pv_per_mb
│   ├── analysis/                   tract_length_variation
│   └── figures/                    build_master_tables, plot_pv_burden, plot_cog_functions, plot_prokka_top_functions
│                                   notebooks/* mirror the .py figures for interactive tweaking
│
├── standards/                      canonical reference data (committed)
│   ├── colors.py                   per-genus colours + plot order
│   ├── default_archaea_genera.txt  18 genera analysed in the report
│   ├── default_controls_species.tsv default control list (>=8/species)
│   └── contaminants.txt            scientific names to drop
│
├── envs/                           conda yamls (used by setup.sh)
│   ├── prokka.yml
│   ├── phasome.yml
│   └── eggnog.yml
│
├── external/                       PhasomeIt is cloned here (gitignored)
├── INPUTS/                         downloaded source data (gitignored)
│   ├── archaea/                    ATB fastas + metadata
│   └── controls/                   control fastas + species map
└── OUTPUTS/                        all derived data (gitignored)
    ├── 01_prokka/{archaea,controls}/
    ├── 02_dataset/{species_calls,genus_gff_lists}/
    ├── 03_phasomeit/{archaea,controls}/
    ├── 04_extraction/
    ├── 05_eggnog/{input,output,db}/
    ├── 06_cleaned/
    ├── 07_quantification/
    ├── 08_analysis/
    ├── 09_figures/
    └── logs/
```

## Threshold note

The original write-up used a per-genus `>=5` strain threshold. This pipeline
applies a `>=5 strains per species` rule (a genus is kept if any of its
species clears the threshold). The threshold lives in
`MIN_PER_SPECIES` in `config.sh`. Setting it to `1` and re-running stage 02
recovers the original per-genus behaviour.

After stage 02, inspect `OUTPUTS/02_dataset/species_calls/species_count_summary.csv`
to see which species and genera dropped out under the threshold before
committing compute to stages 03+.

## Reproducing the report figures

`pipeline/13_figures.job` produces Figures 1–4 as PDFs and PNGs in
`OUTPUTS/09_figures/`. The figures depend only on output CSVs from stages
07, 10, 11 and 12, so you can rerun figures alone (`run_all.sh --only 13`)
after tweaking colours or layout in `standards/colors.py`.

If you'd rather iterate interactively, use the notebook mirrors under
`scripts/figures/notebooks/`. They read the same CSVs as the scripts.

## Citation

This pipeline implements the analysis described in Isaaqignis (2026),
*The cultivable archaeal phasome: a comparative-genomics survey of
simple-sequence-repeat-mediated phase variation across 18 archaeal genera*
(MSc dissertation, University of Manchester).

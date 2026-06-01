# archaeal-phasome

Scripts for identifying putative phase-variable (PPV) loci in archaeal
genomes, with three bacterial control genera as reference points.

The repo contains the SLURM jobs and Python scripts that run Prokka,
PhasomeIt and eggNOG-mapper across the AllTheBacteria archaeal set,
parse PhasomeIt HTML output, count and normalise PPV loci, join
PhasomeIt PV genes to eggNOG annotations, run a within-species
tract-length variation follow-up, and the three Jupyter notebooks that
produce the report figures.

**Threshold change from the original write-up.** The published version
applied a minimum-strain threshold of `>= 5 per genus`. This repo applies
`>= 5 per species`; a genus is kept only if at least one of its species
clears the threshold. The species-call table is regenerated from
AllTheBacteria metadata by `dataset/build_species_calls.py`, not from a
checked-in static CSV.

See [RUN.md](RUN.md) for the step-by-step run order on a CSF-style SLURM
cluster.

## Layout

```
repo/
├── README.md
├── RUN.md                          end-to-end run order with sbatch commands
├── LICENSE
├── config.sh.example               copy to config.sh and edit
├── .gitignore
│
├── standards/
│   ├── colors.py                   per-genus colour map and plot ordering
│   ├── SpeciesCallsArchaea.csv     reference copy of the per-species threshold output
│   └── SpeciesCallsBacteria.tsv    controls species map (input, not generated)
│
├── dataset/                        NEW — replaces the static species-call CSV
│   ├── build_species_calls.py      filter AllTheBacteria metadata, apply per-species >=5 threshold
│   ├── build_genus_gff_lists.py    emit <Genus>_gffs.txt lists from species calls + prokka outputs
│   └── prepare_dataset.job         SLURM wrapper for both
│
├── annotation/
│   ├── prokka_archaea.job          prokka, kingdom Archaea
│   └── prokka_controls.job         prokka, kingdom Bacteria
│
├── phasomeit/
│   ├── setup_archaea.py            stage per-genus folders from <Genus>_gffs.txt
│   ├── setup_controls.py           stage per-genus folders from a sample/species tsv
│   ├── stage_phasomeit_inputs.job  SLURM wrapper for the two setup scripts
│   └── phasomeit_run.job           phasomeit array runner, cutoffs -c 7 6 0 5 5
│
├── parsing/
│   ├── parse_groups.py             phasomeit group html -> 4 csvs
│   ├── parse_tracts.py             phasomeit strain html -> tract rows
│   ├── summarise_runs.py           per-genus run completeness
│   ├── build_full_group_summary.py NEW — domain-labelled group rows for the notebook
│   ├── summarise_prokka_functions.py NEW — per-genus ranking of Prokka-annotated PPV functions
│   └── parse_phasomeit.job         SLURM wrapper for the parsing stage (chains all of the above)
│
├── eggnog/
│   ├── download_eggnog_db.job      NEW — one-time download of the eggnog v5 db
│   ├── eggnog_input_merger.job     concat prokka .faa files with sample-id prefixes
│   ├── run_eggnog_archaea.job      emapper, tax_scope Archaea
│   ├── run_eggnog_bacteria.job     emapper, default tax scope
│   ├── parse_eggnog.py             emapper.annotations -> tidy tsvs with broad roles
│   ├── link_eggnog.py              join phasomeit pv members to eggnog by locus tag
│   └── post_eggnog.job             SLURM wrapper for parse + link
│
├── quantification/
│   ├── pv_per_fasta.py             pv count per genome
│   ├── genome_lengths.py           genome size in bp and mb
│   ├── pv_per_mb.py                per-genome and per-genus pv/mb
│   └── quantify.job                SLURM wrapper that chains the three
│
├── analysis/                       NEW — method follow-up (Section 4.4 of report)
│   ├── tract_length_variation.py   within-species tract-length variance per gene group
│   └── tract_length_variation.job  SLURM wrapper
│
└── figures/
    ├── build_master_tables.ipynb
    ├── plot_pv_burden.ipynb
    └── plot_cog_functions.ipynb
```

## Pipeline summary

1. `dataset/prepare_dataset.job` — filtered species calls (per-species
   threshold), per-genus GFF lists.
2. `annotation/prokka_*.job` — Prokka annotation, run on every fasta.
3. `phasomeit/stage_phasomeit_inputs.job` — stage per-genus .gbk
   directories.
4. `phasomeit/phasomeit_run.job` — PhasomeIt array, one task per genus.
5. `parsing/parse_phasomeit.job` — HTML → CSV (group summary, tract data,
   members, pairwise, prokka function ranking).
6. `eggnog/run_eggnog_*.job` then `eggnog/post_eggnog.job` — annotate and
   link PV members to eggNOG / arCOG categories.
7. `quantification/quantify.job` — PV counts, genome sizes, PV / Mb.
8. `analysis/tract_length_variation.job` — within-species tract-length
   variance (method follow-up from Section 4.4 of the report).
9. `figures/*.ipynb` — Figures 1–4 and the master tables.

## Reproducing the original write-up's numbers

The numbers in the report were produced under the per-genus `>= 5`
threshold. To reproduce them exactly, run `dataset/build_species_calls.py`
with `--keep-singletons` (this disables the species threshold), and add
back the per-genus filter manually from `species_count_summary.csv`. The
in-repo version of `standards/SpeciesCallsArchaea.csv` is the per-species
output and will not exactly match the published numbers.

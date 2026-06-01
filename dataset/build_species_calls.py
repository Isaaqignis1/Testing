#!/usr/bin/env python3
"""build the filtered species-call table for the archaeal dataset

Reads the AllTheBacteria metadata table, drops unknown / non-archaeal
contaminant species, then applies a *per-species* minimum strain threshold.
Genera are retained if any of their species survive the threshold; only the
surviving species' samples are written out.

This replaces the per-genus threshold used in the original write-up
(Section 2.1, Methods) with a stricter per-species rule.

Output: SpeciesCallsArchaea.csv  (sample_id, scientific_name, genus, species)
        species_count_summary.csv  (per-species + per-genus counts before/after)
"""

from __future__ import annotations

import argparse
import csv
from collections import Counter, defaultdict
from pathlib import Path

import pandas as pd


# Confirmed non-archaeal contaminants identified during dataset assembly.
# Match the `scientific_name` column exactly (Genus species).
DEFAULT_CONTAMINANTS: set[str] = {
    "Streptococcus pneumoniae",
    "Streptococcus pseudopneumoniae",
    "Paraclostridium dentum",
    "Lysobacter enzymogenes",
    "Chromohalobacter canadensis",
}

UNKNOWN_LABELS: set[str] = {"unknown", "Unknown", "UNKNOWN", "", "NA", "N/A"}


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--meta", required=True, type=Path,
                   help="AllTheBacteria metadata TSV / CSV "
                        "(needs sample_id and scientific_name columns)")
    p.add_argument("--out-dir", required=True, type=Path)
    p.add_argument("--min-per-species", type=int, default=5,
                   help="minimum strain count required per species "
                        "(default: 5)")
    p.add_argument("--extra-contaminants", type=Path, default=None,
                   help="optional file with one extra contaminant "
                        "scientific name per line")
    p.add_argument("--keep-singletons", action="store_true",
                   help="report-only mode: do NOT drop species below "
                        "the threshold (for diagnostics)")
    return p.parse_args()


def detect_sep(path: Path) -> str:
    return "\t" if path.suffix.lower() in {".tsv", ".tab"} else ","


def load_metadata(path: Path) -> pd.DataFrame:
    df = pd.read_csv(path, sep=detect_sep(path), dtype=str).fillna("")
    rename = {}
    if "sample_accession" in df.columns and "sample_id" not in df.columns:
        rename["sample_accession"] = "sample_id"
    if "Species" in df.columns and "scientific_name" not in df.columns:
        rename["Species"] = "scientific_name"
    df = df.rename(columns=rename)
    if "sample_id" not in df.columns or "scientific_name" not in df.columns:
        raise SystemExit(
            f"{path}: need sample_id (or sample_accession) and "
            f"scientific_name (or Species) columns; found {list(df.columns)}"
        )
    return df[["sample_id", "scientific_name"]].copy()


def split_genus_species(scientific_name: str) -> tuple[str, str]:
    parts = scientific_name.strip().split(maxsplit=1)
    if len(parts) == 0:
        return "", ""
    if len(parts) == 1:
        return parts[0], ""
    return parts[0], parts[1]


def load_extra_contaminants(path: Path | None) -> set[str]:
    if path is None or not path.exists():
        return set()
    return {ln.strip() for ln in path.read_text().splitlines() if ln.strip()}


def main() -> None:
    args = parse_args()
    args.out_dir.mkdir(parents=True, exist_ok=True)

    contaminants = DEFAULT_CONTAMINANTS | load_extra_contaminants(
        args.extra_contaminants
    )

    meta = load_metadata(args.meta)
    n_total = len(meta)

    # --- Stage 1: drop unknowns + contaminants ------------------------------
    is_unknown = meta["scientific_name"].isin(UNKNOWN_LABELS) | meta[
        "scientific_name"
    ].str.lower().str.startswith("unknown")
    is_contam = meta["scientific_name"].isin(contaminants)

    dropped_unknown = is_unknown.sum()
    dropped_contam = is_contam.sum()

    kept = meta.loc[~is_unknown & ~is_contam].copy()
    kept[["genus", "species"]] = kept["scientific_name"].apply(
        lambda s: pd.Series(split_genus_species(s))
    )

    # --- Stage 2: per-species threshold -------------------------------------
    species_counts = (
        kept.groupby(["genus", "species"])
        .size()
        .reset_index(name="n_strains_pre")
    )

    threshold = args.min_per_species
    species_counts["passes_species_threshold"] = species_counts[
        "n_strains_pre"
    ].ge(threshold)

    surviving_species = set(
        species_counts.loc[
            species_counts["passes_species_threshold"],
            ["genus", "species"],
        ].apply(tuple, axis=1)
    )

    if args.keep_singletons:
        final = kept.copy()
        print("WARNING: --keep-singletons set, threshold NOT applied")
    else:
        keep_mask = kept[["genus", "species"]].apply(tuple, axis=1).isin(
            surviving_species
        )
        final = kept.loc[keep_mask].copy()

    # --- Write final species calls (one row per surviving sample) -----------
    final = final[["sample_id", "scientific_name", "genus", "species"]].sort_values(
        ["genus", "species", "sample_id"]
    )
    final_path = args.out_dir / "SpeciesCallsArchaea.csv"
    final.to_csv(final_path, index=False)

    # --- Write count summary ------------------------------------------------
    final_counts = (
        final.groupby(["genus", "species"])
        .size()
        .reset_index(name="n_strains_final")
    )
    summary = species_counts.merge(
        final_counts, on=["genus", "species"], how="left"
    )
    summary["n_strains_final"] = summary["n_strains_final"].fillna(0).astype(int)
    summary = summary.sort_values(["genus", "species"])

    summary_path = args.out_dir / "species_count_summary.csv"
    summary.to_csv(summary_path, index=False)

    # --- Per-genus summary --------------------------------------------------
    genus_summary = (
        final.groupby("genus")
        .agg(
            n_species_kept=("species", "nunique"),
            n_strains_kept=("sample_id", "nunique"),
        )
        .reset_index()
        .sort_values("genus")
    )
    genus_summary_path = args.out_dir / "genus_count_summary.csv"
    genus_summary.to_csv(genus_summary_path, index=False)

    # --- Stdout report ------------------------------------------------------
    print(f"input rows                  : {n_total:,}")
    print(f"dropped (unknown species)   : {dropped_unknown:,}")
    print(f"dropped (contaminants)      : {dropped_contam:,}")
    print(f"after stage 1               : {len(kept):,}")
    print(f"per-species threshold       : >= {threshold} strains")
    print(f"species kept                : {summary['passes_species_threshold'].sum():,} "
          f"/ {len(summary):,}")
    print(f"after stage 2               : {len(final):,} samples")
    print(f"genera retained             : {final['genus'].nunique()}")
    print()
    print(f"wrote {final_path}")
    print(f"wrote {summary_path}")
    print(f"wrote {genus_summary_path}")


if __name__ == "__main__":
    main()

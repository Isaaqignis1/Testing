#!/usr/bin/env python3
"""Build the filtered species-call table for the archaeal dataset.

Reads the ATB-filtered metadata, drops unknown / contaminant species, then
applies a per-species minimum strain threshold. A genus is retained if any
of its species clears the threshold; only surviving species' samples are
written out.

Outputs in --out-dir:
  SpeciesCallsArchaea.csv      sample_id, scientific_name, genus, species
  species_count_summary.csv    per-species pre/post counts + pass flag
  genus_count_summary.csv      per-genus rollup
"""

from __future__ import annotations

import argparse
from pathlib import Path

import pandas as pd


UNKNOWN_LABELS: set[str] = {"unknown", "Unknown", "UNKNOWN", "", "NA", "N/A"}


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--meta",            required=True, type=Path,
                   help="TSV/CSV with sample_id + scientific_name columns")
    p.add_argument("--out-dir",         required=True, type=Path)
    p.add_argument("--min-per-species", type=int, default=5)
    p.add_argument("--contaminants-file", type=Path, default=None,
                   help="one scientific_name per line; default contaminants "
                        "ship in standards/contaminants.txt")
    p.add_argument("--keep-singletons", action="store_true",
                   help="report-only mode: do NOT drop species below threshold")
    return p.parse_args()


def detect_sep(p: Path) -> str:
    return "\t" if p.suffix.lower() in {".tsv", ".tab"} else ","


def load_meta(p: Path) -> pd.DataFrame:
    df = pd.read_csv(p, sep=detect_sep(p), dtype=str).fillna("")
    rename = {}
    if "Sample"           in df.columns: rename["Sample"] = "sample_id"
    if "sample_accession" in df.columns: rename.setdefault("sample_accession", "sample_id")
    if "accession"        in df.columns: rename.setdefault("accession", "sample_id")
    if "Species"          in df.columns: rename["Species"] = "scientific_name"
    df = df.rename(columns=rename)
    miss = {"sample_id", "scientific_name"} - set(df.columns)
    if miss:
        raise SystemExit(f"{p}: missing columns {miss}; found {list(df.columns)}")
    return df[["sample_id", "scientific_name"]].copy()


def load_contaminants(path: Path | None) -> set[str]:
    if path is None or not path.exists():
        return set()
    return {ln.strip() for ln in path.read_text().splitlines()
            if ln.strip() and not ln.strip().startswith("#")}


def split_genus_species(s: str) -> tuple[str, str]:
    parts = s.strip().split(maxsplit=1)
    if not parts:           return "", ""
    if len(parts) == 1:     return parts[0], ""
    return parts[0], parts[1]


def main() -> None:
    a = parse_args()
    a.out_dir.mkdir(parents=True, exist_ok=True)

    contam = load_contaminants(a.contaminants_file)
    meta = load_meta(a.meta)
    n_total = len(meta)

    is_unknown = meta["scientific_name"].isin(UNKNOWN_LABELS) | \
                 meta["scientific_name"].str.lower().str.startswith("unknown")
    is_contam  = meta["scientific_name"].isin(contam)

    n_unknown = int(is_unknown.sum())
    n_contam  = int(is_contam.sum())
    kept = meta.loc[~is_unknown & ~is_contam].copy()
    kept[["genus", "species"]] = kept["scientific_name"].apply(
        lambda s: pd.Series(split_genus_species(s))
    )

    sp_counts = (kept.groupby(["genus", "species"])
                     .size().reset_index(name="n_strains_pre"))
    sp_counts["passes_species_threshold"] = sp_counts["n_strains_pre"].ge(a.min_per_species)
    survive = set(sp_counts.loc[sp_counts["passes_species_threshold"],
                                 ["genus", "species"]].apply(tuple, axis=1))

    if a.keep_singletons:
        final = kept.copy()
        print("WARNING: --keep-singletons set, threshold NOT applied")
    else:
        mask = kept[["genus", "species"]].apply(tuple, axis=1).isin(survive)
        final = kept.loc[mask].copy()

    final = final[["sample_id", "scientific_name", "genus", "species"]] \
        .sort_values(["genus", "species", "sample_id"])
    final_path = a.out_dir / "SpeciesCallsArchaea.csv"
    final.to_csv(final_path, index=False)

    fc = (final.groupby(["genus", "species"]).size()
                .reset_index(name="n_strains_final"))
    summary = sp_counts.merge(fc, on=["genus", "species"], how="left")
    summary["n_strains_final"] = summary["n_strains_final"].fillna(0).astype(int)
    summary.to_csv(a.out_dir / "species_count_summary.csv", index=False)

    gs = (final.groupby("genus")
              .agg(n_species_kept=("species", "nunique"),
                   n_strains_kept=("sample_id", "nunique"))
              .reset_index().sort_values("genus"))
    gs.to_csv(a.out_dir / "genus_count_summary.csv", index=False)

    print(f"input rows                : {n_total:,}")
    print(f"dropped (unknown)         : {n_unknown:,}")
    print(f"dropped (contaminants)    : {n_contam:,}")
    print(f"after stage 1             : {len(kept):,}")
    print(f"per-species threshold     : >= {a.min_per_species}")
    print(f"species kept              : {summary['passes_species_threshold'].sum():,} / {len(summary):,}")
    print(f"final samples             : {len(final):,}")
    print(f"genera retained           : {final['genus'].nunique()}")
    print(f"\nwrote {final_path}")
    print(f"wrote {a.out_dir / 'species_count_summary.csv'}")
    print(f"wrote {a.out_dir / 'genus_count_summary.csv'}")


if __name__ == "__main__":
    main()

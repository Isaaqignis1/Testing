#!/usr/bin/env python3
"""within-species tract-length variation across strains, per gene group

For each (genus, species, group_num) computes the spread of SSR tract
on-lengths across that species' strains. Groups with variation in tract
length across multiple strains are high-confidence PV candidates; groups
with identical tract lengths across all strains carry an SSR but no
detectable signature of switching (Section 4.4 of the report).

Inputs:
  --tracts            phasomeit_tract_data.csv (from parse_tracts.py)
  --species-calls     SpeciesCallsArchaea.csv  (from build_species_calls.py)
  --controls-calls    SpeciesCallsBacteria.tsv (controls species map)

Outputs:
  tract_length_variation_per_group.csv
        one row per (domain, genus, species, group_num); n_strains,
        n_unique_on_lengths, on_length_range, on_length_std,
        is_variable (True if >=2 distinct lengths and >=2 strains carry the SSR)
  tract_length_variation_per_genus.csv
        per-genus rollup; fraction of groups with variable tracts

Strains are linked to species via the species-calls tables on
strain_accession (== sample_id used everywhere else). Rows without a
species call are still counted but reported under species=='unknown'.
"""

from __future__ import annotations

import argparse
from pathlib import Path

import numpy as np
import pandas as pd


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--tracts", required=True, type=Path)
    p.add_argument("--species-calls",  required=True, type=Path,
                   help="SpeciesCallsArchaea.csv")
    p.add_argument("--controls-calls", required=True, type=Path,
                   help="SpeciesCallsBacteria.tsv")
    p.add_argument("--outdir", required=True, type=Path)
    p.add_argument("--pv-classes", nargs="+",
                   default=["intragenic", "upstream"],
                   help="pv_class values to include (default: intragenic, upstream)")
    return p.parse_args()


def load_species_calls(path: Path) -> pd.DataFrame:
    sep = "\t" if path.suffix.lower() in {".tsv", ".tab"} else ","
    df = pd.read_csv(path, sep=sep, dtype=str)
    rename = {}
    if "sample_accession" in df.columns:
        rename["sample_accession"] = "sample_id"
    df = df.rename(columns=rename)
    if "sample_id" not in df.columns:
        raise SystemExit(f"{path}: needs sample_id (or sample_accession)")
    if "species" not in df.columns and "scientific_name" in df.columns:
        df["species"] = df["scientific_name"].str.split(
            n=1
        ).str[1].fillna("")
    if "genus" not in df.columns and "scientific_name" in df.columns:
        df["genus"] = df["scientific_name"].str.split(n=1).str[0]
    return df[["sample_id", "genus", "species"]].fillna("")


def main() -> None:
    args = parse_args()
    args.outdir.mkdir(parents=True, exist_ok=True)

    tracts = pd.read_csv(args.tracts, low_memory=False)
    needed = {"domain", "genus", "strain_accession", "group_num",
              "pv_class", "on_length"}
    miss = needed - set(tracts.columns)
    if miss:
        raise SystemExit(f"{args.tracts}: missing columns {miss}")

    tracts = tracts[tracts["pv_class"].isin(args.pv_classes)].copy()
    tracts["on_length"] = pd.to_numeric(tracts["on_length"], errors="coerce")
    tracts = tracts.dropna(subset=["on_length", "group_num"])
    tracts["group_num"] = tracts["group_num"].astype(int)

    species = pd.concat(
        [
            load_species_calls(args.species_calls).assign(domain="Archaea"),
            load_species_calls(args.controls_calls).assign(domain="Bacteria"),
        ],
        ignore_index=True,
    )

    merged = tracts.merge(
        species, left_on=["domain", "strain_accession"],
        right_on=["domain", "sample_id"], how="left",
    )
    merged["species"] = merged["species"].fillna("unknown").replace("", "unknown")
    # use the genus from the tract data (consistent with the rest of the pipeline)
    merged = merged.drop(columns=["sample_id", "genus_y"], errors="ignore")
    if "genus_x" in merged.columns:
        merged = merged.rename(columns={"genus_x": "genus"})

    # one row per strain x group x tract; collapse to unique on-lengths per strain
    # by taking the modal on_length (most common) — strains sometimes have
    # multiple tracts in the same group, we want the dominant tract length
    strain_group = (
        merged.groupby(["domain", "genus", "species",
                        "strain_accession", "group_num"])["on_length"]
              .agg(lambda s: s.mode().iat[0] if not s.mode().empty else np.nan)
              .reset_index()
    )

    # per-group stats within (domain, genus, species, group_num)
    per_group = (
        strain_group
        .groupby(["domain", "genus", "species", "group_num"])
        .agg(
            n_strains_with_tract=("on_length", "count"),
            n_unique_on_lengths=("on_length", "nunique"),
            on_length_min=("on_length", "min"),
            on_length_max=("on_length", "max"),
            on_length_mean=("on_length", "mean"),
            on_length_std=("on_length", "std"),
        )
        .reset_index()
    )
    per_group["on_length_range"] = (
        per_group["on_length_max"] - per_group["on_length_min"]
    )
    per_group["is_variable"] = (
        (per_group["n_unique_on_lengths"] >= 2)
        & (per_group["n_strains_with_tract"] >= 2)
    )
    per_group = per_group.sort_values(
        ["domain", "genus", "species", "group_num"]
    )
    per_group.to_csv(args.outdir / "tract_length_variation_per_group.csv",
                     index=False)

    # rollup per genus
    per_genus = (
        per_group
        .groupby(["domain", "genus"])
        .agg(
            n_groups=("group_num", "count"),
            n_variable=("is_variable", "sum"),
            mean_unique_on_lengths=("n_unique_on_lengths", "mean"),
            mean_on_length_range=("on_length_range", "mean"),
        )
        .reset_index()
    )
    per_genus["frac_variable"] = (
        per_genus["n_variable"] / per_genus["n_groups"]
    ).round(4)
    per_genus = per_genus.sort_values(
        ["domain", "frac_variable"], ascending=[True, False]
    )
    per_genus.to_csv(args.outdir / "tract_length_variation_per_genus.csv",
                     index=False)

    print(f"wrote tract_length_variation_per_group.csv "
          f"({len(per_group):,} rows)")
    print(f"wrote tract_length_variation_per_genus.csv "
          f"({len(per_genus):,} rows)")
    print()
    print("top genera by fraction of variable groups:")
    print(per_genus.head(10).to_string(index=False))


if __name__ == "__main__":
    main()

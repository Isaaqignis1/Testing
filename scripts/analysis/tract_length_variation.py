#!/usr/bin/env python3
"""within-species tract-length variation across strains, per gene group.

For each (genus, species, group_num) computes the spread of SSR on-lengths
across that species' strains. Variable groups are high-confidence PV
candidates; fixed groups carry an SSR but no detectable switching signature
(Section 4.4 of the report).

Outputs:
  tract_length_variation_per_group.csv
  tract_length_variation_per_genus.csv
"""

from __future__ import annotations

import argparse
from pathlib import Path

import numpy as np
import pandas as pd


def load_species_calls(path: Path) -> pd.DataFrame:
    sep = "\t" if path.suffix.lower() in {".tsv", ".tab"} else ","
    df = pd.read_csv(path, sep=sep, dtype=str)
    if "sample_accession" in df.columns:
        df = df.rename(columns={"sample_accession": "sample_id"})
    if "sample_id" not in df.columns:
        raise SystemExit(f"{path}: needs sample_id")
    if "species" not in df.columns and "scientific_name" in df.columns:
        df["species"] = df["scientific_name"].str.split(n=1).str[1].fillna("")
    if "genus" not in df.columns and "scientific_name" in df.columns:
        df["genus"] = df["scientific_name"].str.split(n=1).str[0]
    return df[["sample_id", "genus", "species"]].fillna("")


def main() -> None:
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--tracts", required=True, type=Path)
    p.add_argument("--species-calls",  required=True, type=Path)
    p.add_argument("--controls-calls", required=True, type=Path)
    p.add_argument("--outdir", required=True, type=Path)
    p.add_argument("--pv-classes", nargs="+", default=["intragenic", "upstream"])
    a = p.parse_args()
    a.outdir.mkdir(parents=True, exist_ok=True)

    tr = pd.read_csv(a.tracts, low_memory=False)
    need = {"domain", "genus", "strain_accession", "group_num", "pv_class", "on_length"}
    miss = need - set(tr.columns)
    if miss: raise SystemExit(f"{a.tracts}: missing {miss}")

    tr = tr[tr["pv_class"].isin(a.pv_classes)].copy()
    tr["on_length"] = pd.to_numeric(tr["on_length"], errors="coerce")
    tr = tr.dropna(subset=["on_length", "group_num"])
    tr["group_num"] = tr["group_num"].astype(int)

    sp = pd.concat([
        load_species_calls(a.species_calls).assign(domain="Archaea"),
        load_species_calls(a.controls_calls).assign(domain="Bacteria"),
    ], ignore_index=True)

    m = tr.merge(sp, left_on=["domain", "strain_accession"],
                  right_on=["domain", "sample_id"], how="left")
    m["species"] = m["species"].fillna("unknown").replace("", "unknown")
    m = m.drop(columns=["sample_id", "genus_y"], errors="ignore")
    if "genus_x" in m.columns:
        m = m.rename(columns={"genus_x": "genus"})

    sg = (m.groupby(["domain", "genus", "species",
                      "strain_accession", "group_num"])["on_length"]
            .agg(lambda s: s.mode().iat[0] if not s.mode().empty else np.nan)
            .reset_index())

    per_g = (sg.groupby(["domain", "genus", "species", "group_num"])
                .agg(n_strains_with_tract=("on_length", "count"),
                     n_unique_on_lengths=("on_length", "nunique"),
                     on_length_min=("on_length", "min"),
                     on_length_max=("on_length", "max"),
                     on_length_mean=("on_length", "mean"),
                     on_length_std=("on_length", "std"))
                .reset_index())
    per_g["on_length_range"] = per_g["on_length_max"] - per_g["on_length_min"]
    per_g["is_variable"] = (per_g["n_unique_on_lengths"] >= 2) & (per_g["n_strains_with_tract"] >= 2)
    per_g = per_g.sort_values(["domain", "genus", "species", "group_num"])
    per_g.to_csv(a.outdir / "tract_length_variation_per_group.csv", index=False)

    per_genus = (per_g.groupby(["domain", "genus"])
                       .agg(n_groups=("group_num", "count"),
                            n_variable=("is_variable", "sum"),
                            mean_unique_on_lengths=("n_unique_on_lengths", "mean"),
                            mean_on_length_range=("on_length_range", "mean"))
                       .reset_index())
    per_genus["frac_variable"] = (per_genus["n_variable"] / per_genus["n_groups"]).round(4)
    per_genus = per_genus.sort_values(["domain", "frac_variable"], ascending=[True, False])
    per_genus.to_csv(a.outdir / "tract_length_variation_per_genus.csv", index=False)

    print(f"wrote per_group  ({len(per_g):,} rows)")
    print(f"wrote per_genus  ({len(per_genus):,} rows)")
    print("\ntop genera by fraction of variable groups:")
    print(per_genus.head(10).to_string(index=False))


if __name__ == "__main__":
    main()

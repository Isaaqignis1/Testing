#!/usr/bin/env python3
"""join pv counts with genome lengths and emit per-genome + per-genus pv/mb."""

from __future__ import annotations

import argparse
from pathlib import Path

import pandas as pd


def load_lengths(p: Path, source: str) -> pd.DataFrame:
    df = pd.read_csv(p)[["fasta_name", "genus", "genome_size_bp",
                          "genome_size_mb"]].copy()
    df["source"] = source
    return df


def main() -> None:
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--pv-counts",       required=True, type=Path)
    p.add_argument("--lengths-main",    required=True, type=Path)
    p.add_argument("--lengths-control", required=True, type=Path)
    p.add_argument("--out-genome",      required=True, type=Path)
    p.add_argument("--out-genus",       required=True, type=Path)
    a = p.parse_args()
    a.out_genome.parent.mkdir(parents=True, exist_ok=True)
    a.out_genus.parent.mkdir(parents=True, exist_ok=True)

    pv = pd.read_csv(a.pv_counts)
    lengths = pd.concat([
        load_lengths(a.lengths_main, "main"),
        load_lengths(a.lengths_control, "control"),
    ], ignore_index=True)

    merged = pv.merge(
        lengths[["source", "fasta_name", "genome_size_bp", "genome_size_mb"]],
        on=["source", "fasta_name"], how="left")
    if "genus_x" in merged.columns:
        merged = merged.rename(columns={"genus_x": "genus"}).drop(columns=["genus_y"])
    merged["ppv_per_mb"] = merged["PV_count"] / merged["genome_size_mb"]
    merged = merged[["source", "genus", "fasta_name",
                     "PV_count", "genome_size_bp", "genome_size_mb",
                     "ppv_per_mb"]].sort_values(["source", "genus", "fasta_name"])
    merged.to_csv(a.out_genome, index=False)
    print(f"per-genome: {a.out_genome} ({len(merged):,})")

    gs = (merged.groupby(["source", "genus"])
                 .agg(n_genomes=("fasta_name", "count"),
                      mean_PV=("PV_count", "mean"),
                      sd_PV=("PV_count", "std"),
                      mean_genome_size_mb=("genome_size_mb", "mean"),
                      mean_ppv_per_mb=("ppv_per_mb", "mean"),
                      sd_ppv_per_mb=("ppv_per_mb", "std"))
                 .round(6).reset_index()
                 .sort_values(["source", "mean_ppv_per_mb"],
                              ascending=[True, False]))
    gs.to_csv(a.out_genus, index=False)
    print(f"per-genus:  {a.out_genus} ({len(gs):,})")


if __name__ == "__main__":
    main()

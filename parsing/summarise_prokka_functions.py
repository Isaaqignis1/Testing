#!/usr/bin/env python3
"""per-genus ranking of top Prokka-annotated PPV gene-group functions

Operates on the group-level summary produced by parse_groups.py.
Outputs three tables under <outdir>/:

  pv_function_counts_by_genus.csv
      one row per (genus, function); groups, sum and average pv_in_gene
  pv_function_counts_total.csv
      one row per function, summed across genera
  pv_genus_summary.csv
      one row per genus, totals + per-genome rates (n_genomes from
      summarise_runs.py output, if provided)

This view retains the bacterial-namesake Prokka labels rather than COG
categories, so it is reported in the repo as a supplementary table
(Methods 2.6 and Discussion 4.2.1 of the report).
"""

from __future__ import annotations

import argparse
from collections import defaultdict
from pathlib import Path

import pandas as pd


def clean_function(func: str) -> str:
    if not isinstance(func, str) or not func.strip():
        return "unknown"
    func = func.strip()
    if func.lower() in {
        "hypothetical protein",
        "putative hypothetical protein",
        "conserved hypothetical protein",
    }:
        return "hypothetical protein"
    return func


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--groups", required=True, type=Path,
                   help="phasomeit_group_summary.csv (parse_groups.py output)")
    p.add_argument("--genus-summary", type=Path, default=None,
                   help="optional summarise_runs.py output for per-genome rates")
    p.add_argument("--outdir", required=True, type=Path)
    p.add_argument("--min-pv-in-gene", type=int, default=1,
                   help="only count groups with pv_in_gene >= this (default 1)")
    return p.parse_args()


def main() -> None:
    args = parse_args()
    args.outdir.mkdir(parents=True, exist_ok=True)

    df = pd.read_csv(args.groups)
    df["function"] = df["likely_function"].map(clean_function)
    for col in ("pv_in_gene", "total_pv", "total_genes"):
        df[col] = pd.to_numeric(df[col], errors="coerce").fillna(0).astype(int)

    df = df.loc[df["pv_in_gene"] >= args.min_pv_in_gene].copy()

    # --- 1. counts by genus + function --------------------------------------
    by_gf = (
        df.groupby(["genus", "function"])
          .agg(
              n_groups=("group_num", "count"),
              sum_pv_in_gene=("pv_in_gene", "sum"),
              sum_total_pv=("total_pv",   "sum"),
              sum_total_genes=("total_genes", "sum"),
          )
          .reset_index()
    )
    by_gf["avg_pv_in_gene"] = (by_gf["sum_pv_in_gene"] / by_gf["n_groups"]).round(4)
    by_gf = by_gf.sort_values(["genus", "sum_pv_in_gene"], ascending=[True, False])
    by_gf.to_csv(args.outdir / "pv_function_counts_by_genus.csv", index=False)

    # --- 2. counts total ----------------------------------------------------
    total = (
        df.groupby("function")
          .agg(
              n_groups=("group_num", "count"),
              n_genera=("genus", "nunique"),
              sum_pv_in_gene=("pv_in_gene", "sum"),
              sum_total_pv=("total_pv",   "sum"),
              sum_total_genes=("total_genes", "sum"),
          )
          .reset_index()
    )
    total["avg_pv_in_gene"] = (total["sum_pv_in_gene"] / total["n_groups"]).round(4)
    total = total.sort_values(["n_groups", "function"], ascending=[False, True])
    total.to_csv(args.outdir / "pv_function_counts_total.csv", index=False)

    # --- 3. per-genus summary ----------------------------------------------
    per_genus = (
        df.groupby("genus")
          .agg(
              pv_groups=("group_num", "count"),
              sum_pv_in_gene=("pv_in_gene", "sum"),
              sum_total_pv=("total_pv",   "sum"),
              sum_total_genes=("total_genes", "sum"),
          )
          .reset_index()
    )

    if args.genus_summary is not None and args.genus_summary.exists():
        gs = pd.read_csv(args.genus_summary)
        col = "n_genomes_gbk" if "n_genomes_gbk" in gs.columns else "n_genomes"
        per_genus = per_genus.merge(
            gs[["genus", col]].rename(columns={col: "n_genomes"}),
            on="genus", how="left",
        )
        per_genus["pv_groups_per_genome"] = (
            per_genus["pv_groups"] / per_genus["n_genomes"]
        ).round(4)
        per_genus["pv_in_gene_per_genome"] = (
            per_genus["sum_pv_in_gene"] / per_genus["n_genomes"]
        ).round(4)
    else:
        per_genus["n_genomes"] = pd.NA
        per_genus["pv_groups_per_genome"] = pd.NA
        per_genus["pv_in_gene_per_genome"] = pd.NA

    per_genus.to_csv(args.outdir / "pv_genus_summary.csv", index=False)

    print(f"wrote pv_function_counts_by_genus.csv ({len(by_gf):,} rows)")
    print(f"wrote pv_function_counts_total.csv    ({len(total):,} rows)")
    print(f"wrote pv_genus_summary.csv            ({len(per_genus):,} rows)")
    print(f"outdir: {args.outdir}")


if __name__ == "__main__":
    main()

#!/usr/bin/env python3
"""per-genus ranking of top Prokka-annotated PPV gene-group functions.

Operates on the group-level summary produced by parse_groups.py. Outputs:
  pv_function_counts_by_genus.csv
  pv_function_counts_total.csv
  pv_genus_summary.csv         (per-genome rates if --genus-summary given)
"""

from __future__ import annotations

import argparse
from pathlib import Path

import pandas as pd


def clean(f):
    if not isinstance(f, str) or not f.strip(): return "unknown"
    f = f.strip()
    if f.lower() in {"hypothetical protein", "putative hypothetical protein",
                     "conserved hypothetical protein"}:
        return "hypothetical protein"
    return f


def main() -> None:
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--groups", required=True, type=Path)
    p.add_argument("--genus-summary", type=Path, default=None)
    p.add_argument("--outdir", required=True, type=Path)
    p.add_argument("--min-pv-in-gene", type=int, default=1)
    a = p.parse_args()
    a.outdir.mkdir(parents=True, exist_ok=True)

    df = pd.read_csv(a.groups)
    df["function"] = df["likely_function"].map(clean)
    for c in ("pv_in_gene", "total_pv", "total_genes"):
        df[c] = pd.to_numeric(df[c], errors="coerce").fillna(0).astype(int)
    df = df.loc[df["pv_in_gene"] >= a.min_pv_in_gene].copy()

    by_gf = (df.groupby(["genus", "function"])
               .agg(n_groups=("group_num", "count"),
                    sum_pv_in_gene=("pv_in_gene", "sum"),
                    sum_total_pv=("total_pv", "sum"),
                    sum_total_genes=("total_genes", "sum"))
               .reset_index())
    by_gf["avg_pv_in_gene"] = (by_gf["sum_pv_in_gene"] / by_gf["n_groups"]).round(4)
    by_gf = by_gf.sort_values(["genus", "sum_pv_in_gene"], ascending=[True, False])
    by_gf.to_csv(a.outdir / "pv_function_counts_by_genus.csv", index=False)

    total = (df.groupby("function")
               .agg(n_groups=("group_num", "count"),
                    n_genera=("genus", "nunique"),
                    sum_pv_in_gene=("pv_in_gene", "sum"),
                    sum_total_pv=("total_pv", "sum"),
                    sum_total_genes=("total_genes", "sum"))
               .reset_index())
    total["avg_pv_in_gene"] = (total["sum_pv_in_gene"] / total["n_groups"]).round(4)
    total = total.sort_values(["n_groups", "function"], ascending=[False, True])
    total.to_csv(a.outdir / "pv_function_counts_total.csv", index=False)

    per_g = (df.groupby("genus")
               .agg(pv_groups=("group_num", "count"),
                    sum_pv_in_gene=("pv_in_gene", "sum"),
                    sum_total_pv=("total_pv", "sum"),
                    sum_total_genes=("total_genes", "sum"))
               .reset_index())
    if a.genus_summary and a.genus_summary.exists():
        gs = pd.read_csv(a.genus_summary)
        col = "n_genomes_gbk" if "n_genomes_gbk" in gs.columns else "n_genomes"
        per_g = per_g.merge(gs[["genus", col]].rename(columns={col: "n_genomes"}),
                            on="genus", how="left")
        per_g["pv_groups_per_genome"]   = (per_g["pv_groups"] / per_g["n_genomes"]).round(4)
        per_g["pv_in_gene_per_genome"]  = (per_g["sum_pv_in_gene"] / per_g["n_genomes"]).round(4)
    else:
        per_g["n_genomes"] = pd.NA
        per_g["pv_groups_per_genome"] = pd.NA
        per_g["pv_in_gene_per_genome"] = pd.NA
    per_g.to_csv(a.outdir / "pv_genus_summary.csv", index=False)

    print(f"wrote 3 files to {a.outdir}")


if __name__ == "__main__":
    main()

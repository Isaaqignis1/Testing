#!/usr/bin/env python3
"""Supplementary figure — top Prokka-annotated PPV gene-group functions
per outlier genus (Methanobacterium_D, Methanococcus, Methanobrevibacter_A,
Methanohalophilus, Nitrosopelagicus).

Input:
  --counts-by-genus  OUTPUTS/04_extraction/archaea/prokka_function_ranking/pv_function_counts_by_genus.csv

Output:
  supp_top_prokka_functions.pdf|png
"""

from __future__ import annotations

import argparse
from pathlib import Path

import matplotlib.pyplot as plt
import pandas as pd


DEFAULT_GENERA = [
    "Methanobacterium_D", "Methanococcus", "Methanobrevibacter_A",
    "Methanohalophilus", "Nitrosopelagicus",
]


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--counts-by-genus", required=True, type=Path)
    p.add_argument("--outdir",          required=True, type=Path)
    p.add_argument("--genera", nargs="+", default=DEFAULT_GENERA)
    p.add_argument("--top-n", type=int, default=10)
    return p.parse_args()


def main() -> None:
    a = parse_args()
    a.outdir.mkdir(parents=True, exist_ok=True)
    df = pd.read_csv(a.counts_by_genus)
    df = df[df["genus"].isin(a.genera)].copy()
    df = df[df["function"] != "hypothetical protein"]
    df = df.sort_values(["genus", "sum_pv_in_gene"], ascending=[True, False])

    fig, axes = plt.subplots(len(a.genera), 1,
                              figsize=(11, 2.2 * len(a.genera)),
                              constrained_layout=True)
    if len(a.genera) == 1:
        axes = [axes]
    for ax, genus in zip(axes, a.genera):
        sub = df[df["genus"] == genus].head(a.top_n)
        ax.barh(sub["function"], sub["sum_pv_in_gene"],
                 color="#3a7bb3", edgecolor="black", linewidth=0.4)
        ax.invert_yaxis()
        ax.set_title(f"{genus} — top {a.top_n} Prokka PPV functions",
                      loc="left", fontsize=10)
        ax.set_xlabel("Sum pv_in_gene")
        ax.grid(axis="x", alpha=0.3)
        ax.tick_params(axis="y", labelsize=8)

    for ext in ("pdf", "png"):
        fig.savefig(a.outdir / f"supp_top_prokka_functions.{ext}", dpi=200)
    plt.close(fig)
    print(f"wrote supplementary figure to {a.outdir}")


if __name__ == "__main__":
    main()

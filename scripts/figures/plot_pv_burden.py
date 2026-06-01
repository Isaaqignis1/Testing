#!/usr/bin/env python3
"""Figure 1 (scatter) + Figure 2 (boxplot) — PV burden across genera.

Inputs:
  --pv-per-mb     OUTPUTS/07_quantification/pv_per_mb_per_genome.csv
                  (one row per genome: source, genus, PV_count, ppv_per_mb)
  --pv-per-genus  OUTPUTS/07_quantification/pv_per_mb_per_genus.csv
                  (one row per genus, used for the scatter)

Outputs in --outdir:
  fig1_scatter_pv_per_strain_vs_pv_per_mb.pdf
  fig1_scatter_pv_per_strain_vs_pv_per_mb.png
  fig2_box_pv_per_mb.pdf
  fig2_box_pv_per_mb.png
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

import matplotlib.pyplot as plt
import matplotlib.lines as mlines
import numpy as np
import pandas as pd


def load_colors(standards_dir: Path):
    sys.path.insert(0, str(standards_dir))
    from colors import color_map, genus_order, control_color_map, control_order  # type: ignore
    return color_map, genus_order, control_color_map, control_order


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--pv-per-mb",    required=True, type=Path)
    p.add_argument("--pv-per-genus", required=True, type=Path)
    p.add_argument("--standards",    required=True, type=Path,
                   help="path to the standards/ directory")
    p.add_argument("--outdir",       required=True, type=Path)
    return p.parse_args()


def main() -> None:
    a = parse_args()
    a.outdir.mkdir(parents=True, exist_ok=True)
    color_map, genus_order, ctrl_map, ctrl_order = load_colors(a.standards)

    per_g = pd.read_csv(a.pv_per_genus)
    per_g = per_g.rename(columns={"mean_PV": "pv_per_strain",
                                    "mean_ppv_per_mb": "pv_per_mb"})
    per_g["domain"] = np.where(per_g["source"].eq("control"), "Bacteria", "Archaea")

    arch = per_g[per_g["domain"] == "Archaea"]
    bact = per_g[per_g["domain"] == "Bacteria"]

    # -------- Figure 1 ----------------------------------------------------
    fig, ax = plt.subplots(figsize=(7, 5))
    for _, r in arch.iterrows():
        ax.scatter(r["pv_per_strain"], r["pv_per_mb"],
                    c=color_map.get(r["genus"], "#888"),
                    s=160, marker="o", edgecolors="black",
                    linewidths=0.6, alpha=0.9, zorder=3)
    for _, r in bact.iterrows():
        ax.scatter(r["pv_per_strain"], r["pv_per_mb"],
                    c=ctrl_map.get(r["genus"], "#444"),
                    s=160, marker="D", edgecolors="black",
                    linewidths=0.6, alpha=0.9, zorder=4)
    for _, r in per_g.iterrows():
        ax.annotate(r["genus"], (r["pv_per_strain"], r["pv_per_mb"]),
                     xytext=(6, 4), textcoords="offset points",
                     fontsize=8, color="#222")
    ax.set_xlabel("Mean PPV per strain")
    ax.set_ylabel("Mean PPV per Mb")
    ax.set_title("Figure 1 — PPV burden by genus")
    ax.grid(alpha=0.3)
    handles = [
        mlines.Line2D([], [], color="black", marker="o", linestyle="None",
                       markersize=10, label="Archaea"),
        mlines.Line2D([], [], color="black", marker="D", linestyle="None",
                       markersize=10, label="Bacteria reference"),
    ]
    ax.legend(handles=handles, loc="upper left", frameon=False)
    fig.tight_layout()
    for ext in ("pdf", "png"):
        fig.savefig(a.outdir / f"fig1_scatter_pv_per_strain_vs_pv_per_mb.{ext}",
                     dpi=200)
    plt.close(fig)

    # -------- Figure 2 — boxplot of per-genome ppv_per_mb -----------------
    pg = pd.read_csv(a.pv_per_mb)
    pg["domain"] = np.where(pg["source"].eq("control"), "Bacteria", "Archaea")

    ordered = [g for g in genus_order if g in pg["genus"].unique()] + \
              [g for g in ctrl_order if g in pg["genus"].unique()]
    pg["__o"] = pg["genus"].apply(lambda g: ordered.index(g) if g in ordered else 9999)
    pg = pg.sort_values("__o")

    fig, ax = plt.subplots(figsize=(12, 5))
    data = [pg.loc[pg["genus"] == g, "ppv_per_mb"].dropna().to_numpy()
            for g in ordered]
    bp = ax.boxplot(data, patch_artist=True, showfliers=True, widths=0.7)
    for patch, g in zip(bp["boxes"], ordered):
        col = color_map.get(g) or ctrl_map.get(g) or "#cccccc"
        patch.set_facecolor(col); patch.set_alpha(0.85)
    ax.set_xticks(range(1, len(ordered) + 1))
    ax.set_xticklabels(ordered, rotation=70, ha="right", fontsize=8)
    ax.set_ylabel("PPV per Mb")
    ax.set_title("Figure 2 — Per-genome PPV/Mb by genus")
    # vertical line separating archaea from bacteria
    n_arch = len([g for g in ordered if g in genus_order])
    if 0 < n_arch < len(ordered):
        ax.axvline(n_arch + 0.5, color="grey", ls="--", lw=0.8)
    ax.grid(axis="y", alpha=0.3)
    fig.tight_layout()
    for ext in ("pdf", "png"):
        fig.savefig(a.outdir / f"fig2_box_pv_per_mb.{ext}", dpi=200)
    plt.close(fig)

    print(f"wrote figures to {a.outdir}")


if __name__ == "__main__":
    main()

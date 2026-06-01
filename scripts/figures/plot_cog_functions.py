#!/usr/bin/env python3
"""Figure 3 (stacked bar + log per-strain) + Figure 4 (heatmap) — COG
functional composition of PPV gene groups across genera.

Inputs:
  --master-groups     OUTPUTS/09_figures/master_groups.csv (build_master_tables.py)
  --strain-counts     OUTPUTS/02_dataset/species_calls/genus_count_summary.csv
                      (used to normalise to per-strain rates)

Outputs:
  fig3a_stacked_broad_role.pdf|png
  fig3b_per_strain_by_broad_role.pdf|png
  fig4a_heatmap_count_per_strain.pdf|png
  fig4b_heatmap_pct.pdf|png
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
import numpy as np
import pandas as pd


BROAD_ORDER = [
    "Surface & Defence",
    "Environmental Response",
    "Adaptive Metabolism",
    "Core cellular Machinery",
    "Poorly Characterised",
    "Other",
]

# Functional order within heatmaps (Figure 4)
COG_FUNCTION_ORDER = [
    "Cell wall/membrane/envelope biogenesis",
    "Cell motility",
    "Defense mechanisms",
    "Signal transduction",
    "Transcription",
    "Inorganic ion transport and metabolism",
    "Carbohydrate transport and metabolism",
    "Secondary metabolite biosynthesis and catabolism",
    "Amino acid transport and metabolism",
    "Energy production and conversion",
    "Coenzyme transport and metabolism",
    "Lipid transport and metabolism",
    "Nucleotide transport and metabolism",
    "Translation",
    "Replication, recombination and repair",
    "Post-translational modification and chaperones",
    "Cell cycle control and division",
    "Intracellular trafficking and secretion",
]


def load_colors(standards: Path):
    sys.path.insert(0, str(standards))
    from colors import color_map, genus_order, control_color_map, control_order  # type: ignore
    return color_map, genus_order, control_color_map, control_order


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--master-groups", required=True, type=Path)
    p.add_argument("--strain-counts", required=True, type=Path)
    p.add_argument("--standards",     required=True, type=Path)
    p.add_argument("--outdir",        required=True, type=Path)
    p.add_argument("--use-filtered",  action="store_true", default=True)
    return p.parse_args()


def main() -> None:
    a = parse_args()
    a.outdir.mkdir(parents=True, exist_ok=True)
    _, genus_order, _, ctrl_order = load_colors(a.standards)
    ordered = genus_order + ctrl_order

    mg = pd.read_csv(a.master_groups, low_memory=False)
    if a.use_filtered and "passes_filter" in mg.columns:
        mg = mg[mg["passes_filter"] == True].copy()

    # Strain counts per genus (for per-strain normalisation)
    sc = pd.read_csv(a.strain_counts)
    sc_col = "n_strains_kept" if "n_strains_kept" in sc.columns else "n_strains"
    strain_n = sc.set_index("genus")[sc_col].to_dict()

    # --- Figure 3a: stacked broad-role composition (proportional) ----------
    have_role = mg.dropna(subset=["broad_role"]).copy()
    counts = (have_role.groupby(["genus", "broad_role"]).size()
                          .reset_index(name="n"))
    pivot = counts.pivot(index="genus", columns="broad_role", values="n").fillna(0)
    pivot = pivot.reindex([g for g in ordered if g in pivot.index])
    pivot_pct = pivot.div(pivot.sum(axis=1), axis=0) * 100
    pivot_pct = pivot_pct.reindex(columns=[c for c in BROAD_ORDER if c in pivot_pct.columns])

    fig, ax = plt.subplots(figsize=(12, 6))
    bottom = np.zeros(len(pivot_pct))
    palette = plt.cm.tab20.colors
    for i, c in enumerate(pivot_pct.columns):
        ax.bar(pivot_pct.index, pivot_pct[c], bottom=bottom,
                color=palette[i % len(palette)], label=c, edgecolor="white", linewidth=0.5)
        bottom += pivot_pct[c].values
    ax.set_ylabel("% of annotated PPV gene groups")
    ax.set_title("Figure 3a — Functional composition by broad COG role")
    ax.set_xticklabels(pivot_pct.index, rotation=70, ha="right", fontsize=8)
    ax.legend(loc="center left", bbox_to_anchor=(1.02, 0.5), frameon=False, fontsize=8)
    fig.tight_layout()
    for ext in ("pdf", "png"):
        fig.savefig(a.outdir / f"fig3a_stacked_broad_role.{ext}", dpi=200)
    plt.close(fig)

    # --- Figure 3b: per-strain counts by broad role (log) ------------------
    per_strain = counts.copy()
    per_strain["n_strains"] = per_strain["genus"].map(strain_n)
    per_strain["count_per_strain"] = per_strain["n"] / per_strain["n_strains"]
    pps = per_strain.pivot(index="genus", columns="broad_role",
                            values="count_per_strain").fillna(0)
    pps = pps.reindex([g for g in ordered if g in pps.index])

    fig, ax = plt.subplots(figsize=(12, 6))
    width = 0.13
    x = np.arange(len(pps))
    for i, c in enumerate([c for c in BROAD_ORDER if c in pps.columns]):
        ax.bar(x + i * width - (len([c for c in BROAD_ORDER if c in pps.columns]) * width / 2),
                pps[c], width, label=c)
    ax.set_yscale("log")
    ax.set_xticks(x)
    ax.set_xticklabels(pps.index, rotation=70, ha="right", fontsize=8)
    ax.set_ylabel("PPV groups per strain (log)")
    ax.set_title("Figure 3b — Per-strain PPV groups by broad COG role")
    ax.legend(loc="center left", bbox_to_anchor=(1.02, 0.5), frameon=False, fontsize=8)
    fig.tight_layout()
    for ext in ("pdf", "png"):
        fig.savefig(a.outdir / f"fig3b_per_strain_by_broad_role.{ext}", dpi=200)
    plt.close(fig)

    # --- Figure 4: heatmaps ------------------------------------------------
    have_fn = mg.dropna(subset=["COG_function"]).copy()
    have_fn = have_fn[have_fn["COG_function"] != "Unknown"]
    fc = have_fn.groupby(["genus", "COG_function"]).size().reset_index(name="n")
    fc["n_strains"] = fc["genus"].map(strain_n)
    fc["count_per_strain"] = fc["n"] / fc["n_strains"]

    cps = (fc.pivot(index="genus", columns="COG_function",
                      values="count_per_strain").fillna(0))
    cps = cps.reindex([g for g in ordered if g in cps.index])
    cps = cps.reindex(columns=[c for c in COG_FUNCTION_ORDER if c in cps.columns])

    def heatmap(matrix: pd.DataFrame, title: str, fname: str, *, log_norm: bool):
        fig, ax = plt.subplots(figsize=(14, 7))
        if log_norm:
            norm = mcolors.LogNorm(vmin=max(matrix.values[matrix.values > 0].min(), 0.001),
                                     vmax=matrix.values.max())
        else:
            norm = None
        im = ax.imshow(matrix.values, aspect="auto", cmap="viridis", norm=norm)
        ax.set_xticks(range(len(matrix.columns)))
        ax.set_xticklabels(matrix.columns, rotation=70, ha="right", fontsize=8)
        ax.set_yticks(range(len(matrix.index)))
        ax.set_yticklabels(matrix.index, fontsize=8)
        ax.set_title(title)
        fig.colorbar(im, ax=ax, label="count per strain")
        fig.tight_layout()
        for ext in ("pdf", "png"):
            fig.savefig(a.outdir / f"{fname}.{ext}", dpi=200)
        plt.close(fig)

    heatmap(cps, "Figure 4a — PPV groups per strain by COG function",
             "fig4a_heatmap_count_per_strain", log_norm=True)

    pct = fc.copy()
    pct["pct"] = pct["n"] / pct.groupby("genus")["n"].transform("sum") * 100
    pct_mat = pct.pivot(index="genus", columns="COG_function", values="pct").fillna(0)
    pct_mat = pct_mat.reindex([g for g in ordered if g in pct_mat.index])
    pct_mat = pct_mat.reindex(columns=[c for c in COG_FUNCTION_ORDER if c in pct_mat.columns])
    heatmap(pct_mat, "Figure 4b — % of PPV gene groups by COG function within genus",
             "fig4b_heatmap_pct", log_norm=False)

    print(f"wrote figures to {a.outdir}")


if __name__ == "__main__":
    main()

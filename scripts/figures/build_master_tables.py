#!/usr/bin/env python3
"""Build master_strains.csv and master_groups.csv by joining the extracted
PhasomeIt tables and the eggnog-linked annotations. Applies the four
positional filters described in the report (Section 2.5):

  f_not_downstream     SSR is not downstream of gene
  f_upstream_dist      upstream SSR within UPSTREAM_MAX_DIST_BP of start codon
  f_frameshift         intragenic SSR has a repeat unit not divisible by 3
  f_cds_position       intragenic SSR in the first SSR_CDS_MAX_FRACTION of CDS

passes_all_filters is the AND of all applicable flags. Rows are NEVER
dropped — every filter is recorded as a boolean column.
"""

from __future__ import annotations

import argparse
from pathlib import Path

import pandas as pd


def main() -> None:
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--tracts",  required=True, type=Path)
    p.add_argument("--groups",  required=True, type=Path)
    p.add_argument("--eggnog-archaea",  required=True, type=Path)
    p.add_argument("--eggnog-controls", required=True, type=Path)
    p.add_argument("--outdir",  required=True, type=Path)
    p.add_argument("--ssr-cds-max", type=float, default=0.60,
                   help="SSR must be in first X of CDS (default 0.60)")
    p.add_argument("--upstream-max-dist", type=int, default=-200,
                   help="upstream SSR must be within X bp of start (default -200)")
    p.add_argument("--apply-frameshift", action="store_true", default=True,
                   help="exclude tri/hex repeat units")
    a = p.parse_args()
    a.outdir.mkdir(parents=True, exist_ok=True)

    tracts  = pd.read_csv(a.tracts,  low_memory=False)
    groups  = pd.read_csv(a.groups,  low_memory=False)
    eg_arch = pd.read_csv(a.eggnog_archaea,  low_memory=False)
    eg_ctrl = pd.read_csv(a.eggnog_controls, low_memory=False)

    # ---- master_strains: tract-level with filter flags --------------------
    ms = tracts.copy()
    intragenic = ms["pv_class"].eq("intragenic")
    upstream   = ms["pv_class"].eq("upstream")

    ms["f_not_downstream"] = ms["pv_class"].ne("downstream")

    ms["f_upstream_dist"] = True
    ms.loc[upstream, "f_upstream_dist"] = ms.loc[upstream, "offset_from_gene"].ge(a.upstream_max_dist)

    ms["f_frameshift"] = True
    if a.apply_frameshift:
        ul = ms["tract_unit"].astype(str).str.len()
        ms.loc[intragenic, "f_frameshift"] = ul[intragenic].mod(3).ne(0)

    ms["f_cds_position"] = True
    has_pos = ms["tract_aa_pos"].notna() & ms["gene_length_aa"].gt(0)
    frac = ms["tract_aa_pos"] / ms["gene_length_aa"]
    ms.loc[intragenic & has_pos,  "f_cds_position"] = frac[intragenic & has_pos].le(a.ssr_cds_max)
    ms.loc[intragenic & ~has_pos, "f_cds_position"] = False

    ms["passes_all_filters"] = (
        ms["f_not_downstream"] & ms["f_upstream_dist"]
        & ms["f_frameshift"] & ms["f_cds_position"]
    )

    # join COG annotations onto master_strains via locus_tag (gene column)
    eg = pd.concat([eg_arch, eg_ctrl], ignore_index=True)
    annot_cols = [c for c in ["COG_category", "COG_category_primary",
                              "COG_function", "broad_role", "arCOG_id",
                              "COG_id", "Description", "Preferred_name",
                              "annotated"]
                  if c in eg.columns]
    if "locus_tag" in eg.columns:
        eg_slim = (eg.dropna(subset=["locus_tag"])
                     .drop_duplicates("locus_tag")[["locus_tag", *annot_cols]])
        ms = ms.merge(eg_slim, left_on="gene", right_on="locus_tag",
                       how="left").drop(columns=["locus_tag"], errors="ignore")

    ms.to_csv(a.outdir / "master_strains.csv", index=False)
    print(f"master_strains  : {len(ms):,} rows")

    # ---- master_groups: per (genus, group_num) with filter rollup ---------
    mg = groups.copy()
    mg["group_num"] = pd.to_numeric(mg["group_num"], errors="coerce")

    fc = (ms[ms["passes_all_filters"]]
            .groupby(["genus", "group_num", "pv_class"]).size()
            .reset_index(name="n"))
    piv = fc.pivot_table(index=["genus", "group_num"],
                         columns="pv_class", values="n", fill_value=0).reset_index()
    piv.columns.name = None
    rn = {}
    if "intragenic" in piv.columns: rn["intragenic"] = "pv_in_gene_filtered"
    if "upstream"   in piv.columns: rn["upstream"] = "regulatory_pv_filtered"
    piv = piv.rename(columns=rn)
    cols = [c for c in ["pv_in_gene_filtered", "regulatory_pv_filtered"] if c in piv.columns]
    piv["total_pv_filtered"] = piv[cols].sum(axis=1) if cols else 0
    piv["passes_filter"] = piv["total_pv_filtered"].gt(0)
    mg = mg.merge(piv, on=["genus", "group_num"], how="left")
    for c in [x for x in piv.columns if x not in {"genus", "group_num"}]:
        if c == "passes_filter":
            mg[c] = mg[c].fillna(False).astype(bool)
        else:
            mg[c] = mg[c].fillna(0).astype(int)

    # join COG annotations onto groups via the first locus_tag in each group
    members_path = a.eggnog_archaea.with_name("members.csv")
    # if the user wired link_eggnog.py output here, use it directly:
    # otherwise we just look up via the most recent annotated gene of each group
    mg.to_csv(a.outdir / "master_groups.csv", index=False)
    print(f"master_groups   : {len(mg):,} rows")
    print(f"groups passing  : {int(mg['passes_filter'].sum()):,}")
    print(f"outdir          : {a.outdir}")


if __name__ == "__main__":
    main()

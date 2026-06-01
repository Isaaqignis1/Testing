#!/usr/bin/env python3
"""build the per-genus group summary table consumed by build_master_tables.

Reads phasomeit_group_summary.csv (output of parse_groups.py) for the
archaeal and bacterial sets and emits phasomeit_full_group_summary.csv
with one row per (genus, group_num), adding a `domain` column.
"""

from __future__ import annotations

import argparse
from pathlib import Path

import pandas as pd


KEEP = ["genus", "file", "page_title", "group_num", "group_name",
        "likely_function", "pv_in_gene", "total_pv", "total_genes",
        "regulatory_pv", "has_intragenic_pv", "has_any_pv"]


def load_one(p: Path, domain: str, skip_empty: bool) -> pd.DataFrame:
    if not p.exists():
        if skip_empty:
            print(f"  skip missing: {p}")
            return pd.DataFrame(columns=KEEP + ["domain"])
        raise SystemExit(f"missing: {p}")
    df = pd.read_csv(p)
    if df.empty and skip_empty:
        print(f"  skip empty:   {p}")
        return pd.DataFrame(columns=KEEP + ["domain"])
    keep = [c for c in KEEP if c in df.columns]
    df = df[keep].copy(); df["domain"] = domain
    return df


def main() -> None:
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--archaea",     required=True, type=Path)
    p.add_argument("--bacteria",    required=True, type=Path)
    p.add_argument("--out",         required=True, type=Path)
    p.add_argument("--skip-empty",  action="store_true")
    a = p.parse_args()
    a.out.parent.mkdir(parents=True, exist_ok=True)

    arch = load_one(a.archaea,  "Archaea",  a.skip_empty)
    bact = load_one(a.bacteria, "Bacteria", a.skip_empty)
    out = pd.concat([arch, bact], ignore_index=True)
    out.to_csv(a.out, index=False)
    print(f"wrote {a.out} ({len(out):,} rows)")
    for d, g in out.groupby("domain"):
        print(f"  {d}: {len(g):,} groups, {g['genus'].nunique()} genera")


if __name__ == "__main__":
    main()

#!/usr/bin/env python3
"""build the per-genus group summary table consumed by build_master_tables.ipynb

Reads phasomeit_group_summary.csv (output of parse_groups.py) and emits
phasomeit_full_group_summary.csv with one row per (genus, group_num),
dropping the dynamic strain_* columns and adding a `domain` column so the
notebook can concatenate archaeal and bacterial sets.

Two inputs are required (one per domain). If you have only one, pass the
same path twice with --skip-empty.
"""

from __future__ import annotations

import argparse
from pathlib import Path

import pandas as pd


KEEP_COLS = [
    "genus", "file", "page_title",
    "group_num", "group_name", "likely_function",
    "pv_in_gene", "total_pv", "total_genes",
    "regulatory_pv", "has_intragenic_pv", "has_any_pv",
]


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--archaea", required=True, type=Path,
                   help="phasomeit_group_summary.csv for the archaeal run")
    p.add_argument("--bacteria", required=True, type=Path,
                   help="phasomeit_group_summary.csv for the bacterial controls")
    p.add_argument("--out", required=True, type=Path)
    p.add_argument("--skip-empty", action="store_true",
                   help="silently treat a missing or empty input file as no rows")
    return p.parse_args()


def load_one(path: Path, domain: str, skip_empty: bool) -> pd.DataFrame:
    if not path.exists():
        if skip_empty:
            print(f"  skipped (missing): {path}")
            return pd.DataFrame(columns=KEEP_COLS + ["domain"])
        raise SystemExit(f"missing input: {path}")
    df = pd.read_csv(path)
    if df.empty and skip_empty:
        print(f"  skipped (empty):   {path}")
        return pd.DataFrame(columns=KEEP_COLS + ["domain"])
    keep = [c for c in KEEP_COLS if c in df.columns]
    df = df[keep].copy()
    df["domain"] = domain
    return df


def main() -> None:
    args = parse_args()
    args.out.parent.mkdir(parents=True, exist_ok=True)

    print("loading inputs")
    arch = load_one(args.archaea,  "Archaea",  args.skip_empty)
    bact = load_one(args.bacteria, "Bacteria", args.skip_empty)

    combined = pd.concat([arch, bact], ignore_index=True)
    combined.to_csv(args.out, index=False)

    print(f"wrote {args.out} ({len(combined):,} rows)")
    for d, g in combined.groupby("domain"):
        print(f"  {d}: {len(g):,} groups across {g['genus'].nunique()} genera")


if __name__ == "__main__":
    main()

#!/usr/bin/env python3
"""join phasomeit pv members to eggnog annotations by locus tag."""

from __future__ import annotations

import argparse
from pathlib import Path

import pandas as pd


ANNOT = ["COG_category", "COG_category_primary", "COG_function", "broad_role",
         "arCOG_id", "COG_id", "Description", "Preferred_name", "annotated"]


def main() -> None:
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--members", required=True, type=Path)
    p.add_argument("--eggnog",  required=True, type=Path)
    p.add_argument("--out",     required=True, type=Path)
    a = p.parse_args()
    a.out.parent.mkdir(parents=True, exist_ok=True)

    members = pd.read_csv(a.members)
    eg = pd.read_csv(a.eggnog, sep="\t")
    slim = (eg.dropna(subset=["locus_tag"])
              .drop_duplicates(subset=["locus_tag"])
              [["locus_tag", *ANNOT]])
    joined = members.merge(slim, on="locus_tag", how="left")
    joined.to_csv(a.out, index=False)
    matched = joined["annotated"].fillna(False).astype(bool).sum()
    print(f"wrote {len(joined):,} rows to {a.out}")
    print(f"  matched   {matched:,}")
    print(f"  unmatched {len(joined) - matched:,}")


if __name__ == "__main__":
    main()

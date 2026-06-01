#!/usr/bin/env python3
"""PV count per fasta across archaeal main + bacterial control sets."""

from __future__ import annotations

import argparse
import re
from pathlib import Path

import pandas as pd


def load_mapping(genus: str, gbk_dir: Path) -> dict[str, str]:
    f = gbk_dir / genus / f"{genus}_gbks.txt"
    if not f.exists():
        raise FileNotFoundError(f"missing mapping: {f}")
    mp = {}
    for i, line in enumerate(f.read_text().splitlines(), start=1):
        p = line.strip()
        if not p: continue
        mp[f"strain_{i}"] = Path(p).stem
    return mp


def main() -> None:
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--main-base",    required=True, type=Path)
    p.add_argument("--main-csv",     required=True, type=Path)
    p.add_argument("--control-base", required=True, type=Path)
    p.add_argument("--control-csv",  required=True, type=Path)
    p.add_argument("--out",          required=True, type=Path)
    a = p.parse_args()
    a.out.parent.mkdir(parents=True, exist_ok=True)

    datasets = [
        {"source": "main",    "gbk_dir": a.main_base,    "csv": a.main_csv},
        {"source": "control", "gbk_dir": a.control_base, "csv": a.control_csv},
    ]
    results = []
    for d in datasets:
        if not d["csv"].exists():
            print(f"skipping {d['source']}: missing {d['csv']}"); continue
        df = pd.read_csv(d["csv"])
        pv_cols = [c for c in df.columns if c.endswith("_pv_locus")]
        for genus, sub in df.groupby("genus"):
            try: mp = load_mapping(genus, d["gbk_dir"])
            except FileNotFoundError as e:
                print(f"  {e}"); continue
            for col in pv_cols:
                m = re.search(r"strain_(\d+)_pv_locus$", col)
                if not m: continue
                k = f"strain_{m.group(1)}"
                if k not in mp: continue
                results.append({
                    "source": d["source"], "genus": genus,
                    "fasta_name": mp[k],
                    "PV_count": int(sub[col].notna().sum()),
                })
    pd.DataFrame(results).to_csv(a.out, index=False)
    print(f"wrote {len(results):,} rows to {a.out}")


if __name__ == "__main__":
    main()

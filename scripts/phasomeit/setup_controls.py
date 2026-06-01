#!/usr/bin/env python3
"""stage per-genus phasomeit input folders for the bacterial controls.

Reads a controls species TSV (sample_accession, scientific_name) and groups
each sample's Prokka .gbk into a per-genus directory under --out-root.
"""

from __future__ import annotations

import argparse
import shutil
from pathlib import Path

import pandas as pd


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--prokka-dir", required=True, type=Path)
    p.add_argument("--mapping",    required=True, type=Path,
                   help="tsv with sample_accession + scientific_name")
    p.add_argument("--out-root",   required=True, type=Path)
    return p.parse_args()


def detect_sep(p: Path) -> str:
    return "\t" if p.suffix.lower() in {".tsv", ".tab"} else ","


def main() -> None:
    a = parse_args()
    a.out_root.mkdir(parents=True, exist_ok=True)

    df = pd.read_csv(a.mapping, sep=detect_sep(a.mapping)).fillna("")
    if "sample_accession" not in df.columns:
        for k in ("sample_id", "Sample", "accession"):
            if k in df.columns:
                df = df.rename(columns={k: "sample_accession"})
                break
    if "scientific_name" not in df.columns:
        for k in ("Species", "species"):
            if k in df.columns:
                df = df.rename(columns={k: "scientific_name"})
                break
    if not {"sample_accession", "scientific_name"}.issubset(df.columns):
        raise SystemExit("mapping needs sample_accession + scientific_name")

    df["genus"] = df["scientific_name"].str.split(n=1).str[0]
    total_copy = total_miss = 0

    for genus, sub in df.groupby("genus"):
        gdir = a.out_root / genus
        gdir.mkdir(parents=True, exist_ok=True)
        gbk_list = gdir / f"{genus}_gbks.txt"
        miss = []
        paths = []
        for _, row in sub.iterrows():
            sid = row["sample_accession"]
            src = a.prokka_dir / sid / f"{sid}.gbk"
            paths.append(str(src))
            if not src.exists():
                miss.append(str(src))
                continue
            shutil.copy2(src, gdir / src.name)
            total_copy += 1
        gbk_list.write_text("\n".join(paths) + "\n")
        if miss:
            (gdir / "missing_gbks.txt").write_text("\n".join(miss) + "\n")
            total_miss += len(miss)
        print(f"{genus}: listed={len(sub)} missing={len(miss)} copied={len(sub) - len(miss)}")

    print(f"\ntotal copied {total_copy}, missing {total_miss}")


if __name__ == "__main__":
    main()

#!/usr/bin/env python3
"""Produce per-genus <Genus>_gffs.txt lists from a species-calls table.

For each retained sample, locates the Prokka <sample_id>/<sample_id>.gff and
groups the paths by genus. Samples without a matching Prokka output go to
missing_samples.tsv (rather than being silently dropped).
"""

from __future__ import annotations

import argparse
from collections import defaultdict
from pathlib import Path

import pandas as pd


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--species-calls", required=True, type=Path)
    p.add_argument("--prokka-dir",    required=True, type=Path)
    p.add_argument("--out-dir",       required=True, type=Path)
    return p.parse_args()


def main() -> None:
    a = parse_args()
    a.out_dir.mkdir(parents=True, exist_ok=True)

    calls = pd.read_csv(a.species_calls)
    needed = {"sample_id", "genus"}
    if not needed.issubset(calls.columns):
        raise SystemExit(f"{a.species_calls}: need columns {needed}")

    grouped: dict[str, list[str]] = defaultdict(list)
    missing: list[tuple[str, str]] = []

    for _, row in calls.iterrows():
        sid = str(row["sample_id"]).strip()
        genus = str(row["genus"]).strip()
        if not sid or not genus:
            continue
        gff = a.prokka_dir / sid / f"{sid}.gff"
        if gff.exists():
            grouped[genus].append(str(gff))
        else:
            missing.append((sid, genus))

    for genus, paths in sorted(grouped.items()):
        out = a.out_dir / f"{genus}_gffs.txt"
        out.write_text("\n".join(paths) + "\n", encoding="utf-8")
        print(f"{genus}: {len(paths)} gffs -> {out.name}")

    if missing:
        miss_path = a.out_dir / "missing_samples.tsv"
        with miss_path.open("w") as fh:
            fh.write("sample_id\tgenus\n")
            for s, g in missing:
                fh.write(f"{s}\t{g}\n")
        print(f"\nmissing prokka outputs: {len(missing)} -> {miss_path.name}")

    print(f"\noutput dir: {a.out_dir}")


if __name__ == "__main__":
    main()

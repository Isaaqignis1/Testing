#!/usr/bin/env python3
"""produce per-genus <Genus>_gffs.txt lists from a species-calls table

Reads SpeciesCallsArchaea.csv (output of build_species_calls.py) and the
Prokka output root. For each retained sample, locates the matching
<sample_id>/<sample_id>.gff and groups the paths by genus. One text file
per genus is written, which is then consumed by phasomeit/setup_archaea.py.

Samples without a corresponding Prokka .gff are written to a missing_samples
file rather than silently dropped.
"""

from __future__ import annotations

import argparse
from collections import defaultdict
from pathlib import Path

import pandas as pd


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--species-calls", required=True, type=Path,
                   help="SpeciesCallsArchaea.csv produced by "
                        "build_species_calls.py")
    p.add_argument("--prokka-dir", required=True, type=Path,
                   help="Prokka output root, one subdir per sample_id")
    p.add_argument("--out-dir", required=True, type=Path,
                   help="directory to write <Genus>_gffs.txt files into")
    return p.parse_args()


def main() -> None:
    args = parse_args()
    args.out_dir.mkdir(parents=True, exist_ok=True)

    calls = pd.read_csv(args.species_calls)
    needed = {"sample_id", "genus"}
    if not needed.issubset(calls.columns):
        raise SystemExit(
            f"{args.species_calls}: need columns {needed}; "
            f"found {list(calls.columns)}"
        )

    grouped: dict[str, list[str]] = defaultdict(list)
    missing: list[tuple[str, str]] = []

    for _, row in calls.iterrows():
        sid = str(row["sample_id"]).strip()
        genus = str(row["genus"]).strip()
        if not sid or not genus:
            continue
        gff = args.prokka_dir / sid / f"{sid}.gff"
        if gff.exists():
            grouped[genus].append(str(gff))
        else:
            missing.append((sid, genus))

    for genus, paths in sorted(grouped.items()):
        out = args.out_dir / f"{genus}_gffs.txt"
        out.write_text("\n".join(paths) + "\n", encoding="utf-8")
        print(f"{genus}: {len(paths)} gffs -> {out.name}")

    if missing:
        miss_path = args.out_dir / "missing_samples.tsv"
        with miss_path.open("w") as fh:
            fh.write("sample_id\tgenus\n")
            for sid, genus in missing:
                fh.write(f"{sid}\t{genus}\n")
        print(f"\nmissing prokka outputs: {len(missing)} -> {miss_path.name}")
    else:
        print("\nall samples have a prokka .gff")

    print(f"\noutput dir: {args.out_dir}")


if __name__ == "__main__":
    main()

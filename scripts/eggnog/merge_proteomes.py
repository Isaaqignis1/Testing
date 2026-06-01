#!/usr/bin/env python3
"""Concatenate per-genome Prokka .faa files into one FASTA per domain,
prepending the sample id to every header so proteins trace back to their
genome."""

from __future__ import annotations

import argparse
from pathlib import Path


def merge(prokka_root: Path, out_fa: Path) -> int:
    out_fa.parent.mkdir(parents=True, exist_ok=True)
    n = 0
    with out_fa.open("w") as fout:
        for faa in sorted(prokka_root.glob("*/*.faa")):
            sample = faa.parent.name
            with faa.open() as fh:
                for line in fh:
                    if line.startswith(">"):
                        fout.write(f">{sample}_{line[1:]}")
                        n += 1
                    else:
                        fout.write(line)
    return n


def main() -> None:
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--prokka-archaea",  required=True, type=Path)
    p.add_argument("--prokka-controls", required=True, type=Path)
    p.add_argument("--out-archaea",     required=True, type=Path)
    p.add_argument("--out-controls",    required=True, type=Path)
    a = p.parse_args()

    n_arch = merge(a.prokka_archaea,  a.out_archaea)
    n_bact = merge(a.prokka_controls, a.out_controls)
    print(f"archaea : {n_arch:,} sequences -> {a.out_archaea}")
    print(f"controls: {n_bact:,} sequences -> {a.out_controls}")


if __name__ == "__main__":
    main()

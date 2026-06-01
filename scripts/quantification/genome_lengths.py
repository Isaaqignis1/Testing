#!/usr/bin/env python3
"""per-genome size and contig count from fasta files under a genus tree."""

from __future__ import annotations

import argparse
import csv
from pathlib import Path


def fasta_len(p: Path) -> tuple[int, int]:
    total = 0; seqs = 0
    with p.open() as fh:
        for line in fh:
            line = line.strip()
            if not line: continue
            if line.startswith(">"): seqs += 1
            else: total += len(line)
    return total, seqs


def main() -> None:
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--base", required=True, type=Path)
    p.add_argument("--out",  required=True, type=Path)
    a = p.parse_args()
    a.out.parent.mkdir(parents=True, exist_ok=True)

    rows = []
    for gd in sorted(a.base.iterdir()):
        if not gd.is_dir(): continue
        genus = gd.name
        # accept .fasta/.fna/.fa
        for fa in sorted(list(gd.glob("*.fasta")) + list(gd.glob("*.fna")) + list(gd.glob("*.fa"))):
            bp, n = fasta_len(fa)
            rows.append({
                "fasta_name": fa.stem, "genus": genus,
                "genome_size_bp": bp, "genome_size_mb": bp / 1_000_000,
                "n_contigs": n, "fasta_path": str(fa),
            })

    with a.out.open("w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=["fasta_name", "genus", "genome_size_bp",
                                            "genome_size_mb", "n_contigs", "fasta_path"])
        w.writeheader(); w.writerows(rows)
    print(f"wrote {len(rows):,} rows to {a.out}")


if __name__ == "__main__":
    main()

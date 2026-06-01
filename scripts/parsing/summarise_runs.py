#!/usr/bin/env python3
"""per-genus phasomeit run completeness summary."""

from __future__ import annotations

import argparse
import csv
import re
from pathlib import Path


def main() -> None:
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--root", required=True, type=Path)
    p.add_argument("--out",  required=True, type=Path)
    a = p.parse_args()

    rows = []
    for gd in sorted(p for p in a.root.iterdir() if p.is_dir()):
        genus = gd.name
        n_genomes = len(list(gd.glob("*.gbk")))
        sd = gd / "summary_tracts"
        ok = sd.exists()
        n_strains = len(list((sd / "strains").glob("*.html"))) if (sd / "strains").exists() else 0
        n_html = 0
        n_unique = 0
        gdir = sd / "groups"
        if gdir.exists():
            skip = {"index.html", "tracts.html", "Not in group.html"}
            files = [f.name for f in gdir.glob("*.html") if f.name not in skip]
            n_html = len(files)
            uniq = set()
            for n in files:
                gid = n[:-5].split("_")[0]
                if re.fullmatch(r"\d+", gid): uniq.add(gid)
            n_unique = len(uniq)
        rows.append({
            "genus": genus, "run_success": ok,
            "n_genomes_gbk": n_genomes,
            "n_strain_pages": n_strains,
            "n_group_html_files": n_html,
            "n_unique_groups": n_unique,
        })

    a.out.parent.mkdir(parents=True, exist_ok=True)
    with a.out.open("w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=["genus", "run_success", "n_genomes_gbk",
                                            "n_strain_pages", "n_group_html_files",
                                            "n_unique_groups"])
        w.writeheader(); w.writerows(rows)
    print(f"wrote {a.out}")


if __name__ == "__main__":
    main()

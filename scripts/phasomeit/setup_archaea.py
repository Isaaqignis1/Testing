#!/usr/bin/env python3
"""stage per-genus phasomeit input folders from <Genus>_gffs.txt lists.

Copies each genus's Prokka .gbk files into a per-genus directory ready for
PhasomeIt. The .gff lists produced by build_genus_gff_lists.py are converted
to .gbk paths (Prokka writes both).
"""

from __future__ import annotations

import argparse
import re
import shutil
from pathlib import Path


def genus_from_filename(p: Path) -> str:
    n = p.stem
    if n.endswith("_gffs"):
        return n[:-5]
    return re.sub(r"_gff(s)?$", "", n)


def gff_to_gbk(s: str) -> str:
    return re.sub(r"\.gff$", ".gbk", s)


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--in-dir",   required=True, type=Path)
    p.add_argument("--out-root", required=True, type=Path)
    return p.parse_args()


def main() -> None:
    a = parse_args()
    if not a.in_dir.exists():
        raise SystemExit(f"input folder not found: {a.in_dir}")
    a.out_root.mkdir(parents=True, exist_ok=True)

    txts = sorted(a.in_dir.glob("*_gffs.txt"))
    if not txts:
        raise SystemExit(f"no '*_gffs.txt' files in {a.in_dir}")

    total_copy = total_miss = 0
    for txt in txts:
        genus = genus_from_filename(txt)
        gdir = a.out_root / genus
        gdir.mkdir(parents=True, exist_ok=True)

        gbk_list = gdir / f"{genus}_gbks.txt"
        miss = []
        copied = 0
        paths = []
        for ln in txt.read_text().splitlines():
            ln = ln.strip()
            if not ln:
                continue
            paths.append(gff_to_gbk(ln))

        gbk_list.write_text("\n".join(paths) + ("\n" if paths else ""))
        for p in paths:
            src = Path(p)
            if not src.exists():
                miss.append(str(src))
                continue
            shutil.copy2(src, gdir / src.name)
            copied += 1
        if miss:
            (gdir / "missing_gbks.txt").write_text("\n".join(miss) + "\n")
        total_copy += copied
        total_miss += len(miss)
        print(f"{genus}: listed={len(paths)} copied={copied} missing={len(miss)}")

    print(f"\ntotal copied {total_copy}, missing {total_miss}")
    print(f"out: {a.out_root}")


if __name__ == "__main__":
    main()

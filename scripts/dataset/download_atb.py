#!/usr/bin/env python3
"""Download AllTheBacteria archaeal assemblies from OSF.

Two files, two URLs (both in the nehtp OSF node):
  1. species_calls.tsv.gz   (sample, species)   - default nypmt
  2. atb.archaea.assembly.202407.tar.xz  (all 818 fastas)  - default argkq

Flow:
  1. wget species_calls.tsv.gz                       (a few KB)
  2. wget the tarball into the cache                 (~154 MB)
  3. tar xf into --outdir-archaea (flat layout)
  4. Write metadata TSV (sample_id, scientific_name, genus)
  5. Write fasta_list.txt with the basenames

Idempotent: existing tarball / extracted files are kept.

Controls are NOT downloaded by this script. The OSF archaea node has no
controls. Supply controls separately if you want them; otherwise the
pipeline runs archaea-only.
"""

from __future__ import annotations

import argparse
import csv
import gzip
import subprocess
import sys
import tarfile
from pathlib import Path


DEFAULT_ARCHAEA_TARBALL_URL = "https://osf.io/download/argkq/"
DEFAULT_SPECIES_CALLS_URL   = "https://osf.io/download/nypmt/"

FASTA_EXTS = (".fa", ".fasta", ".fna",
              ".fa.gz", ".fasta.gz", ".fna.gz")


def log(m: str) -> None:
    print(m, flush=True)


def download(url: str, dest: Path) -> Path:
    """wget url -> dest, skip if already present and non-empty."""
    dest.parent.mkdir(parents=True, exist_ok=True)
    if dest.exists() and dest.stat().st_size > 0:
        log(f"  already cached: {dest.name}  ({dest.stat().st_size/1e6:.1f} MB)")
        return dest
    log(f"  downloading {url}")
    log(f"           -> {dest}")
    r = subprocess.run(
        ["wget", "--tries=3", "--timeout=300", "-O", str(dest), url],
        check=False,
    )
    if r.returncode != 0 or not dest.exists() or dest.stat().st_size == 0:
        raise SystemExit(f"download failed: {url}")
    return dest


def parse_args():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--outdir-archaea", required=True, type=Path,
                   help="where extracted FASTAs are placed (flat)")
    p.add_argument("--metadata-out",   required=True, type=Path,
                   help="TSV: sample_id, scientific_name, genus")
    p.add_argument("--fasta-list-out", required=True, type=Path,
                   help="text file listing extracted FASTA basenames")
    p.add_argument("--tar-cache",      required=True, type=Path,
                   help="where to keep the downloaded tarball + species_calls")
    p.add_argument("--archaea-tarball-url", default=DEFAULT_ARCHAEA_TARBALL_URL)
    p.add_argument("--species-calls-url",   default=DEFAULT_SPECIES_CALLS_URL)
    p.add_argument("--keep-tarball",   action="store_true",
                   help="keep the .tar.xz after extraction (default keeps it)")
    p.add_argument("--skip-download",  action="store_true",
                   help="for testing - assume files already in cache")
    return p.parse_args()


def main():
    a = parse_args()
    a.tar_cache.mkdir(parents=True, exist_ok=True)
    a.outdir_archaea.mkdir(parents=True, exist_ok=True)

    # 1. species_calls.tsv.gz
    species_path = a.tar_cache / "species_calls.tsv.gz"
    if not a.skip_download:
        download(a.species_calls_url, species_path)

    # 2. parse it
    rows = []
    with gzip.open(species_path, "rt") as fh:
        r = csv.DictReader(fh, delimiter="\t")
        for row in r:
            rows.append(row)
    log(f"  species_calls rows: {len(rows)}")

    # 3. metadata.tsv
    a.metadata_out.parent.mkdir(parents=True, exist_ok=True)
    n_meta = 0
    with a.metadata_out.open("w", newline="") as fh:
        w = csv.writer(fh, delimiter="\t")
        w.writerow(["sample_id", "scientific_name", "genus"])
        for row in rows:
            sample = (row.get("Sample") or row.get("sample") or "").strip()
            sn     = (row.get("Species") or row.get("species") or "").strip()
            if not sample or not sn:
                continue
            genus = sn.split(maxsplit=1)[0]
            w.writerow([sample, sn, genus])
            n_meta += 1
    log(f"  wrote {a.metadata_out} ({n_meta} samples)")

    # 4. tarball
    tar_path = a.tar_cache / "atb.archaea.assembly.tar.xz"
    if not a.skip_download:
        download(a.archaea_tarball_url, tar_path)

    # 5. extract, flattened
    log(f"  extracting {tar_path} -> {a.outdir_archaea}")
    n_extracted = 0
    with tarfile.open(tar_path, "r:xz") as tf:
        for ti in tf:
            if not ti.isfile():
                continue
            base = Path(ti.name).name
            low = base.lower()
            if not any(low.endswith(ext) for ext in FASTA_EXTS):
                continue
            # flatten path inside outdir
            ti.name = base
            tf.extract(ti, path=a.outdir_archaea)
            n_extracted += 1
    log(f"  extracted: {n_extracted} fasta files")

    # 6. fasta_list.txt
    a.fasta_list_out.parent.mkdir(parents=True, exist_ok=True)
    fastas = sorted(
        p.name for p in a.outdir_archaea.iterdir()
        if p.is_file() and not p.name.startswith(".")
    )
    a.fasta_list_out.write_text("\n".join(fastas) + ("\n" if fastas else ""))
    log(f"  wrote {a.fasta_list_out} ({len(fastas)} entries)")

    log("done.")


if __name__ == "__main__":
    main()

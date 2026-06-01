#!/usr/bin/env python3
"""Download AllTheBacteria assemblies for a target genus list (+ controls).

AllTheBacteria packages assemblies in xzipped tar archives hosted on OSF.
The index file file_list.all.latest.tsv.gz lists, for every sample, the
tarball it lives in and the path inside that tarball.

This script:
  1. Downloads the index file_list.all.latest.tsv.gz from OSF.
  2. Filters rows whose species_sylph genus is in the target genus list.
  3. For controls, picks the first N samples per scientific_name.
  4. Groups results by tar_xz_url, downloads each unique tarball into the
     cache, and extracts only the required FASTAs.
  5. Writes filtered metadata + fasta_list.txt + controls species map.

Idempotent: existing tarballs and extracted FASTAs are kept.
"""

from __future__ import annotations

import argparse
import csv
import gzip
import os
import shutil
import subprocess
import sys
import tarfile
from collections import defaultdict
from pathlib import Path
from urllib.request import urlopen, Request


DEFAULT_INDEX_URL = "https://osf.io/download/4yv85/"


def log(msg: str) -> None:
    print(msg, flush=True)


def download(url: str, dest: Path) -> Path:
    """Download url to dest if missing or zero-byte. Returns dest."""
    dest.parent.mkdir(parents=True, exist_ok=True)
    if dest.exists() and dest.stat().st_size > 0:
        return dest
    log(f"  downloading {url} -> {dest.name}")
    rc = subprocess.run(
        ["wget", "--quiet", "--tries=3", "--timeout=120",
         "-O", str(dest), url],
        check=False,
    ).returncode
    if rc != 0 or not dest.exists() or dest.stat().st_size == 0:
        try:
            req = Request(url, headers={"User-Agent": "archaea-pv-analysis/1.0"})
            with urlopen(req) as r, dest.open("wb") as out:
                shutil.copyfileobj(r, out)
        except Exception as exc:
            raise SystemExit(f"download failed: {url}\n  {exc}")
    return dest


def read_genus_list(p: Path) -> set:
    return {ln.strip() for ln in p.read_text().splitlines()
            if ln.strip() and not ln.strip().startswith("#")}


def read_controls_tsv(p: Path):
    sep = "\t" if p.suffix.lower() in {".tsv", ".tab"} else ","
    with p.open() as fh:
        return [r for r in csv.DictReader(fh, delimiter=sep)]


def stream_index(path: Path):
    with gzip.open(path, "rt") as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        for r in reader:
            yield r


def genus_of(species: str) -> str:
    return species.strip().split(maxsplit=1)[0] if species else ""


def extract_from_tarball(tarball: Path, members, outdir: Path):
    """Extract only the given member paths from tarball into outdir.
    Returns (n_extracted, n_missing)."""
    outdir.mkdir(parents=True, exist_ok=True)
    needed = set(members)
    n_ok = 0
    with tarfile.open(tarball, "r:xz") as tf:
        for ti in tf:
            name = ti.name
            base = os.path.basename(name)
            if name in needed:
                match = name
            elif base in needed:
                match = base
            else:
                continue
            ti.name = base
            tf.extract(ti, path=outdir)
            n_ok += 1
            needed.discard(match)
    return n_ok, len(needed)


def parse_args():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--archaea-genera-file", required=True, type=Path)
    p.add_argument("--controls-tsv",       type=Path, default=None)
    p.add_argument("--outdir-archaea",     required=True, type=Path)
    p.add_argument("--outdir-controls",    required=True, type=Path)
    p.add_argument("--metadata-out",       required=True, type=Path)
    p.add_argument("--controls-out",       required=True, type=Path)
    p.add_argument("--tar-cache",          required=True, type=Path)
    p.add_argument("--index-url",          default=DEFAULT_INDEX_URL)
    p.add_argument("--keep-tarballs",      action="store_true")
    p.add_argument("--skip-archaea",       action="store_true")
    p.add_argument("--skip-controls",      action="store_true")
    p.add_argument("--max-per-species",    type=int, default=None,
                   help="cap archaeal downloads per species, for quick tests")
    return p.parse_args()


def main():
    a = parse_args()
    a.outdir_archaea.mkdir(parents=True, exist_ok=True)
    a.outdir_controls.mkdir(parents=True, exist_ok=True)
    a.tar_cache.mkdir(parents=True, exist_ok=True)

    target_genera = read_genus_list(a.archaea_genera_file)
    log(f"target genera: {len(target_genera)}")

    # 1. ATB index
    index_path = a.tar_cache / "file_list.all.latest.tsv.gz"
    download(a.index_url, index_path)
    log(f"  index file: {index_path.stat().st_size/1e6:.1f} MB")

    # 2. filter
    arch_rows = []
    controls_buckets = defaultdict(list)
    controls_names = set()
    if a.controls_tsv and not a.skip_controls:
        for cspec in read_controls_tsv(a.controls_tsv):
            controls_names.add(cspec["scientific_name"])

    per_species_count = defaultdict(int)
    for r in stream_index(index_path):
        species = (r.get("species_sylph") or "").strip()
        sample = (r.get("sample") or "").strip()
        if not species or not sample:
            continue
        g = genus_of(species)
        if g in target_genera:
            if a.max_per_species:
                if per_species_count[species] >= a.max_per_species:
                    continue
                per_species_count[species] += 1
            arch_rows.append({**r, "scientific_name": species, "genus": g})
        if species in controls_names:
            controls_buckets[species].append({**r, "scientific_name": species, "genus": g})

    log(f"  archaeal samples matched: {len(arch_rows)}")
    for cn, rows in controls_buckets.items():
        log(f"  controls {cn}: {len(rows)} candidates")

    a.metadata_out.parent.mkdir(parents=True, exist_ok=True)
    with a.metadata_out.open("w", newline="") as fh:
        w = csv.writer(fh, delimiter="\t")
        w.writerow(["sample_id", "scientific_name", "genus", "tar_xz",
                    "filename_in_tar_xz", "tar_xz_url"])
        for r in arch_rows:
            w.writerow([r["sample"], r["scientific_name"], r["genus"],
                        r.get("tar_xz", ""), r.get("filename_in_tar_xz", ""),
                        r.get("tar_xz_url", "")])
    log(f"  wrote {a.metadata_out}")

    # 3. controls pick
    controls_picks = []
    if a.controls_tsv and not a.skip_controls:
        for cspec in read_controls_tsv(a.controls_tsv):
            sn = cspec["scientific_name"]
            need = int(cspec.get("min_strains", "8"))
            cands = controls_buckets.get(sn, [])
            if len(cands) < need:
                log(f"  WARN: {sn}: only {len(cands)} cands (wanted {need})")
            controls_picks.extend(cands[:need])

    if controls_picks:
        a.controls_out.parent.mkdir(parents=True, exist_ok=True)
        with a.controls_out.open("w", newline="") as fh:
            w = csv.writer(fh, delimiter="\t")
            w.writerow(["sample_accession", "scientific_name"])
            for r in controls_picks:
                w.writerow([r["sample"], r["scientific_name"]])
        log(f"  wrote {a.controls_out}")

    # 4. download + extract by tarball
    def fasta_basename(filename_in_tar: str) -> str:
        return Path(filename_in_tar).name

    def process_set(rows, outdir, list_path, skip, label):
        if skip or not rows:
            log(f"  skipping {label} ({len(rows)} rows)")
            return
        by_tar = defaultdict(list)
        for r in rows:
            url = r.get("tar_xz_url", "").strip()
            if url:
                by_tar[url].append(r)

        with list_path.open("w") as listfh:
            for url, group in by_tar.items():
                tar_name = group[0].get("tar_xz") or Path(url).name + ".tar.xz"
                tar_path = a.tar_cache / tar_name
                try:
                    download(url, tar_path)
                except SystemExit as e:
                    log(f"  {label}: tarball failed: {e}")
                    continue
                members = [r["filename_in_tar_xz"] for r in group]
                n_ok, n_miss = extract_from_tarball(tar_path, members, outdir)
                log(f"  {label}: {tar_name}: {n_ok}/{len(members)} (missing {n_miss})")
                for r in group:
                    bn = fasta_basename(r["filename_in_tar_xz"])
                    if (outdir / bn).exists():
                        listfh.write(bn + "\n")
                if not a.keep_tarballs:
                    try:
                        tar_path.unlink()
                    except FileNotFoundError:
                        pass
        log(f"  wrote {list_path}")

    process_set(arch_rows, a.outdir_archaea,
                a.outdir_archaea.parent / "fasta_list.txt",
                a.skip_archaea, "archaea")
    process_set(controls_picks, a.outdir_controls,
                a.outdir_controls.parent / "fasta_list.txt",
                a.skip_controls, "controls")

    log("done.")


if __name__ == "__main__":
    main()

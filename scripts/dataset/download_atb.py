#!/usr/bin/env python3
"""Download AllTheBacteria assemblies for a target genus list and controls.

The AllTheBacteria release tree exposes per-batch metadata and tarballs at
the EBI FTP mirror. This script:

  1. Downloads the species-call table (atb metadata).
  2. Filters rows whose Genus is in --archaea-genera-file.
  3. Resolves each sample's batch + filename and downloads the fasta.
  4. (Optional) Downloads bacterial controls listed in --controls-tsv —
     for each scientific_name, pulls the first N matching assemblies.
  5. Writes fasta_list.txt files used by the Prokka array jobs.

This is bandwidth-heavy. Plan for several GB and an hour+ on the first run.
Idempotent: existing files are skipped via size check.

Defaults are wired to the AllTheBacteria v0.2 release. If ATB moves, override
--atb-base-url.
"""

from __future__ import annotations

import argparse
import csv
import gzip
import io
import shutil
import sys
from collections import defaultdict
from pathlib import Path
from urllib.parse import urljoin
from urllib.request import urlopen, Request


DEFAULT_ATB_BASE = (
    "https://ftp.ebi.ac.uk/pub/databases/AllTheBacteria/Releases/0.2/"
)
DEFAULT_METADATA_URL = (
    DEFAULT_ATB_BASE + "metadata/species_calls.tsv.gz"
)
DEFAULT_ASSEMBLY_INDEX_URL = (
    DEFAULT_ATB_BASE + "assembly/file_list.all.latest.tsv.gz"
)


def log(msg: str) -> None:
    print(msg, flush=True)


def fetch(url: str, dest: Path, gzipped: bool = False) -> Path:
    """Download url to dest if not present; return dest path.
    If gzipped, also decompress to dest.with_suffix('') and return that."""
    dest.parent.mkdir(parents=True, exist_ok=True)
    if not dest.exists() or dest.stat().st_size == 0:
        log(f"  downloading {url}")
        req = Request(url, headers={"User-Agent": "archaea-pv-analysis/1.0"})
        with urlopen(req) as r, dest.open("wb") as out:
            shutil.copyfileobj(r, out)
    if gzipped:
        out = dest.with_suffix("")
        if not out.exists() or out.stat().st_size == 0:
            log(f"  decompressing {dest.name}")
            with gzip.open(dest, "rb") as g, out.open("wb") as f:
                shutil.copyfileobj(g, f)
        return out
    return dest


def read_genus_list(path: Path) -> set[str]:
    g = set()
    for ln in path.read_text().splitlines():
        ln = ln.strip()
        if ln and not ln.startswith("#"):
            g.add(ln)
    return g


def read_controls_tsv(path: Path) -> list[dict]:
    out: list[dict] = []
    sep = "\t" if path.suffix.lower() in {".tsv", ".tab"} else ","
    with path.open() as fh:
        reader = csv.DictReader(fh, delimiter=sep)
        for r in reader:
            r = {k.strip(): (v.strip() if v else "") for k, v in r.items()}
            out.append(r)
    return out


def load_species_calls(tsv: Path) -> list[dict]:
    """Yield rows from the ATB species_calls.tsv as dicts."""
    rows = []
    with tsv.open() as fh:
        reader = csv.reader(fh, delimiter="\t")
        header = next(reader)
        header = [h.strip() for h in header]
        for parts in reader:
            if not parts:
                continue
            d = dict(zip(header, parts))
            rows.append(d)
    return rows


def normalize_genus(scientific_name: str) -> str:
    parts = scientific_name.strip().split(maxsplit=1)
    return parts[0] if parts else ""


def build_atb_assembly_url(sample_id: str, file_index: dict[str, str],
                            base_url: str) -> str | None:
    """Resolve a sample to its fasta URL using the ATB file index."""
    rel = file_index.get(sample_id)
    if rel is None:
        return None
    return urljoin(base_url + "assembly/", rel)


def load_file_index(path: Path) -> dict[str, str]:
    """ATB file_list.tsv columns: sample_id\trelative_path"""
    idx: dict[str, str] = {}
    with path.open() as fh:
        for ln in fh:
            ln = ln.strip()
            if not ln:
                continue
            cols = ln.split("\t")
            if len(cols) < 2:
                continue
            sample, rel = cols[0], cols[1]
            idx[sample] = rel
    return idx


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--archaea-genera-file", required=True, type=Path,
                   help="text file, one GTDB genus per line")
    p.add_argument("--controls-tsv", type=Path, default=None,
                   help="TSV with scientific_name + min_strains columns")
    p.add_argument("--outdir-archaea", required=True, type=Path)
    p.add_argument("--outdir-controls", required=True, type=Path)
    p.add_argument("--metadata-out", required=True, type=Path,
                   help="path to write the filtered archaea metadata TSV")
    p.add_argument("--atb-base-url", default=DEFAULT_ATB_BASE)
    p.add_argument("--metadata-url", default=DEFAULT_METADATA_URL)
    p.add_argument("--assembly-index-url", default=DEFAULT_ASSEMBLY_INDEX_URL)
    p.add_argument("--skip-archaea", action="store_true")
    p.add_argument("--skip-controls", action="store_true")
    return p.parse_args()


def main() -> None:
    args = parse_args()
    args.outdir_archaea.mkdir(parents=True, exist_ok=True)
    args.outdir_controls.mkdir(parents=True, exist_ok=True)

    cache = args.outdir_archaea.parent / "_cache"
    cache.mkdir(exist_ok=True)

    # ---- 1. metadata + file index --------------------------------------
    log("[1/4] fetching ATB metadata and file index")
    meta_gz  = cache / "species_calls.tsv.gz"
    meta_tsv = fetch(args.metadata_url, meta_gz, gzipped=True)
    idx_gz   = cache / "file_list.tsv.gz"
    idx_tsv  = fetch(args.assembly_index_url, idx_gz, gzipped=True)
    file_index = load_file_index(idx_tsv)
    log(f"     file index: {len(file_index):,} samples")

    rows = load_species_calls(meta_tsv)
    log(f"     species calls: {len(rows):,} samples")

    # ---- 2. filter to target archaeal genera --------------------------
    target = read_genus_list(args.archaea_genera_file)
    log(f"[2/4] filtering to {len(target)} target genera")

    def get(row: dict, *keys: str) -> str:
        for k in keys:
            if k in row and row[k]:
                return row[k]
        return ""

    arch_rows: list[dict] = []
    for r in rows:
        sci = get(r, "Species", "scientific_name", "species")
        genus = normalize_genus(sci)
        if genus in target:
            sample = get(r, "Sample", "sample_id", "sample_accession", "accession")
            r["sample_id"] = sample
            r["scientific_name"] = sci
            r["genus"] = genus
            arch_rows.append(r)
    log(f"     {len(arch_rows):,} archaeal samples matched")

    args.metadata_out.parent.mkdir(parents=True, exist_ok=True)
    with args.metadata_out.open("w", newline="") as fh:
        w = csv.writer(fh, delimiter="\t")
        w.writerow(["sample_id", "scientific_name", "genus"])
        for r in arch_rows:
            w.writerow([r["sample_id"], r["scientific_name"], r["genus"]])
    log(f"     wrote {args.metadata_out}")

    # ---- 3. download archaeal fastas -----------------------------------
    if not args.skip_archaea:
        log("[3/4] downloading archaeal fastas")
        n_ok = n_skip = n_miss = 0
        list_path = args.outdir_archaea.parent / "fasta_list.txt"
        with list_path.open("w") as listfh:
            for r in arch_rows:
                sid = r["sample_id"]
                url = build_atb_assembly_url(sid, file_index, args.atb_base_url)
                if url is None:
                    n_miss += 1
                    continue
                fname = Path(url).name
                dest = args.outdir_archaea / fname
                try:
                    fetch(url, dest)
                    listfh.write(fname + "\n")
                    if dest.stat().st_size > 0:
                        n_ok += 1
                    else:
                        n_miss += 1
                except Exception as exc:
                    log(f"     {sid}: {exc}")
                    n_miss += 1
        log(f"     ok={n_ok}  missing={n_miss}")
        log(f"     fasta list: {list_path}")
    else:
        log("[3/4] skipping archaea (--skip-archaea)")

    # ---- 4. controls ---------------------------------------------------
    if args.skip_controls or args.controls_tsv is None:
        log("[4/4] skipping controls")
        return

    log("[4/4] downloading controls")
    controls_spec = read_controls_tsv(args.controls_tsv)
    by_species: dict[str, list[dict]] = defaultdict(list)
    for r in rows:
        sci = get(r, "Species", "scientific_name", "species")
        by_species[sci].append({**r, "sample_id":
                                 get(r, "Sample", "sample_id",
                                     "sample_accession", "accession"),
                                 "scientific_name": sci})

    ctrl_list_path = args.outdir_controls.parent / "fasta_list.txt"
    ctrl_species_path = args.outdir_controls.parent / "controls_species.tsv"
    with ctrl_list_path.open("w") as listfh, ctrl_species_path.open("w") as spfh:
        spfh.write("sample_accession\tscientific_name\n")
        for spec in controls_spec:
            sn = spec["scientific_name"]
            need = int(spec.get("min_strains", "8"))
            candidates = by_species.get(sn, [])
            if len(candidates) < need:
                log(f"     WARN: {sn}: only {len(candidates)} candidates "
                    f"(wanted {need})")
            picks = candidates[:need]
            for r in picks:
                sid = r["sample_id"]
                url = build_atb_assembly_url(sid, file_index, args.atb_base_url)
                if url is None:
                    continue
                fname = Path(url).name
                dest = args.outdir_controls / fname
                try:
                    fetch(url, dest)
                    listfh.write(fname + "\n")
                    spfh.write(f"{sid}\t{sn}\n")
                except Exception as exc:
                    log(f"     {sid}: {exc}")
            log(f"     {sn}: downloaded {len(picks)}")

    log(f"     fasta list: {ctrl_list_path}")
    log(f"     species map: {ctrl_species_path}")


if __name__ == "__main__":
    main()

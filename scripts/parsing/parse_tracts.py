#!/usr/bin/env python3
"""extract ssr tract rows from phasomeit html, combining archaeal and bacterial sets."""

from __future__ import annotations

import argparse
import csv
import re
import statistics
from collections import Counter
from html import unescape as html_unescape
from pathlib import Path


def strip_tags(h: str) -> str:
    h = re.sub(r"(?is)<script.*?>.*?</script>", " ", h)
    h = re.sub(r"(?is)<style.*?>.*?</style>",   " ", h)
    h = re.sub(r"(?is)<[^>]+>", " ", h)
    return re.sub(r"\s+", " ", html_unescape(h)).strip()


def safe_int(v):
    try: return int(str(v).strip())
    except (ValueError, TypeError): return None


def norm(s: str) -> str:
    return re.sub(r"[^a-z0-9]", "", s.lower())


def extract_strain_accession(h: str) -> str:
    m = re.search(r"Filepath:\s*([^\s<]+)", h)
    if m:
        stem = Path(m.group(1).strip()).stem
        if stem and stem.lower() != "strain": return stem
    m = re.search(r"Strain name:\s*([^<\n\r]+)", h)
    if m:
        n = m.group(1).strip()
        if n and n.lower() != "strain": return n
    return ""


def parse_tract_type(raw: str):
    m = re.match(r"^([A-Za-z]+)(\d+)$", raw.strip())
    if m: return m.group(1).upper(), int(m.group(2))
    return raw.strip(), None


def pv_class_from_offset(o):
    if o is None: return "unknown"
    if o == 0:    return "intragenic"
    if o < 0:     return "upstream"
    return "downstream"


def group_num_from_cell(c: str):
    m = re.search(r"groups/(\d+)\.html", c)
    if m: return int(m.group(1))
    return safe_int(strip_tags(c))


MEMBER_BLOCK_RE = re.compile(
    r'<p\s+id="([^"]+)">.*?white-space\s*:\s*pre[^>]*>(.*?)</div>',
    re.DOTALL | re.IGNORECASE)
LINE_OFS_RE = re.compile(r"^\+(\d+)\s+")


def count_residues(b: str) -> int:
    total = 0
    for line in b.split("<br />"):
        p = re.sub(r"<[^>]+>", "", line)
        p = html_unescape(p)
        p = LINE_OFS_RE.sub("", p)
        total += len(p.replace(" ", "").replace("\n", "").replace("\r", ""))
    return total


def extract_tract_aa_pos(b: str):
    glen = count_residues(b) or None
    pos = None
    for line in b.split("<br />"):
        if "#0033dd" not in line.lower(): continue
        plain = html_unescape(re.sub(r"<[^>]+>", "", line))
        m = LINE_OFS_RE.match(plain) or LINE_OFS_RE.match(plain.lstrip())
        ofs = int(m.group(1)) if m else 0
        before = re.split(r"<span[^>]*#0033dd", line, maxsplit=1)[0]
        before = html_unescape(re.sub(r"<[^>]+>", "", before))
        before = LINE_OFS_RE.sub("", before)
        pos = ofs + len(before.replace(" ", "")) + 1
        break
    return pos, glen


def parse_group_html(h: str) -> dict:
    out = {}
    for m in MEMBER_BLOCK_RE.finditer(h):
        locus, sb = m.groups()
        pos, glen = extract_tract_aa_pos(sb)
        out[locus] = {"tract_aa_pos": pos, "gene_length_aa_seq": glen}
    return out


def build_group_lookup(gd: Path) -> dict:
    gdir = gd / "summary_tracts" / "groups"
    if not gdir.is_dir(): return {}
    lookup = {}
    rx = re.compile(r"^(\d+)\.html$")
    for f in sorted(gdir.glob("*.html")):
        m = rx.match(f.name)
        if not m: continue
        gn = int(m.group(1))
        try: h = f.read_text(errors="ignore")
        except Exception: continue
        members = parse_group_html(h)
        if not members: continue
        positions = [v["tract_aa_pos"] for v in members.values() if v["tract_aa_pos"] is not None]
        lens      = [v["gene_length_aa_seq"] for v in members.values() if v["gene_length_aa_seq"] is not None]
        lookup[gn] = {
            "median_tract_aa_pos": statistics.median(positions) if positions else None,
            "median_gene_length_aa_seq": statistics.median(lens) if lens else None,
            "members": members,
        }
    return lookup


STRAIN_COLS = {
    "tractno": "tract_no", "contig": "contig", "location": "location_bp",
    "tract": "tract_type_raw", "onlength": "on_length_raw",
    "genegroup": "group_num_raw", "offsetfromgene": "offset_raw",
    "gene": "gene", "length": "gene_length_aa", "function": "function",
}
STRAIN_REQ = {"tractno", "tract", "location", "genegroup"}


def parse_strain_html(html, domain, genus, fname, glookup):
    rows = []
    acc = extract_strain_accession(html)
    tables = re.findall(r"(?is)<table\b.*?>.*?</table>", html)
    for tbl in tables:
        rr = re.findall(r"(?is)<tr\b.*?>.*?</tr>", tbl)
        if len(rr) < 2: continue
        hcells = re.findall(r"(?is)<t[dh]\b.*?>(.*?)</t[dh]>", rr[0])
        hn = [norm(strip_tags(c)) for c in hcells]
        if not STRAIN_REQ.issubset(set(hn)): continue
        cm = {STRAIN_COLS[n]: i for i, n in enumerate(hn) if n in STRAIN_COLS}
        for row in rr[1:]:
            craw = re.findall(r"(?is)<t[dh]\b.*?>(.*?)</t[dh]>", row)
            ctxt = [strip_tags(c) for c in craw]
            def gt(k):
                i = cm.get(k)
                if i is None or i >= len(ctxt): return ""
                return ctxt[i].strip()
            def gr(k):
                i = cm.get(k)
                if i is None or i >= len(craw): return ""
                return craw[i]
            raw_t = gt("tract_type_raw")
            unit, rlen = parse_tract_type(raw_t)
            ons = gt("on_length_raw")
            on = safe_int(ons) if ons.lower() != "none" else None
            gn = group_num_from_cell(gr("group_num_raw"))
            ofs = safe_int(gt("offset_raw"))
            cls = pv_class_from_offset(ofs)
            gl = safe_int(gt("gene_length_aa"))
            locus = gt("gene")
            tpos = gls = None
            if gn is not None and gn in glookup:
                g = glookup[gn]
                if locus in g["members"]:
                    tpos = g["members"][locus]["tract_aa_pos"]
                    gls  = g["members"][locus]["gene_length_aa_seq"]
                else:
                    tpos = g["median_tract_aa_pos"]
                    gls  = g["median_gene_length_aa_seq"]
            rows.append({
                "domain": domain, "genus": genus, "strain_file": fname,
                "strain_accession": acc,
                "tract_no": safe_int(gt("tract_no")), "contig": gt("contig"),
                "location_bp": safe_int(gt("location_bp")),
                "tract_type": raw_t, "tract_unit": unit, "tract_repeat_len": rlen,
                "on_length": on, "group_num": gn,
                "offset_from_gene": ofs, "pv_class": cls,
                "gene": locus, "gene_length_aa": gl,
                "function": gt("function"),
                "tract_aa_pos": tpos, "gene_length_aa_seq": gls,
            })
        if rows: break
    return rows


def crawl(base: Path, domain: str, all_rows: list) -> tuple[int, int]:
    if not base.exists():
        print(f"  base not found: {base}")
        return 0, 0
    ng = nf = 0
    for gd in sorted(p for p in base.iterdir() if p.is_dir()):
        sd = gd / "summary_tracts" / "strains"
        if not sd.is_dir(): continue
        genus = gd.name
        print(f"  {genus} ...", end=" ", flush=True)
        gl = build_group_lookup(gd)
        sf = sorted(sd.glob("*.html"))
        gr = 0
        for hf in sf:
            try: html = hf.read_text(errors="ignore")
            except Exception as e:
                print(f"\n    cannot read {hf}: {e}"); continue
            rs = parse_strain_html(html, domain, genus, hf.name, gl)
            all_rows.extend(rs)
            gr += len(rs); nf += 1
        ng += 1
        print(f"{len(sf)} files, {gr:,} tracts, {len(gl)} groups indexed")
    return ng, nf


FIELDS = [
    "domain", "genus", "strain_file", "strain_accession",
    "tract_no", "contig", "location_bp",
    "tract_type", "tract_unit", "tract_repeat_len", "on_length",
    "group_num", "offset_from_gene", "pv_class",
    "gene", "gene_length_aa", "function",
    "tract_aa_pos", "gene_length_aa_seq",
]


def main() -> None:
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--archaea-base",  required=True, type=Path)
    p.add_argument("--controls-base", required=True, type=Path)
    p.add_argument("--outdir",        required=True, type=Path)
    a = p.parse_args()
    a.outdir.mkdir(parents=True, exist_ok=True)

    rows = []
    print("archaea");  ag, af = crawl(a.archaea_base,  "Archaea",  rows)
    print("\nbacteria"); cg, cf = crawl(a.controls_base, "Bacteria", rows)

    out = a.outdir / "phasomeit_tract_data.csv"
    with out.open("w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=FIELDS, extrasaction="ignore")
        w.writeheader(); w.writerows(rows)
    print(f"\nwrote {out}")
    print(f"total rows {len(rows):,}")
    print(f"  archaea  {sum(1 for r in rows if r['domain']=='Archaea'):,}")
    print(f"  bacteria {sum(1 for r in rows if r['domain']=='Bacteria'):,}")
    for k, v in sorted(Counter(r["pv_class"] for r in rows).items()):
        n_pos = sum(1 for r in rows if r["pv_class"] == k and r["tract_aa_pos"] is not None)
        print(f"  {k:15s} {v:6,}  (with pos: {n_pos:,})")


if __name__ == "__main__":
    main()

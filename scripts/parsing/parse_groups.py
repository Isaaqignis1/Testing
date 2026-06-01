#!/usr/bin/env python3
"""parse phasomeit group html pages into group_summary, annotated_functions,
members and pairwise CSVs."""

from __future__ import annotations

import argparse
import csv
import re
from html import unescape
from pathlib import Path


SKIP_FILES = {"index.html", "tracts.html", "Not in group.html"}

COLOUR_STATUS = {
    "#00b000": "intragenic_pv",
    "#ff9900": "regulatory_pv",
    "#cccccc": "non_pv_homologue",
}

SUMMARY_HEADERS = {"group", "name", "likelyfunction",
                   "pvingene", "totalpv", "totalgenes"}

COLOUR_RE     = re.compile(r"background-color\s*:\s*(#[0-9a-fA-F]{6})", re.I)
HREF_RE       = re.compile(r"href\s*=\s*[\"']#([^\"']+)[\"']")
PLAIN_LOCI_RE = re.compile(r"\b([A-Z]{2,}_\d{5})\b")
SUBFILE_RE    = re.compile(r"^\d+_\d+\.html$")
GROUP_RE      = re.compile(r"^(\d+)\.html$")


def unescape_text(html: str) -> str:
    html = re.sub(r"(?is)<script.*?>.*?</script>", " ", html)
    html = re.sub(r"(?is)<style.*?>.*?</style>",   " ", html)
    html = re.sub(r"(?is)<.*?>", " ", html)
    return re.sub(r"\s+", " ", unescape(html)).strip()


def norm(t: str) -> str:
    return re.sub(r"[^a-z0-9]+", "", t.lower())


def find_tables(html: str) -> list[str]:
    return re.findall(r"(?is)<table\b.*?>.*?</table>", html)


def parse_table(t: str) -> list[list[str]]:
    rows = []
    for r in re.findall(r"(?is)<tr\b.*?>.*?</tr>", t):
        cells = re.findall(r"(?is)<t[dh]\b.*?>(.*?)</t[dh]>", r)
        rows.append([unescape_text(c) for c in cells])
    return rows


def parse_table_raw(t: str) -> list[list[str]]:
    rows = []
    for r in re.findall(r"(?is)<tr\b.*?>.*?</tr>", t):
        cells = re.findall(r"(?is)<t[dh]\b.*?>(.*?)</t[dh]>", r)
        rows.append(cells)
    return rows


def safe_int(x):
    try: return int(str(x).strip())
    except Exception: return ""


def decode_strain_cell(raw: str) -> dict:
    out = {"pv_status": "absent", "pv_locus": "", "all_loci": ""}
    cm = COLOUR_RE.search(raw)
    if cm:
        c = cm.group(1).lower()
        out["pv_status"] = COLOUR_STATUS.get(c, f"unknown:{c}")
        hm = HREF_RE.search(raw)
        if hm: out["pv_locus"] = hm.group(1)
    all_loci = re.findall(r'href=["\']#([A-Z]+_\d+)["\']', raw)
    plain    = PLAIN_LOCI_RE.findall(unescape_text(raw))
    seen, combined = set(), []
    for l in all_loci + plain:
        if l not in seen:
            seen.add(l); combined.append(l)
    out["all_loci"] = ",".join(combined)
    return out


def parse_summary_table(t, genus, file, page_title):
    raw  = parse_table_raw(t)
    text = parse_table(t)
    if len(text) < 2: return []
    sm = {}
    strain_cols = []
    for i, h in enumerate(text[0]):
        nh = norm(h)
        if   nh == "group":          sm["group_num"] = i
        elif nh == "name":           sm["group_name"] = i
        elif nh == "likelyfunction": sm["likely_function"] = i
        elif nh == "pvingene":       sm["pv_in_gene"] = i
        elif nh == "totalpv":        sm["total_pv"] = i
        elif nh == "totalgenes":     sm["total_genes"] = i
        elif i >= 6:                 strain_cols.append((i, h.strip()))
    if len(sm) < 5: return []
    out_rows = []
    for rr, tr in zip(raw[1:], text[1:]):
        def g(k):
            i = sm.get(k)
            if i is None or i >= len(tr): return ""
            return tr[i]
        pvg = safe_int(g("pv_in_gene")); tp = safe_int(g("total_pv")); tg = safe_int(g("total_genes"))
        reg = ""
        if isinstance(pvg, int) and isinstance(tp, int): reg = tp - pvg
        row = {
            "genus": genus, "file": file, "page_title": page_title,
            "group_num":      g("group_num"),
            "group_name":     g("group_name"),
            "likely_function":g("likely_function"),
            "pv_in_gene":     pvg,
            "total_pv":       tp,
            "total_genes":    tg,
            "regulatory_pv":  reg,
            "has_intragenic_pv": 1 if isinstance(pvg, int) and pvg > 0 else 0,
            "has_any_pv":        1 if isinstance(tp, int) and tp > 0 else 0,
        }
        for ci, lbl in strain_cols:
            raw_cell = rr[ci] if ci < len(rr) else ""
            d = decode_strain_cell(raw_cell)
            sl = re.sub(r"[^a-z0-9]+", "_", lbl.lower()).strip("_")
            row[f"strain_{sl}_status"]   = d["pv_status"]
            row[f"strain_{sl}_pv_locus"] = d["pv_locus"]
            row[f"strain_{sl}_all_loci"] = d["all_loci"]
        out_rows.append(row)
    return out_rows


def parse_annotated_functions(t, genus, file, page_title):
    text = parse_table(t)
    if len(text) < 2: return []
    headers = [norm(h) for h in text[0]]
    def col(ns):
        for n in ns:
            if n in headers: return headers.index(n)
        return None
    fc = col(["function", "likelyfunction", "annotatedfunctions"])
    cc = col(["occuring", "occurring", "count", "n"])
    sc = col(["score"])
    if fc is None: return []
    out = []
    for r in text[1:]:
        if not any(r): continue
        out.append({
            "genus": genus, "file": file, "page_title": page_title,
            "annotated_function": r[fc] if fc is not None and fc < len(r) else "",
            "occuring":           r[cc] if cc is not None and cc < len(r) else "",
            "score":              r[sc] if sc is not None and sc < len(r) else "",
        })
    return out


def is_annot_table(rows):
    if not rows: return False
    h = {norm(x) for x in rows[0]}
    return ("function" in h or "annotatedfunctions" in h) and \
           ("occuring" in h or "occurring" in h or "count" in h)


MEMBER_BLOCK_RE = re.compile(
    r'<p\s+id="([^"]+)">\s*<b>[^:]+:\s*[^<]+</b>\s*'
    r'<a\s+href="([^"]+)">\[tract entry\]</a>.*?'
    r'Function:\s*(.*?)</div>'
    r'.*?white-space\s*:\s*pre[^>]*>(.*?)</div>',
    re.DOTALL,
)


def extract_tract_info(seq_html: str) -> dict:
    spans = re.findall(
        r"<span[^>]*background-color\s*:\s*#0033dd[^>]*>(.*?)</span>",
        seq_html, re.I)
    residues = "".join(spans)
    pos = ""
    for line in seq_html.split("<br />"):
        if "#0033dd" not in line.lower(): continue
        m = re.search(r"\+(\d+)\s+", line)
        if m:
            ofs = int(m.group(1))
            seq_start = re.search(r"\+\d+\s+", line)
            sp = line[seq_start.end():] if seq_start else line
            before = re.split(r"<span", sp)[0]
            before_plain = re.sub(r"<[^>]+>", "", before)
            pos = ofs + len(before_plain) + 1
        break
    return {"tract_residues": residues, "tract_length_aa": len(residues), "tract_aa_pos": pos}


def parse_members(html, genus, file, page_title, group_num):
    out = []
    for m in MEMBER_BLOCK_RE.finditer(html):
        locus, url, fn, sb = m.groups()
        ti = extract_tract_info(sb)
        sm = re.search(r"/strains/(\d+)\.html", url)
        out.append({
            "genus": genus, "file": file, "page_title": page_title,
            "group_num": group_num, "locus_tag": locus,
            "strain_file_idx": sm.group(1) if sm else "",
            "function": re.sub(r"\s+", " ", fn).strip(),
            **ti,
        })
    return out


PAIRWISE_RE = re.compile(
    r'(\w+)\s+vs:\s+(\w+)\s+in\s+<a\s+href="([^"]+)">(.*?)</a>[^<]*<br/>'
    r'\s*<span[^>]*>Gene length:\s*(\d+)bp\s*/\s*(\d+)aa\s+PV:\s*(Yes|No)</span>'
    r'\s*<br\s*/>\s*<span[^>]*>Function:\s*(.*?)</span>'
    r'.*?Score:\s*([\d.]+)\s+bits:\s*([\d.]+)\s+e-value:\s*(\S+)<br\s*/>'
    r'\s*length:\s*(\d+)\s+gaps:\s*(\d+)\s+id:\s*(\d+)\s+'
    r'positives:\s*(\d+)\s+coverage:\s*([\d.]+)\s+query coverage\s+([\d.]+)',
    re.DOTALL,
)


def parse_pairwise(html, genus, file, page_title, group_num):
    out = []
    for m in PAIRWISE_RE.finditer(html):
        (q, t, _u, sn, bp, aa, pv, fn, sc, b, e, l, gp, idn, po, cov, qcov) = m.groups()
        out.append({
            "genus": genus, "file": file, "page_title": page_title,
            "group_num": group_num,
            "query_locus": q, "target_locus": t,
            "strain_name": sn.strip(),
            "gene_length_bp": safe_int(bp), "gene_length_aa": safe_int(aa),
            "pv_status": pv, "function": fn.strip(),
            "blast_score": sc, "blast_bits": b, "blast_evalue": e,
            "align_length": safe_int(l), "align_gaps": safe_int(gp),
            "align_identity": safe_int(idn), "align_positives": safe_int(po),
            "align_coverage": cov, "align_query_cov": qcov,
        })
    return out


def write_csv(path: Path, fields: list[str], rows: list[dict]) -> None:
    if not rows:
        with path.open("w", newline="") as fh:
            csv.DictWriter(fh, fieldnames=fields, extrasaction="ignore").writeheader()
        return
    keys = list(dict.fromkeys(fields + [k for r in rows for k in r if k not in fields]))
    with path.open("w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=keys, extrasaction="ignore", restval="")
        w.writeheader(); w.writerows(rows)


def main() -> None:
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--base",   required=True, type=Path)
    p.add_argument("--outdir", required=True, type=Path)
    a = p.parse_args()
    a.outdir.mkdir(parents=True, exist_ok=True)

    groups, annots, members, pairs = [], [], [], []
    for gd in sorted(x for x in a.base.iterdir() if x.is_dir()):
        gdir = gd / "summary_tracts" / "groups"
        if not gdir.is_dir(): continue
        for hf in sorted(gdir.glob("*.html")):
            if hf.name in SKIP_FILES: continue
            try: html = hf.read_text(errors="ignore")
            except Exception: continue
            genus = gd.name
            gn = (re.match(r"^(\d+)", hf.name).group(1)
                  if re.match(r"^(\d+)", hf.name) else "")
            tm = re.search(r"(?is)<title>(.*?)</title>", html)
            pt = unescape_text(tm.group(1)) if tm else ""
            if SUBFILE_RE.match(hf.name):
                pairs.extend(parse_pairwise(html, genus, hf.name, pt, gn))
                continue
            if not GROUP_RE.match(hf.name): continue
            for tbl in find_tables(html):
                rows = parse_table(tbl)
                if not rows: continue
                hs = {norm(h) for h in rows[0]}
                if len(hs & SUMMARY_HEADERS) >= 5:
                    groups.extend(parse_summary_table(tbl, genus, hf.name, pt))
                elif is_annot_table(rows):
                    annots.extend(parse_annotated_functions(tbl, genus, hf.name, pt))
            members.extend(parse_members(html, genus, hf.name, pt, gn))

    write_csv(a.outdir / "phasomeit_group_summary.csv",
              ["genus", "file", "page_title", "group_num", "group_name",
               "likely_function", "pv_in_gene", "total_pv", "total_genes",
               "regulatory_pv", "has_intragenic_pv", "has_any_pv"], groups)
    write_csv(a.outdir / "phasomeit_annotated_functions.csv",
              ["genus", "file", "page_title", "annotated_function",
               "occuring", "score"], annots)
    write_csv(a.outdir / "phasomeit_members.csv",
              ["genus", "file", "page_title", "group_num", "locus_tag",
               "strain_file_idx", "function", "tract_residues",
               "tract_length_aa", "tract_aa_pos"], members)
    write_csv(a.outdir / "phasomeit_pairwise.csv",
              ["genus", "file", "page_title", "group_num", "query_locus",
               "target_locus", "strain_name", "gene_length_bp",
               "gene_length_aa", "pv_status", "function", "blast_score",
               "blast_bits", "blast_evalue", "align_length", "align_gaps",
               "align_identity", "align_positives", "align_coverage",
               "align_query_cov"], pairs)
    print(f"out: {a.outdir}")
    print(f"  groups   {len(groups):,}")
    print(f"  annots   {len(annots):,}")
    print(f"  members  {len(members):,}")
    print(f"  pairs    {len(pairs):,}")


if __name__ == "__main__":
    main()

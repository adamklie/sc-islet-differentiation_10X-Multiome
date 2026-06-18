#!/usr/bin/env python
"""Searchable/sortable index across all per-cluster constituent sheets.

Chimerism score (Spearman-based, consistent with the sheets):
  - n_distinct_tf : # distinct top-TomTom TF symbols among constituents
  - profile_min_rho : min pairwise SPEARMAN between constituent marginalization
                      profiles across 22 cell types (low = constituents diverge)
Flag CHIMERIC if n_distinct_tf>=2 AND profile_min_rho<0.5.
"""
from __future__ import annotations
import json, csv, html
from pathlib import Path
from itertools import combinations
import numpy as np
import pandas as pd
from scipy.stats import spearmanr

CONST = Path("/cellar/users/aklie/data/datasets/sc-islet-differentiation_10X-Multiome/results/3_single_task_models/motifs/gimme_cluster/marginalization/constituent")
TT = CONST/"tomtom"
manifest = json.load(open(CONST/"manifest.json"))
marg = pd.read_csv(CONST/"marginalization_matrix.csv", index_col=0)


def top_tf(name):
    p = TT/f"{name}.tomtom.txt"
    if not p.is_file():
        return None
    best = None
    with open(p) as f:
        for row in csv.DictReader(f, delimiter="\t"):
            tid, ev = row.get("Target_ID"), row.get("E-value")
            if not tid or ev in (None, ""):
                continue
            try: ev = float(ev)
            except ValueError: continue
            if best is None or ev < best[1]:
                best = (tid.split("_")[0], ev)
    return best[0] if best else None


rows = []
for cl, info in manifest.items():
    sid = info["short_id"]
    names = [c["name"] for c in info["constituents"]]
    tfs = [t for t in (top_tf(n) for n in names) if t]
    distinct = sorted(set(tfs))
    present = [n for n in names if n in marg.index]
    min_rho = np.nan
    if len(present) >= 2:
        rs = []
        for a, b in combinations(present, 2):
            va, vb = marg.loc[a].astype(float).values, marg.loc[b].astype(float).values
            if va.std() > 0 and vb.std() > 0:
                rs.append(spearmanr(va, vb)[0])
        if rs:
            min_rho = float(np.nanmin(rs))
    chimeric = (len(distinct) >= 2) and (not np.isnan(min_rho)) and (min_rho < 0.5)
    rows.append(dict(short_id=sid, cluster=cl, consensus_tf=info["tf_name"],
                     n_constituents=len(names), n_distinct_tf=len(distinct),
                     distinct_tfs=distinct, profile_min_rho=min_rho, chimeric=chimeric,
                     sheet=f"sheets/{sid}_{cl}.html"))

df = pd.DataFrame(rows)
df["_st"] = df["short_id"].str.extract(r"ST(\d+)").astype(float)
df = df.sort_values(["chimeric", "n_distinct_tf", "_st"], ascending=[False, False, True]).drop(columns="_st")
df.assign(distinct_tfs=df["distinct_tfs"].apply(lambda x: ";".join(x))).to_csv(CONST/"chimerism_summary.tsv", sep="\t", index=False)

def chips(tfs): return " ".join(f'<span class="chip">{html.escape(t)}</span>' for t in tfs)
trs = []
for _, r in df.iterrows():
    flag = '<span class="flag yes">chimeric</span>' if r["chimeric"] else '<span class="flag no">—</span>'
    mr = "n/a" if pd.isna(r["profile_min_rho"]) else f'{r["profile_min_rho"]:.2f}'
    trs.append(f'<tr data-chimeric="{int(r["chimeric"])}"><td><a href="{r["sheet"]}">{html.escape(r["short_id"])}</a></td>'
               f'<td class="mono">{html.escape(r["cluster"])}</td><td><b>{html.escape(str(r["consensus_tf"]))}</b></td>'
               f'<td class="num">{r["n_constituents"]}</td><td class="num">{r["n_distinct_tf"]}</td>'
               f'<td class="num">{mr}</td><td>{flag}</td><td>{chips(r["distinct_tfs"])}</td></tr>')
n_chim = int(df["chimeric"].sum())
HTML = f"""<!doctype html><html><head><meta charset="utf-8"><title>Constituent chimerism index</title>
<style>
body{{font-family:-apple-system,Arial,sans-serif;max-width:1200px;margin:18px auto;color:#222}}
h1{{border-bottom:2px solid #333;padding-bottom:6px}}
.meta{{background:#fffae0;border-left:3px solid #f0c040;padding:8px 12px;font-size:13px;margin:10px 0}}
table{{border-collapse:collapse;width:100%;font-size:13px}}
th,td{{border-bottom:1px solid #e3e3e3;padding:6px 9px;text-align:left}}
th{{background:#f4f4f4;cursor:pointer;position:sticky;top:0}}th:hover{{background:#e8e8e8}}
td.num,td.mono{{font-family:ui-monospace,monospace}}td.num{{text-align:right}}
tr[data-chimeric="1"]{{background:#fff3f3}}
.chip{{display:inline-block;background:#e8f0fe;padding:1px 7px;margin:1px;border-radius:4px;font-family:monospace;font-size:11px}}
.flag.yes{{background:#b2182b;color:#fff;padding:1px 8px;border-radius:8px;font-size:11px;font-weight:bold}}.flag.no{{color:#bbb}}
.controls{{margin:10px 0}}input{{padding:5px;font-size:14px}}a{{color:#1e6091;text-decoration:none}}a:hover{{text-decoration:underline}}
</style></head><body>
<h1>ChromBPNet constituent-motif chimerism index</h1>
<div class="meta">{len(df)} multi-constituent clusters · <b>{n_chim} flagged chimeric</b>
(≥2 distinct constituent TomTom TFs AND min pairwise marginalization-profile Spearman ρ &lt; 0.5).
Each ST links to its constituent sheet. Click a header to sort; use the box to filter.</div>
<div class="controls"><input id="q" placeholder="filter: ST, cluster, TF…" oninput="flt()" size="40">
<label><input type="checkbox" id="chim" onchange="flt()"> chimeric only</label></div>
<table id="t"><thead><tr>
<th onclick="srt(0)">ST</th><th onclick="srt(1)">cluster</th><th onclick="srt(2)">consensus TF</th>
<th onclick="srt(3)">#constit</th><th onclick="srt(4)">#distinct TF</th><th onclick="srt(5)">min ρ</th>
<th onclick="srt(6)">flag</th><th>constituent TFs</th></tr></thead><tbody>
{''.join(trs)}
</tbody></table>
<script>
function flt(){{const q=document.getElementById('q').value.toLowerCase(),co=document.getElementById('chim').checked;
document.querySelectorAll('#t tbody tr').forEach(tr=>{{tr.style.display=((!q||tr.innerText.toLowerCase().includes(q))&&(!co||tr.dataset.chimeric==='1'))?'':'none';}});}}
function srt(c){{const tb=document.querySelector('#t tbody'),rows=[...tb.rows],num=(c>=3&&c<=5),
asc=tb.dataset.sc!=String(c)||tb.dataset.asc!=='1';
rows.sort((a,b)=>{{let x=a.cells[c].innerText.trim(),y=b.cells[c].innerText.trim();
if(num){{x=parseFloat(x)||-1e9;y=parseFloat(y)||-1e9;return asc?x-y:y-x;}}return asc?x.localeCompare(y):y.localeCompare(x);}});
rows.forEach(r=>tb.appendChild(r));tb.dataset.sc=c;tb.dataset.asc=asc?'1':'0';}}
</script></body></html>"""
(CONST/"index.html").write_text(HTML)
print(f"Wrote {CONST/'index.html'} ({len(df)} clusters, {n_chim} chimeric)")
print(df.drop(columns=["sheet","distinct_tfs"]).head(20).to_string(index=False))

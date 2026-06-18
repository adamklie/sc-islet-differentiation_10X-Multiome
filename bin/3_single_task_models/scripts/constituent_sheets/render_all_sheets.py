#!/usr/bin/env python
"""Render ONE multi-constituent cluster's constituent sheet (approved ST36 layout).
Run per-cluster via SLURM array. Reuses bundle primitives (logos, TomTom logos,
expression heatmap) + custom assembly: clustered summary heatmap with row logos,
marg panel, TomTom table (with query RC), full-archetype candidates (expression-
prefiltered) with clustered expression heatmap, Spearman-ordered scatters.
"""
from __future__ import annotations
import sys, csv, io, base64, html as _html, argparse, json
from pathlib import Path
import numpy as np
import pandas as pd
from scipy.stats import pearsonr, spearmanr
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.offsetbox import OffsetImage, AnnotationBbox
import matplotlib.image as mpimg
import seaborn as sns
from scipy.cluster.hierarchy import linkage, leaves_list

REVIEW = Path("/cellar/users/aklie/projects/islet_organoid_differentiation/sc-islet-differentiation_10X-Multiome/scratch/2026_05_25_motif_review")
sys.path.insert(0, str(REVIEW))
from lib.annotation_report import (
    _load_meme_pwms, _parse_tomtom, _aligned_rows, _pwm_logo_b64,
    _render_db_target_logo, _set_current_db_pwms, _color_seq, DEFAULT_VIERSTRA_MEME)
from lib.tf_block import _heatmap_png, load_context, _resolve_candidates

# Self-contained constituent-sheet workspace under the `constituent` summary dir:
# query.meme + manifest.json + marginalization_matrix*.csv are written there by the
# marginalization pipeline; tomtom/ + sheets/ are produced by the sheet scripts.
CONST = Path("/cellar/users/aklie/data/datasets/sc-islet-differentiation_10X-Multiome/results/3_single_task_models/motifs/gimme_cluster/marginalization/constituent")
INP = ANA = CONST          # query.meme, manifest.json, marginalization_matrix*.csv
TT = CONST/"tomtom"
SHEETS = CONST/"sheets"
META_TSV = "/cellar/users/aklie/projects/islet_organoid_differentiation/sc-islet-differentiation_10X-Multiome/config/cell_type_metadata.tsv"
CPM_CSV = "/cellar/users/aklie/data/datasets/sc-islet-differentiation_10X-Multiome/results/2_process_data/pseudobulk_all_vierstra_tf_cpm.csv"
# Legacy islet-development panel — now only a LAST-RESORT fallback for the summary
# heatmap x-axis. The x-axis is built per-cluster from each cluster's own candidate TFs
# (see build_cluster_sheet); FOCUS is used only if a cluster has no expressed candidates.
FOCUS = ["PDX1","ISL1","NKX6-1","NKX2-2","NEUROD1","LMX1A","LMX1B","NOTO","ALX3",
         "HOXB4","HOXA2","BARX1","MIXL1","HESX1","ARX","MAFB","PAX4","PAX6","FOXA2"]
EVAL_TABLE = 1.0
MIN_EXPR_CPM = 5.0
SCATTER_K = 24
SUMMARY_TF_CAP = 24   # max TFs on the summary-heatmap x-axis (ranked by max |rho| if exceeded)


def tier(e):
    e = float(e)
    return ("excellent" if e < 1e-6 else "strong" if e < 1e-3 else
            "moderate" if e < 5e-2 else "weak" if e < 5e-1 else "poor")


def fig_b64(fig, dpi=110):
    buf = io.BytesIO(); fig.savefig(buf, format="png", dpi=dpi, bbox_inches="tight"); plt.close(fig)
    return "data:image/png;base64,"+base64.b64encode(buf.getvalue()).decode()


def b64_to_arr(b64):
    return mpimg.imread(io.BytesIO(base64.b64decode(b64)))


def cluster_rows(df):
    if df.shape[0] < 3:
        return df
    X = np.log1p(df.values.astype(float)); keep = np.where(X.std(axis=1) > 0)[0]
    if len(keep) < 3:
        return df
    try:
        Z = linkage(X[keep], method="average", metric="euclidean")
        order = list(keep[leaves_list(Z)]) + [i for i in range(df.shape[0]) if i not in keep]
        return df.iloc[order]
    except Exception:
        return df


def marg_panel(raw, z, ct_order, ct_colors):
    fig, ax = plt.subplots(2, 1, figsize=(11, 3.2), dpi=95, sharex=True, gridspec_kw={"hspace": 0.35})
    x = np.arange(len(ct_order)); cols = [ct_colors.get(c, "#888") for c in ct_order]
    ax[0].bar(x, raw, color=cols, width=0.85); ax[0].axhline(0, color="#888", lw=0.4)
    ax[0].set_title("marg raw", fontsize=9, loc="left", color="#444")
    m = float(np.max(np.abs(z))) or 1.0
    ax[1].bar(x, z, color=cols, width=0.85); ax[1].axhline(0, color="#888", lw=0.4)
    ax[1].set_ylim(-m*1.05, m*1.05); ax[1].set_title("marg z-score", fontsize=9, loc="left", color="#444")
    ax[1].set_xticks(x); ax[1].set_xticklabels(ct_order, rotation=90, fontsize=7)
    for a in ax: a.spines["top"].set_visible(False); a.spines["right"].set_visible(False); a.tick_params(axis="y", labelsize=7)
    return fig_b64(fig, dpi=95)


def scatter_grid(marg_row, rows, cpm, cts, ct_colors):
    n = len(rows)
    if n == 0:
        return ""
    ncols = min(4, n); nrows = int(np.ceil(n / ncols))
    if n == 1:
        panel_w, legend_w = 4.2, 2.8; fig_w = panel_w + legend_w; legend_ncol = 2; rect_right = panel_w / fig_w
    else:
        fig_w = 3.7 * ncols; legend_ncol = 1; rect_right = 0.86
    fig, axes = plt.subplots(nrows, ncols, figsize=(fig_w, 3.0 * nrows), squeeze=False)
    colors = [ct_colors.get(ct, "#888888") for ct in cts]
    y = marg_row.loc[cts].astype(float).values
    for i, (tf, pr, pp, sr, sp) in enumerate(rows):
        ax = axes[i // ncols][i % ncols]
        x = np.log1p(cpm[tf].loc[cts].astype(float).values)
        ax.scatter(x, y, s=52, c=colors, edgecolors="black", linewidths=0.4, alpha=0.95)
        ax.set_title(f"{tf}   ρ={sr:+.2f}   r={pr:+.2f}", fontsize=12)
        ax.axhline(0, color="grey", lw=0.5)
        ax.set_xlabel("log(CPM+1)", fontsize=10.5); ax.set_ylabel("marg Δ", fontsize=10.5)
        ax.tick_params(labelsize=9)
    for j in range(n, nrows * ncols):
        axes[j // ncols][j % ncols].set_axis_off()
    handles = [Line2D([0], [0], marker="o", color="w", markerfacecolor=ct_colors.get(ct, "#888"),
                      markeredgecolor="black", markeredgewidth=0.4, markersize=8, label=ct) for ct in cts]
    fig.legend(handles=handles, loc="center right", bbox_to_anchor=(1.0, 0.5), fontsize=9,
               frameon=False, ncol=legend_ncol, title="cell type", title_fontsize=10)
    fig.tight_layout(rect=(0, 0, rect_right, 1))
    return fig_b64(fig)


def tomtom_table(name, db_pwms, fwd_b64, rev_b64, evalue_threshold):
    hits = _parse_tomtom(TT/f"{name}.tomtom.txt", evalue_threshold=evalue_threshold, min_rows=1, max_rows=8)
    if not hits:
        return '<p class="no-hits">No TomTom hits.</p>'
    common_len = max(len(hits[0]["query_consensus"]),
                     max(max(0, h["offset"]) + len(h["target_consensus"]) for h in hits))
    qlogo = f'<img class="seq-logo" src="data:image/png;base64,{fwd_b64}">' if fwd_b64 else f'<span class="seq">{_color_seq(hits[0]["query_consensus"])}</span>'
    qrc = f'<img class="seq-logo" src="data:image/png;base64,{rev_b64}">' if rev_b64 else "—"
    rows = [f'<tr class="query-row"><td>—</td><td><b>this motif</b></td><td>—</td><td>—</td>'
            f'<td>{qlogo}</td><td>{qrc}</td><td><span class="seq">{_color_seq(hits[0]["query_consensus"])}</span></td></tr>']
    for j, h in enumerate(hits):
        _, t_row = _aligned_rows(h["query_consensus"], h["target_consensus"], h["offset"], h["orientation"])
        if h["target"] in db_pwms:
            fwd = _render_db_target_logo(h["target"], "+", 0, None); rc = _render_db_target_logo(h["target"], "-", 0, None)
            al = _render_db_target_logo(h["target"], h["orientation"], max(0, h["offset"]), common_len)
            fwd_h = f'<img class="seq-logo" src="data:image/png;base64,{fwd}">' if fwd else "—"
            rc_h = f'<img class="seq-logo" src="data:image/png;base64,{rc}">' if rc else "—"
            al_h = f'<img class="seq-logo" src="data:image/png;base64,{al}"><br><span class="seq">{_color_seq(t_row)}</span>' if al else f'<span class="seq">{_color_seq(t_row)}</span>'
        else:
            fwd_h = rc_h = "—"; al_h = f'<span class="seq">{_color_seq(t_row)}</span>'
        rows.append(f'<tr><td>{j+1}</td><td>{_html.escape(h["target"])}</td><td>{h["evalue"]:.2e}</td>'
                    f'<td>{h["qvalue"]:.2e}</td><td>{fwd_h}</td><td>{rc_h}</td>'
                    f'<td>{al_h}<br><small>strand {h["orientation"]}, offset {h["offset"]}</small></td></tr>')
    return ('<table class="hits-tbl"><thead><tr><th>rank</th><th>target</th><th>e-value</th><th>q-value</th>'
            '<th>forward PWM</th><th>reverse complement</th><th>alignment to query</th></tr></thead>'
            f'<tbody>{"".join(rows)}</tbody></table>')


def cluster_card(name, ctx):
    hits = _parse_tomtom(TT/f"{name}.tomtom.txt", evalue_threshold=1.0, min_rows=1, max_rows=20)
    if not hits:
        return ""
    cands, cid, info, extra = _resolve_candidates(hits, ctx)
    top = hits[0]["target"]
    if cid is None:
        return f'<p class="tf-note">Top hit <code>{_html.escape(top)}</code> not in Vierstra archetype DB.</p>'
    return (f'<table class="tf-cluster-card">'
            f'<tr><td><b>Top TomTom hit</b></td><td><code>{_html.escape(top)}</code></td></tr>'
            f'<tr><td><b>Primary Vierstra cluster</b></td><td>{cid} · {info.get("Name","?")} · '
            f'DBD: {info.get("DBD","?")} · {info.get("Cluster_size","?")} member motifs</td></tr></table>')


def build_cluster_sheet(cl, info, shared):
    marg, margz, cpm, expressed, ct_order, ct_colors, ctx, db_pwms, pwms = shared
    name_by_sid = {}
    rows_def = []
    sid = info["short_id"]
    if info.get("consensus_name"):
        rows_def.append((f"{sid}_consensus", info["consensus_name"], "cluster consensus", info["n_patterns"]))
    for k, c in enumerate(info["constituents"], 1):
        rows_def.append((f"{sid}_c{k}_{c['celltype']}", c["name"], c["celltype"], 1))
    for r in rows_def:
        name_by_sid[r[0]] = r[1]

    def corr(motif, tf):
        y = marg.loc[motif, ct_order].astype(float).values
        x = np.log1p(cpm[tf].astype(float).values)
        if np.std(x) == 0 or np.std(y) == 0:
            return np.nan, np.nan, np.nan, np.nan
        pr, pp = pearsonr(x, y); sr, sp = spearmanr(x, y)
        return pr, pp, sr, sp

    # summary clustered heatmap (motifs x cluster-relevant TFs) with row logos.
    # x-axis is built PER CLUSTER from the union of each constituent's expressed Vierstra
    # archetype candidate TFs — so each cluster shows its own relevant TFs rather than a
    # global islet panel. Falls back to data-driven top-|rho| TFs, then to FOCUS.
    def _cand_tfs(motif_name):
        h = _parse_tomtom(TT/f"{motif_name}.tomtom.txt", evalue_threshold=1.0, min_rows=1, max_rows=20)
        if not h:
            return []
        cands, _cid, _info, _extra = _resolve_candidates(h, ctx)
        return [t for t in cands if t in expressed]

    focus, best_abs = [], {}
    for r in rows_def:
        for tf in _cand_tfs(r[1]):
            rho = corr(r[1], tf)[2]
            best_abs[tf] = max(best_abs.get(tf, 0.0), abs(rho) if rho == rho else 0.0)
            if tf not in focus:
                focus.append(tf)
    if len(focus) > SUMMARY_TF_CAP:
        focus = sorted(focus, key=lambda t: best_abs.get(t, 0.0), reverse=True)[:SUMMARY_TF_CAP]
    if not focus:   # fallback 1: data-driven top-|rho| expressed TFs across constituents
        scores = {}
        for r in rows_def:
            for tf in (t for t in cpm.columns if t in expressed):
                rho = corr(r[1], tf)[2]
                if rho == rho:
                    scores[tf] = max(scores.get(tf, 0.0), abs(rho))
        focus = [t for t, _ in sorted(scores.items(), key=lambda kv: kv[1], reverse=True)[:SUMMARY_TF_CAP]]
    if not focus:   # fallback 2: legacy islet panel
        focus = [t for t in FOCUS if t in cpm.columns]
    print(f"[{info['short_id']} {cl}] summary-heatmap x-axis ({len(focus)} TFs): {focus}", file=sys.stderr)
    H = pd.DataFrame(index=[r[0] for r in rows_def], columns=focus, dtype=float)
    for sid_ in H.index:
        for tf in focus:
            H.loc[sid_, tf] = corr(name_by_sid[sid_], tf)[2]
    Hc = H.fillna(0.0)
    g = sns.clustermap(Hc, cmap="RdBu_r", vmin=-1, vmax=1, center=0, annot=True, fmt="+.2f",
                       annot_kws={"size": 9.5}, figsize=(0.62*len(focus)+4.0, 0.66*len(H)+3.2),
                       dendrogram_ratio=(0.12, 0.18), linewidths=0.4, linecolor="#eee",
                       cbar_pos=(0.02, 0.83, 0.025, 0.13), cbar_kws={"label": "Spearman ρ"})
    ax = g.ax_heatmap; ax.set_xlabel(""); ax.set_ylabel(""); ax.tick_params(labelsize=10)
    plt.setp(ax.get_xticklabels(), rotation=90); ax.set_yticks([])
    trans = ax.get_yaxis_transform()
    for kk, orig in enumerate(g.dendrogram_row.reordered_ind):
        sid_ = H.index[orig]
        img = b64_to_arr(_pwm_logo_b64(pwms[name_by_sid[sid_]], "+"))
        ax.add_artist(AnnotationBbox(OffsetImage(img, zoom=0.85), (1.012, kk+0.5), xycoords=trans,
                                     box_alignment=(0, 0.5), frameon=False, annotation_clip=False))
        ax.annotate(sid_, (1.20, kk+0.5), xycoords=trans, va="center", ha="left", fontsize=9, annotation_clip=False)
    g.fig.suptitle(f"{sid} ({cl}, {info['tf_name']}) — constituent marg vs TF expression, clustered (Spearman ρ)",
                   fontsize=12, y=1.02)
    heat = fig_b64(g.fig)

    sids = [r[0] for r in rows_def]
    sections = []
    for idx, (s_id, name, src, npat) in enumerate(rows_def):
        fwd_b64 = _pwm_logo_b64(pwms[name], "+"); rev_b64 = _pwm_logo_b64(pwms[name], "-")
        th = _parse_tomtom(TT/f"{name}.tomtom.txt", evalue_threshold=1.0, min_rows=1, max_rows=1)
        top_tf = th[0]["target"].split("_")[0] if th else "?"
        qual = tier(th[0]["evalue"]) if th else "poor"
        panel = marg_panel(marg.loc[name, ct_order].astype(float).values,
                           margz.loc[name, ct_order].astype(float).values, ct_order, ct_colors)
        tbl = tomtom_table(name, db_pwms, fwd_b64, rev_b64, EVAL_TABLE)
        hits_full = _parse_tomtom(TT/f"{name}.tomtom.txt", evalue_threshold=1.0, min_rows=1, max_rows=20)
        cands, cid, info_c, extra = _resolve_candidates(hits_full, ctx) if hits_full else ([], None, {}, [])
        expr_avail = [t for t in cands if t in expressed]; expr_missing = [t for t in cands if t not in expressed]
        crows = [(tf,)+corr(name, tf) for tf in expr_avail]
        crows_sp = sorted([c for c in crows if not np.isnan(c[3])], key=lambda c: c[3], reverse=True)
        allc = [(tf,)+corr(name, tf) for tf in cpm.columns if tf in expressed]
        allc = [c for c in allc if not np.isnan(c[3])]
        top_gw = sorted(allc, key=lambda c: c[3], reverse=True)[:8]
        card = cluster_card(name, ctx)
        heatmap_html = scatter_html = ""
        if expr_avail:
            expr_sub = cluster_rows(cpm[expr_avail].T)
            hm = _heatmap_png(expr_sub, ct_colors)
            heatmap_html = (f'<div class="lbl">TF expression across cell types — Vierstra archetype candidates '
                            f'({len(expr_avail)} expressed TFs), TF rows clustered</div>'
                            f'<img class="wide" src="data:image/png;base64,{hm}">')
        if crows_sp:
            sel = sorted(crows_sp, key=lambda c: abs(c[3]), reverse=True)[:SCATTER_K]
            sel = sorted(sel, key=lambda c: c[3], reverse=True)
            note = f" — top {SCATTER_K} of {len(crows_sp)} by |ρ|" if len(crows_sp) > SCATTER_K else ""
            sc = scatter_grid(marg.loc[name, ct_order], sel, cpm, ct_order, ct_colors)
            scatter_html = (f'<div class="lbl">Marginalization Δ vs expression — archetype candidates ordered by '
                            f'Spearman ρ (high→low){note}; each panel shows ρ and Pearson r</div><img class="wide" src="{sc}">')
        chips = " ".join(f'<span class="tf-chip {"missing" if t in expr_missing else ""}">{t}{" (no expr)" if t in expr_missing else ""}</span>'
                         for t in cands) or '<span class="dim">none</span>'

        def corr_tbl(rs):
            tr = "".join(f"<tr><td><b>{t}</b></td><td>{sr:+.3f}</td><td>{sp:.1e}</td><td>{pr:+.3f}</td></tr>"
                         for (t, pr, pp, sr, sp) in rs)
            return ('<table class="corr-tbl"><thead><tr><th>TF</th><th>Spearman ρ</th><th>p</th><th>Pearson r</th></tr></thead>'
                    f'<tbody>{tr}</tbody></table>')
        nav = ((f'<a href="#{sids[idx-1]}">← prev</a>' if idx > 0 else '<span class="dim">← prev</span>')
               + ' · <a href="#top">top</a> · '
               + (f'<a href="#{sids[idx+1]}">next →</a>' if idx < len(rows_def)-1 else '<span class="dim">next →</span>'))
        sections.append(f"""
<section id="{s_id}" class="page">
  <div class="page-header"><div class="page-title">
    <span class="sid">{s_id}</span><span class="tf">{top_tf}</span>
    <span class="quality q-{qual}">{qual}</span><span class="src">{_html.escape(src)}</span></div>
    <div class="page-meta">page {idx+1}/{len(rows_def)} &nbsp; {nav}</div></div>
  <div class="logos">
    <div class="logo-block"><div class="logo-label">forward</div><img src="data:image/png;base64,{fwd_b64}"></div>
    <div class="logo-block"><div class="logo-label">reverse complement</div><img src="data:image/png;base64,{rev_b64}"></div></div>
  <h3>Marginalization</h3><div class="dist-panel"><img src="{panel}"></div>
  <h3>TomTom hits (Vierstra v2.0beta)</h3>{tbl}
  <h3>Candidate TFs (Vierstra archetype, expression-prefiltered)</h3>
  {card}<div class="chips">{chips}</div>
  <div class="two-col">
    <div><div class="lbl">Candidate-TF correlations (Spearman↓, top 15)</div>{corr_tbl(crows_sp[:15])}</div>
    <div><div class="lbl">Top expression-correlated TFs (database-wide, Spearman↓)</div>{corr_tbl(top_gw)}</div></div>
  {heatmap_html}{scatter_html}
</section>""")

    css = CSS
    out = SHEETS/f"{sid}_{cl}.html"
    out.write_text(f"""<!doctype html><html><head><meta charset="utf-8"><title>{sid} {cl} constituents</title><style>{css}</style></head><body>
<a id="top"></a><h1>{sid} ({cl}, {info['tf_name']}) — constituent-motif report</h1>
<div class="meta">Is this cluster chimeric? Per constituent: marginalization across 22 cell types, TomTom hits, Vierstra archetype candidates (expressed ≥{MIN_EXPR_CPM:g} CPM), clustered expression heatmap, and marg-vs-expression scatters ordered by Spearman ρ.</div>
<h2>Summary — constituent × TF correlation (clustered)</h2><img class="wide" src="{heat}">
<div class="sticky"><input id="search" class="search-box" oninput="flt()" placeholder="Jump to motif: type ID, TF, cell type…"></div>
{''.join(sections)}
<script>function flt(){{const q=document.getElementById('search').value.toLowerCase();document.querySelectorAll('.page').forEach(p=>{{p.style.display=(!q||p.innerText.toLowerCase().includes(q))?'':'none';}});}}</script>
</body></html>""")
    return out


CSS = """
body{font-family:-apple-system,Arial,sans-serif;max-width:1180px;margin:0 auto;padding:16px 20px;color:#222;line-height:1.45}
h1{border-bottom:2px solid #333;padding-bottom:6px}h2{margin-top:8px}
h3{margin-top:16px;margin-bottom:6px;color:#444;border-bottom:1px solid #eee;padding-bottom:2px}
.meta{background:#fffae0;border-left:3px solid #f0c040;padding:8px 12px;font-size:13px;margin:10px 0}
.page{page-break-after:always;padding:18px 20px;border:1px solid #ccc;border-radius:6px;margin-bottom:22px;background:#fafafa}
.page-header{display:flex;flex-wrap:wrap;align-items:baseline;justify-content:space-between;gap:8px;padding-bottom:8px;border-bottom:2px solid #333}
.page-title{display:flex;gap:12px;flex-wrap:wrap;align-items:baseline}
.sid{font-family:ui-monospace,monospace;font-weight:bold;font-size:20px;color:#0a4068;background:#e6f0f8;padding:4px 12px;border-radius:4px}
.tf{font-size:20px;font-weight:bold}.src{color:#777;font-size:13px}
.quality{font-size:12px;padding:2px 10px;border-radius:10px;color:#fff;font-weight:bold}
.q-excellent{background:#1b7837}.q-strong{background:#5aae61}.q-moderate{background:#d9d35e;color:#333}.q-weak{background:#ef8a62}.q-poor{background:#b2182b}
.page-meta{font-size:12px;color:#555}.page-meta a{color:#1e6091;text-decoration:none}.dim{color:#999;font-size:12px}
.logos{display:flex;gap:18px;margin-top:10px}.logo-block{flex:1;text-align:center}.logo-label{font-size:11px;color:#777;margin-bottom:4px}
.logo-block img{max-width:100%;height:auto;max-height:90px;border:1px solid #ddd;background:#fff;padding:4px;border-radius:3px}
.dist-panel img{max-width:100%}.wide{max-width:100%}
.hits-tbl{border-collapse:collapse;width:100%;font-size:12px;margin-top:6px}
.hits-tbl th,.hits-tbl td{border-bottom:1px solid #eee;padding:4px 8px}.hits-tbl th{background:#f4f4f4;text-align:left}
.query-row td{background:#f0f7ff;font-weight:bold}.seq-logo{height:28px;vertical-align:middle;image-rendering:pixelated}
.seq{font-family:ui-monospace,monospace;font-weight:bold;letter-spacing:1px;font-size:12px}.no-hits{color:#b22;font-style:italic}
.chips{margin:6px 0}.tf-chip{display:inline-block;background:#e8f0fe;padding:2px 8px;margin:2px;border-radius:4px;font-family:monospace;font-size:11px}
.tf-chip.missing{background:#f4f4f4;color:#888;font-style:italic}
.corr-tbl{border-collapse:collapse;font-size:12px;margin-top:4px}.corr-tbl th{background:#eee;padding:4px 8px;text-align:left}.corr-tbl td{padding:3px 8px;border-bottom:1px solid #eee}
.tf-cluster-card{background:#f0f6ff;padding:8px 12px;border-radius:5px;font-size:13px;border-collapse:collapse;width:100%}.tf-cluster-card td{padding:3px 8px}
.two-col{display:flex;gap:24px;flex-wrap:wrap;margin-top:8px}.lbl{font-weight:600;font-size:13px;color:#444;margin:12px 0 4px}
.search-box{width:100%;padding:6px;margin:8px 0;font-size:14px}.sticky{position:sticky;top:0;background:#fff;padding:6px 0;border-bottom:1px solid #ccc;z-index:10;margin-bottom:14px}
.tf-note{color:#777;font-size:12px}
@media print{.sticky,.page-meta,.meta{display:none}.page{border:none;background:#fff}}
"""


def main():
    ap = argparse.ArgumentParser(); ap.add_argument("--cluster"); ap.add_argument("--index", type=int)
    args = ap.parse_args()
    SHEETS.mkdir(parents=True, exist_ok=True)
    manifest = json.load(open(CONST/"manifest.json"))
    clusters = list(manifest.keys())
    if args.index is not None:
        cl = clusters[args.index - 1]
    elif args.cluster:
        cl = args.cluster
    else:
        raise SystemExit("need --cluster or --index")
    db_pwms = _load_meme_pwms(DEFAULT_VIERSTRA_MEME); _set_current_db_pwms(db_pwms)
    pwms = _load_meme_pwms(str(INP/"query.meme"))
    marg = pd.read_csv(ANA/"marginalization_matrix.csv", index_col=0)
    margz = pd.read_csv(ANA/"marginalization_matrix_zscore.csv", index_col=0)
    cpm = pd.read_csv(CPM_CSV).set_index("cell_type")
    meta = pd.read_csv(META_TSV, sep="\t").sort_values("display_order")
    ct_order = [c for c in meta["cell_type"].tolist() if c in marg.columns and c in cpm.index]
    ct_colors = dict(zip(meta["cell_type"], meta["color"])); cpm = cpm.loc[ct_order]
    expressed = {tf for tf in cpm.columns if float(cpm[tf].max()) >= MIN_EXPR_CPM}
    ctx = load_context(xlsx_path="/cellar/users/aklie/projects/islet_organoid_differentiation/sc-islet-differentiation_10X-Multiome/ref/motif_annotations.xlsx",
                       pseudobulk_csv=CPM_CSV, metadata_tsv=META_TSV,
                       alias_json="/cellar/users/aklie/data/datasets/sc-islet-differentiation_10X-Multiome/results/2_process_data/vierstra_alias_map.json",
                       cluster_viewer_base="#")
    shared = (marg, margz, cpm, expressed, ct_order, ct_colors, ctx, db_pwms, pwms)
    out = build_cluster_sheet(cl, manifest[cl], shared)
    print(f"Wrote {out} ({out.stat().st_size:,} bytes, {len(manifest[cl]['constituents'])+1} pages)")


if __name__ == "__main__":
    main()

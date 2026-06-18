#!/usr/bin/env python
"""Split the constituent query.meme into per-motif MEME files + per-cluster motif
lists (for the TomTom array). One query per file (multi-query -text only annotates
the first motif).

Reads:  .../gimme_cluster/marginalization/constituent/{query.meme, manifest.json}
Writes: .../constituent/per_motif/<name>.meme
        .../constituent/cluster_motifs/<cluster>.txt
        .../constituent/clusters.txt
"""
import json
from pathlib import Path

CONST = Path("/cellar/users/aklie/data/datasets/sc-islet-differentiation_10X-Multiome/results/3_single_task_models/motifs/gimme_cluster/marginalization/constituent")
PM = CONST/"per_motif"; PM.mkdir(parents=True, exist_ok=True)
CM = CONST/"cluster_motifs"; CM.mkdir(parents=True, exist_ok=True)
HEADER = ("MEME version 4\n\nALPHABET= ACGT\n\nstrands: + -\n\n"
          "Background letter frequencies\nA 0.25 C 0.25 G 0.25 T 0.25\n\n")

# parse global meme into {name: block_text}
blocks = {}
cur, buf = None, []
for line in open(CONST/"query.meme"):
    if line.startswith("MOTIF "):
        if cur:
            blocks[cur] = "".join(buf)
        cur = line.split()[1]; buf = [line]
    elif cur:
        buf.append(line)
if cur:
    blocks[cur] = "".join(buf)
for name, body in blocks.items():
    (PM/f"{name}.meme").write_text(HEADER + body if not body.startswith("MOTIF") else HEADER + body)

manifest = json.load(open(CONST/"manifest.json"))
clusters = list(manifest.keys())
(CONST/"clusters.txt").write_text("\n".join(clusters) + "\n")
for cl, info in manifest.items():
    names = [c["name"] for c in info["constituents"]]
    if info.get("consensus_name"):
        names.append(info["consensus_name"])
    (CM/f"{cl}.txt").write_text("\n".join(names) + "\n")
print(f"per-motif MEMEs: {len(blocks)}")
print(f"clusters: {len(clusters)}")
print(f"example cluster motif list ({clusters[0]}):")
print((CM/f'{clusters[0]}.txt').read_text().strip())

"""Render motif walls for the dual-track gimme cluster catalog.

For each motif in ``motifs_all_peaks/combined_mapped.meme`` (149 motifs:
62 pos_ + 87 neg_), build a LogoPlottingInput and pass to
MotifCompendium's ``motif_collection_html`` to produce one interactive
HTML with the pos and neg tracks side by side.

This is purely a visualization of the GIMME CATALOG as-shipped — useful
for the user's diagnostic ("the neg track of gimme cluster didn't seem to
have any neg motifs in the final output") since the values stored are PPMs
(non-negative by construction; sign info is only in the name prefix).
"""

from __future__ import annotations

import os
from pathlib import Path
import numpy as np

CUDA_ENV = "/cellar/users/aklie/opt/miniconda3/envs/motifcompendium-gpu"
os.environ["CUDA_PATH"] = CUDA_ENV
os.environ["CPLUS_INCLUDE_PATH"] = f"{CUDA_ENV}/include"
os.environ["CUPY_NVRTC_OPTIONS"] = f"-I{CUDA_ENV}/include"

import cupy  # noqa
_ = cupy.ones((10,), dtype=cupy.float16); del _

from MotifCompendium.utils.plotting import LogoPlottingInput
from MotifCompendium.utils.visualization import motif_collection_html

BASE = Path("/cellar/users/aklie/data/datasets/sc-islet-differentiation_10X-Multiome")
CATALOG = BASE / "results/4_multi_task_model/crested/motifs_all_peaks/combined_mapped.meme"
OUT_HTML = BASE / "results/4_multi_task_model/crested/motifs_all_peaks/catalog_walls.html"


def parse_meme(path: Path) -> list[tuple[str, np.ndarray]]:
    """Yield (motif_name, PWM array of shape (L, 4)) from a MEME catalog."""
    motifs: list[tuple[str, np.ndarray]] = []
    cur_name: str | None = None
    in_matrix = False
    rows: list[list[float]] = []
    with open(path) as fh:
        for raw in fh:
            line = raw.rstrip()
            if line.startswith("MOTIF"):
                if cur_name and rows:
                    motifs.append((cur_name, np.asarray(rows, dtype=np.float32)))
                cur_name = line.split()[1]
                rows = []
                in_matrix = False
            elif line.startswith("letter-probability matrix"):
                in_matrix = True
            elif in_matrix:
                if not line.strip():
                    in_matrix = False
                    continue
                try:
                    parts = [float(x) for x in line.split()]
                except ValueError:
                    in_matrix = False
                    continue
                if len(parts) >= 4:
                    rows.append(parts[:4])
                else:
                    in_matrix = False
        if cur_name and rows:
            motifs.append((cur_name, np.asarray(rows, dtype=np.float32)))
    return motifs


def main():
    print(f"Parsing {CATALOG}...")
    motifs = parse_meme(CATALOG)
    print(f"  parsed {len(motifs)} motifs")

    pos_inputs: list[LogoPlottingInput] = []
    neg_inputs: list[LogoPlottingInput] = []
    for name, pwm in motifs:
        # Use a distinct bgcolor per track so the user sees pos/neg at a glance
        if name.startswith("pos_"):
            pos_inputs.append(LogoPlottingInput(motif=pwm, name=name, bgcolor="#e7f3e7"))
        elif name.startswith("neg_"):
            neg_inputs.append(LogoPlottingInput(motif=pwm, name=name, bgcolor="#fce5e5"))
        else:
            # shouldn't happen — every catalog motif is sign-prefixed
            print(f"  WARN: motif without pos_/neg_ prefix: {name}")

    print(f"  pos track: {len(pos_inputs)} motifs")
    print(f"  neg track: {len(neg_inputs)} motifs")

    groups = {
        f"POS track ({len(pos_inputs)} motifs, green) — activating motifs from modisco pos_patterns": pos_inputs,
        f"NEG track ({len(neg_inputs)} motifs, pink) — modisco neg_patterns; "
        f"NOTE values are PPMs (sign-agnostic); see motif_compendium_all_peaks/qc/ for posneg_inverted diagnostic":
            neg_inputs,
    }

    print(f"Rendering → {OUT_HTML}")
    motif_collection_html(groups, str(OUT_HTML))
    sz = OUT_HTML.stat().st_size
    print(f"  done: {sz:,} bytes ({sz/1e6:.1f} MB)")


if __name__ == "__main__":
    main()

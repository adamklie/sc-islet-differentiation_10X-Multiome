
#!/usr/bin/env python

# -------------------------
# Imports
# -------------------------
import os, sys, csv, random, hashlib
import torch
import heapq
import numpy as np
import pandas as pd
import seaborn as sns
from tqdm.auto import tqdm
from scipy.stats import spearmanr, pearsonr
import pyfaidx
import matplotlib.pyplot as plt

import seqpro as sp
from bpnetlite.bpnet import BPNet, CountWrapper
from tangermeme.io import extract_loci, read_meme
from tangermeme.predict import predict
from tangermeme.deep_lift_shap import deep_lift_shap
from tangermeme.plot import plot_logo
from tangermeme.seqlet import recursive_seqlets
from tangermeme.annotate import annotate_seqlets
from tangermeme.utils import _cast_as_tensor, random_one_hot

# -------------------------
# Helpers
# -------------------------
def create_run_hash():
    random_str = ''.join(random.choices('abcdefghijklmnopqrstuvwxyz0123456789', k=10))
    return hashlib.md5(random_str.encode()).hexdigest()[:10]

def log_softmax_profile(y_profile):
    z = y_profile.shape
    y_profile = y_profile.reshape(y_profile.shape[0], -1)
    y_profile = torch.nn.functional.log_softmax(y_profile, dim=-1)
    return y_profile.reshape(*z)

def intersects_with_peak(row, peaks):
    """Return True if enhancer row overlaps with any peak."""
    chrom = row['chrom']
    start = row['start']
    end = row['end']
    chrom_peaks = peaks[peaks["chrom"] == chrom]
    for _, peak in chrom_peaks.iterrows():
        if not (end < peak["start"] or start > peak["end"]):
            return True
    return False

def get_peak_coordinates(row, peaks):
    """Return first overlapping peak coordinates for enhancer row."""
    chrom = row['chrom']
    start = row['start']
    end = row['end']
    chrom_peaks = peaks[peaks["chrom"] == chrom]
    for _, peak in chrom_peaks.iterrows():
        if not (end < peak["start"] or start > peak["end"]):
            return f"{chrom}:{peak['start']}-{peak['end']}"
    return None

def overlap_len(a_start, a_end, b_start, b_end):
    return max(0, min(a_end, b_end) - max(a_start, b_start))

def greedy_match(curr_design_df, baseline_df):
    """Greedy matching of seqlets between design and baseline."""
    cand = (
        curr_design_df.reset_index().merge(
            baseline_df.reset_index(),
            on="annotation",
            suffixes=("_d","_b")
        )
    )
    if cand.empty:
        return [], set(), set()

    cand["overlap"] = cand.apply(
        lambda r: overlap_len(r["start_d"], r["end_d"], r["start_b"], r["end_b"]),
        axis=1
    )
    cand = cand[cand["overlap"] > 0]
    if cand.empty:
        return [], set(), set()

    cand["start_gap"] = (cand["start_d"] - cand["start_b"]).abs()
    cand = cand.sort_values(["overlap","start_gap"], ascending=[False,True])

    used_d, used_b, matches = set(), set(), []
    for _, r in cand.iterrows():
        di, bi = r["index_d"], r["index_b"]
        if di not in used_d and bi not in used_b:
            used_d.add(di); used_b.add(bi)
            matches.append((di, bi))
    return matches, used_d, used_b

def screen_insertions(model, shape, seq_start, seq_end, y=None,
                      loss=torch.nn.MSELoss(reduction='none'),
                      max_iter=100, args=None, n_best=1,
                      alphabet=['A', 'C', 'G', 'T'], 
                      batch_size=32, func=None, additional_func_kwargs={},
                      dtype=None, device='cuda', random_state=None,
                      verbose=False, background=None):
    """
    Screen random insertions between seq_start and seq_end.

    Runs for exactly `max_iter` iterations (each with `batch_size` samples),
    and returns the top-k sequences with lowest loss.

    Parameters
    ----------
    model : torch.nn.Module
        Model for predictions.
    shape : tuple
        Shape of the full one-hot sequence (C, L).
    seq_start : int
        Start coordinate (inclusive) of the insertion span.
    seq_end : int
        End coordinate (exclusive) of the insertion span.
    y : torch.Tensor
        Target output for the loss function.
    max_iter : int
        Number of batches to screen.
    n_best : int
        How many best sequences to keep.
    background : torch.Tensor
        Template sequence (C, L). Insertions replace seq_start:seq_end.

    Returns
    -------
    X : torch.Tensor, shape=(n_best, C, L)
        Top-k screened sequences.
    """

    if y is None:
        y = [9_999_999]
    y = _cast_as_tensor(y)

    if background is None:
        background = torch.zeros(shape, dtype=torch.float32)

    pq = []  # priority queue for top-k

    for iteration in tqdm(range(max_iter)):
        # generate random insertions of length (seq_end - seq_start)
        insert_shape = (batch_size, shape[0], seq_end - seq_start)
        if func is None:
            X_insert = random_one_hot(insert_shape, random_state=random_state,
                                      **additional_func_kwargs)
        else:
            X_insert = func(insert_shape, random_state=random_state,
                            **additional_func_kwargs)

        # place insertions into background template
        X = background.unsqueeze(0).repeat(batch_size, 1, 1).clone()
        X[:, :, seq_start:seq_end] = X_insert

        # model predictions
        y_hat = predict(model, X, dtype=dtype, device=device, args=args)

        # negative loss so bigger = better
        current_loss = -loss(y, y_hat).numpy(force=True).sum(axis=1)

        for i in range(batch_size):
            l = current_loss[i]
            entry = [l, X[i]]

            if len(pq) < n_best:
                heapq.heappush(pq, entry)
            elif l > pq[0][0]:
                heapq.heappushpop(pq, entry)

            if verbose:
                print(f"[iter {iteration}] adding seq with score {l:.4f}")

        if random_state is not None:
            random_state += 1

    # return top-k
    Xs = [heapq.heappop(pq)[1] for _ in range(len(pq))]
    X = torch.flip(torch.stack(Xs), dims=(0,))
    return X

# -------------------------
# Config
# -------------------------
print("Setting up environment and loading data...")

os.chdir("/cellar/users/aklie/projects/igvf/sc-islet-design")
path_out = "/cellar/users/aklie/projects/igvf/sc-islet-design/scratch/2025_10_04"

# Genome + peaks
genome_fasta = "/cellar/users/aklie/data/ref/genomes/hg38/hg38.fa"
genome = pyfaidx.Fasta(genome_fasta)

path_peaks = "/cellar/users/aklie/projects/igvf/sc-islet-design/scratch/2025_09_05/peaks_with_predictions.tsv"
peaks = pd.read_csv(path_peaks, sep="\t")
peaks = peaks[peaks["log_counts_pred"] >= 0]
print(f"Loaded {len(peaks)} peaks")

# Motifs
path_motifs = "/cellar/users/aklie/data/datasets/sc-islet-differentiation_10X-Multiome/models/motifs/meme/combined.meme"
motifs = read_meme(path_motifs)
motif_names = np.array(list(motifs.keys()))
print(f"Loaded {len(motifs)} motifs")

path_annot = "/cellar/users/aklie/data/datasets/sc-islet-differentiation_10X-Multiome/models/motifs/tfs_initial.txt"
tmp = pd.read_csv(path_annot, header=None, sep='\t', names=["motif_id", "annotated_tf", "evalue"])
motif_id_mp = tmp.set_index("motif_id")["annotated_tf"]

# Model
path_model = "/cellar/users/aklie/data/datasets/sc-islet-differentiation_10X-Multiome/models/D15_ENP/fold_0/chrombpnet/0.5/models/chrombpnet_nobias.h5"
model = BPNet.from_chrombpnet(path_model).cuda().eval()
count_wrapper = CountWrapper(model).cuda().eval()
print("Model loaded and moved to GPU")

# MPRA scaffold
mpra_fasta = "/cellar/users/aklie/projects/igvf/sc-islet-design/scratch/2025_07_05/MPRA_integrated_sequence.fa"
mpra_seqs = pyfaidx.Fasta(mpra_fasta)
MPRA_integrated_sequence_1 = str(mpra_seqs["MPRA_integrated_sequence_1"][:].seq).upper()
print(f"Loaded MPRA scaffold, length {len(MPRA_integrated_sequence_1)}")

# Spans file
path_seqlet_spans = "/cellar/users/aklie/projects/igvf/sc-islet-design/scratch/2025_09_25/merged_seqlet_spans.tsv"
seqlet_span_midpoints = pd.read_csv(path_seqlet_spans, sep="\t")["midpoint"].values[:5]
print(f"Loaded seqlet span midpoints: {seqlet_span_midpoints}")

# -------------------------
# Parameters
# -------------------------
flank_choices = list(range(4, 11))   # 4 to 10 inclusive
span_choices = seqlet_span_midpoints
dtype = torch.float64
k = 5
max_iter = 100
batch_size = 64
seq_len = 2114
target_len = 1000
ism_window = 270
plot_window = 270
seq_center = 2114 // 2
seq_start = seq_center - ism_window // 2
seq_end = seq_center + ism_window // 2
target_center = 1000 // 2
target_start = target_center - ism_window // 2
target_end = target_center + ism_window // 2
plot_seq_center = 1057
plot_seq_start = seq_center - plot_window // 2
plot_seq_end = seq_center + plot_window // 2
plot_target_center = 500
plot_target_start = target_center - plot_window // 2
plot_target_end = target_center + plot_window // 2
print("\nParameters set.")
print(f"\tseq_start: {seq_start}, seq_end: {seq_end}")
print(f"\ttarget_start: {target_start}, target_end: {target_end}")
print(f"\tplot_seq_start: {plot_seq_start}, plot_seq_end: {plot_seq_end}")
print(f"\tplot_target_start: {plot_target_start}, plot_target_end: {plot_target_end}")

# -------------------------
# Candidate
# -------------------------
path_fasta = "/cellar/users/aklie/projects/igvf/sc-islet-design/scratch/2025_07_05/ISL1_CRE.fa"
seqs = pyfaidx.Fasta(path_fasta)
seq_headers = list(seqs.keys())
isl1_enhancer_df = pd.DataFrame(seq_headers, columns=['header'])
isl1_enhancer_df['chrom'] = isl1_enhancer_df['header'].str.split(':').str[0]
isl1_enhancer_df['start'] = isl1_enhancer_df['header'].str.split(':').str[1].str.split('-').str[0].astype(int)
isl1_enhancer_df['end'] = isl1_enhancer_df['header'].str.split(':').str[1].str.split('-').str[1].astype(int)
isl1_enhancer_df['intersects_peak'] = isl1_enhancer_df.apply(intersects_with_peak, peaks=peaks, axis=1)
isl1_enhancer_df['peak_coordinates'] = isl1_enhancer_df.apply(get_peak_coordinates, peaks=peaks, axis=1)
isl1_enhancer_df["peak_chrom"] = isl1_enhancer_df["peak_coordinates"].str.split(":").str[0]
isl1_enhancer_df["peak_start"] = isl1_enhancer_df["peak_coordinates"].str.split(":").str[1].str.split("-").str[0].astype(int)
isl1_enhancer_df["peak_end"] = isl1_enhancer_df["peak_coordinates"].str.split(":").str[1].str.split("-").str[1].astype(int)
peaks_selected = peaks.merge(isl1_enhancer_df[["peak_chrom", "peak_start", "peak_end"]], 
                                       left_on=["chrom", "start", "end"], 
                                       right_on=["peak_chrom", "peak_start", "peak_end"], 
                                       how="left", suffixes=("", "_enhancer"))
peaks_selected['CRE_annotation'] = np.where(peaks_selected['peak_chrom'].notnull(), 'Enhancer', np.nan)
peaks_selected['CRE_annotation'] = peaks_selected.groupby('chrom').cumcount() + 1
peaks_selected['CRE_annotation'] = np.where(peaks_selected['peak_chrom'].notnull(), 'CRE_' + peaks_selected['CRE_annotation'].astype(str), np.nan)
peaks_selected = peaks_selected.drop(columns=["peak_chrom", "peak_start", "peak_end"])
deciles = np.percentile(peaks_selected["log_counts_pred"].dropna(), [0, 5, 15, 25, 35, 45, 55, 65, 75, 85, 95, 100])
peaks_selected['decile'] = pd.cut(peaks_selected['log_counts_pred'], bins=deciles, labels=False)
design_res = {}
CRE_ = peaks_selected[peaks_selected["CRE_annotation"].notnull()].iloc[2]
idx_ = peaks_selected[peaks_selected["CRE_annotation"].notnull()].index[2]
chrom_ = CRE_["chrom"]
start_ = CRE_["start"]
end_ = CRE_["end"]
summit_ = CRE_["peak"]
decile_ = CRE_["decile"]
preds_ = peaks_selected["log_counts_pred"].values[idx_]
ref_seq_start = start_ + summit_ - seq_len // 2
ref_seq_end = ref_seq_start + seq_len
design_res["reference_locus"] = f"{chrom_}:{ref_seq_start}-{ref_seq_end}"
ref_seq = genome[chrom_][ref_seq_start:ref_seq_end].seq.upper()
design_res["reference_sequence"] = ref_seq
X_ref = torch.tensor(sp.ohe([ref_seq], alphabet=sp.DNA).transpose(0, 2, 1), dtype=torch.uint8)
y_ref = predict(
    model=count_wrapper.type(torch.float32),
    X=X_ref,
    verbose=False,
)
design_res["reference_prediction"] = y_ref.item()
deciles_above = deciles[int(decile_)+1:]
plot_coordinate_start = start_ + summit_ - plot_window // 2
plot_coordinate_end = plot_coordinate_start + plot_window
print(f"Deciles above {deciles[int(decile_)]}: {deciles_above}")
print(f"Processing CRE: {chrom_}:{start_}-{end_}: {preds_}")
print(f"Plotting coordinates: {chrom_}:{plot_coordinate_start}-{plot_coordinate_end}")

# Output dir
run_hash = create_run_hash()
dir_name = f"ISL1_{chrom_}:{start_}-{end_}_screen_design_{run_hash}"
os.makedirs(os.path.join(path_out, dir_name), exist_ok=True)
print(f"Created output directory: {os.path.join(path_out, dir_name)}")

# Target
y_bars = [torch.zeros(1, dtype=torch.float32) + d for d in deciles_above][-1]
print(f"Design targets set to {y_bars}")

# -------------------------
# Sweep
# -------------------------
if __name__ == "__main__":
    run_hash = create_run_hash()
    print(f"\nStarting sweep with run_hash={run_hash}")
    count = 0
    for midpoint in seqlet_span_midpoints:
        for fl in flank_choices:
            print(f"\nDesign {count+1}/{len(seqlet_span_midpoints)*len(flank_choices)}: span={midpoint}, flanks={fl}")

            # Copy dict
            design_res_curr = design_res.copy()
            
            # Edit positions
            edit_start = midpoint - fl
            edit_end = midpoint + fl

            # Create directory to save
            subdir_name = f"start-{edit_start}_end-{edit_end}"
            os.makedirs(os.path.join(path_out, dir_name, subdir_name), exist_ok=True)

            # Create boolean tensor mask where the flanks around midpoint are 0 and the rest are 1
            msk = torch.ones(1, seq_len, dtype=torch.bool)
            msk[0, edit_start:edit_end] = 0
            msk = msk.reshape(-1)
            print(f"Masking positions {edit_start}-{edit_end}: {msk.sum()}")

            # Run screen
            X_hat = []
            for y_bar in tqdm(y_bars, desc="Running screening for each decile"):
                print(f"Running screening for decile target: {y_bar.item()} "
                    f"in positions ({edit_start}-{edit_end})")
                X_best = screen_insertions(
                    model=count_wrapper,
                    shape=X_ref.shape[1:],            # drop batch dim
                    seq_start=edit_start,
                    seq_end=edit_end,
                    y=y_bar,
                    n_best=k,                        # return top-k hits
                    max_iter=max_iter,               # number of batches to search
                    batch_size=batch_size,
                    device="cuda",
                    verbose=False,
                    background=X_ref.squeeze(0)       # template sequence
                )

                # shape: (k, 4, L)
                X_hat.append(X_best)

            # stack into (len(y_bars), k, 4, L)
            X_hat = torch.stack(X_hat)
            print(f"X_hat shape: {X_hat.shape}")
            
            # Create identifiers for each of X_hat based on ISL1_{chrom_}:{start_}-{end_}_screen_run_{run}_decile_{decile}
            identifiers = []
            for decile in range(X_hat.shape[0]):
                for t, seq in enumerate(range(X_hat.shape[1])):
                    identifiers.append(f"ISL1_{chrom_}:{start_}-{end_}_screen_target-{round(y_bars[decile].item(), 2)}_pos-{midpoint}_flanks-{fl}_{t+1}")
            print(f"First ID: {identifiers[0]}")

            # Stack the everything but the last two dims
            X_optim = X_hat.view(-1, X_hat.shape[-2], X_hat.shape[-1])
            print(f"X_optim shape: {X_optim.shape}")

            # Prdict and plot
            y_hat = []
            for j, y_bar in enumerate(y_bars):
                y0 = y_ref.numpy(force=True)[0, 0]
                X_hat_ = X_optim[j*k:(j+1)*k]
                ax = plt.subplot(111)
                y_hat_ = predict(
                    model=count_wrapper.type(torch.float32),
                    X=X_hat_,
                    verbose=True,
                )
                plt.hist(y_hat_, facecolor='0.5', edgecolor='0.5')
                y_hat.append(y_hat_)
                plt.plot([y0, y0], [0, 5], label="Baseline prediction", c='c')
                plt.plot([y_bar, y_bar], [0, 5], label="Design target", c='m', linestyle=(0, (2, 2)))
                plt.title(f"ISM for {chrom_}:{start_}-{end_} with target {round(y_bar.item(), 2)}")
                plt.legend(loc=(1.01, 0.5))
                plt.xlabel("Predicted Log Counts")
                plt.ylabel("# Sequences")
                sns.despine(left=True)
                plt.grid(False)
                plt.setp(ax.spines.values(), linewidth=2, color='0.5')
                plt.tight_layout()
                plt.savefig(os.path.join(path_out, dir_name, subdir_name, f"ISL1_{chrom_}:{start_}-{end_}_screen_design_target_{round(y_bar.item(), 2)}_{run_hash}.pdf"), dpi=300)
                plt.show()
            y_hat = torch.stack(y_hat).squeeze()
            y_hat = y_hat.view(-1)

            # Grab actual string sequences
            optim_seqs = sp.decode_ohe(X_optim.detach().cpu().numpy().astype(np.uint8), alphabet=sp.DNA, ohe_axis=1)
            optim_seqs = [b''.join(optim_seq).decode('utf-8') for optim_seq in optim_seqs]

            # Get attribution scores for the original sequence
            X_ref_attr_ = deep_lift_shap(
                count_wrapper.type(dtype), 
                X=X_ref.type(dtype),
                n_shuffles=20,
                batch_size=16,
                verbose=False,
                warning_threshold=1e-4
            )
            print(f"X_ref_attr_ shape: {X_ref_attr_.shape}")

            # Get attribution scores for the original sequence
            X_optim_attr_ = deep_lift_shap(
                count_wrapper.type(dtype), 
                X=X_optim.type(dtype),
                n_shuffles=20,
                batch_size=16,
                verbose=False,
                warning_threshold=1e-4
            )
            print(f"X_optim_attr_ shape: {X_optim_attr_.shape}")

            # Save attribution scores to npz
            np.savez_compressed(os.path.join(path_out, dir_name, subdir_name, f"ISL1_{chrom_}:{start_}-{end_}_screen_X_ref_attr_pos-{midpoint}_flanks-{fl}_{run_hash}.npz"), X_ref_attr=X_ref_attr_.numpy(force=True))
            np.savez_compressed(os.path.join(path_out, dir_name, subdir_name, f"ISL1_{chrom_}:{start_}-{end_}_screen_X_optim_attr_pos-{midpoint}_flanks-{fl}_{run_hash}.npz"), X_optim_attr=X_optim_attr_.numpy(force=True))

            # Call and annotate seqlets
            X_attr_ = torch.cat([X_ref_attr_, X_optim_attr_], dim=0)
            total_attr_mp = pd.Series(index=["ISL1_chr5:51382751-51383709_ref_sequence"] + identifiers, data=torch.abs(X_attr_.detach()).sum(dim=(-2, -1)).numpy())
            X_ = torch.cat([X_ref.cpu(), X_optim.cpu()], dim=0)
            seqlets = recursive_seqlets(X_attr_[:, :, seq_start:seq_end].sum(dim=1).detach(), threshold=0.05, additional_flanks=5)
            motif_idxs = annotate_seqlets(X_[:, :, seq_start:seq_end].detach().cpu(), seqlets, path_motifs)[0][:, 0]
            seqlets['start'] += seq_start
            seqlets['end'] += seq_start
            motifs_ = motif_id_mp.loc[motif_names[motif_idxs]]
            seqlets["annotation"] = motifs_.values
            seqlets["id"] = [(["ISL1_chr5:51382751-51383709_ref_sequence"] + identifiers)[idx] for idx in seqlets["example_idx"]]

            # calculate the proportion of total attribution for each seqlet
            seqlets["attribution_proportion"] = seqlets["attribution"].abs() / seqlets["id"].map(total_attr_mp)
            print(f"Seqlets found: {len(seqlets)}")

            # Save
            seqlets[["id", "start", "end", "attribution", "p-value", "annotation"]].sort_values(["id", "start"]).to_csv(
                os.path.join(path_out, dir_name, subdir_name, f"ISL1_{chrom_}:{start_}-{end_}_screen_seqlets_pos-{midpoint}_flanks-{fl}_{run_hash}.tsv"), 
                sep="\t", 
                index=False
            )

            ref_seqlets = seqlets[seqlets["id"].str.contains("_ref_sequence")].copy()
            designs_seqlets  = seqlets[~seqlets["id"].str.contains("_ref_sequence")].copy()
            
            # Build delta-seqlets once (designs_seqlets vs ref_seqlets)
            rows = []
            for design_id, curr_design_df in designs_seqlets.groupby("id", sort=False):
                # work on copies to avoid SettingWithCopy surprises
                d = curr_design_df[["start","end","attribution","attribution_proportion","p-value","annotation"]].copy()
                b = ref_seqlets[      ["start","end","attribution","attribution_proportion","p-value","annotation"]].copy()

                # ensure numeric
                d["start"] = d["start"].astype(int);  d["end"] = d["end"].astype(int);  d["attribution"] = d["attribution"].astype(float)
                b["start"] = b["start"].astype(int);  b["end"] = b["end"].astype(int);  b["attribution"] = b["attribution"].astype(float)

                matches, used_d, used_b = greedy_match(d, b)

                # matched → increase/decrease
                for di, bi in matches:
                    drow, brow = d.loc[di], b.loc[bi]
                    delta_attr = float(drow["attribution"]) - float(brow["attribution"])
                    change_type = "increase" if delta_attr > 0 else ("decrease" if delta_attr < 0 else "increase")
                    rows.append({
                        "id": design_id,
                        "start": int(drow["start"]),
                        "end": int(drow["end"]),
                        "attribution": float(drow["attribution"]),
                        "p-value": float(drow["p-value"]),
                        "attribution_proportion": float(drow["attribution_proportion"]),
                        "delta_attribution": float(delta_attr),
                        "type": change_type,
                        "annotation": drow["annotation"],
                    })

                # additions → present in design only
                if len(d) > 0 and len(used_d) > 0:
                    d_unmatched = d.drop(index=list(used_d))
                else:
                    d_unmatched = d if len(d) else d.iloc[0:0]
                for _, drow in d_unmatched.iterrows():
                    rows.append({
                        "id": design_id,
                        "start": int(drow["start"]),
                        "end": int(drow["end"]),
                        "attribution": float(drow["attribution"]),
                        "p-value": float(drow["p-value"]),
                        "attribution_proportion": float(drow["attribution_proportion"]),
                        "delta_attribution": float(drow["attribution"]),
                        "type": "addition",
                        "annotation": drow["annotation"],
                    })

                # removals → present in ref only
                if len(b) > 0 and len(used_b) > 0:
                    b_unmatched = b.drop(index=list(used_b))
                else:
                    b_unmatched = b if len(b) else b.iloc[0:0]
                for _, brow in b_unmatched.iterrows():
                    rows.append({
                        "id": design_id,
                        "start": int(brow["start"]),
                        "end": int(brow["end"]),
                        "attribution": float("nan"),
                        "p-value": float("nan"),
                        "attribution_proportion": float("nan"),
                        "delta_attribution": -float(brow["attribution"]),
                        "type": "removal",
                        "annotation": brow["annotation"],
                    })

            delta_seqlets_all = pd.DataFrame(rows)
            delta_seqlets_all["type"] = pd.Categorical(
                delta_seqlets_all["type"],
                categories=["increase","decrease","addition","removal"],
                ordered=True,
            )
            
            
            # Save
            delta_seqlets_all.sort_values(["id", "start"]).to_csv(
                os.path.join(path_out, dir_name, subdir_name, f"ISL1_{chrom_}:{start_}-{end_}_screen_delta_seqlets_pos-{midpoint}_flanks-{fl}_{run_hash}.tsv"), 
                sep="\t", 
                index=False
            )
            
            ## Quantify edits
            
            # Diffs is a list of lists of tuples and we need to match up each list with id
            diffs = []
            for j in range(len(optim_seqs)):
                diffs_ = [(i, ref_seq[i], optim_seqs[j][i]) for i in range(len(ref_seq)) if ref_seq[i] != optim_seqs[j][i]]
                diffs.append(diffs_)

            # Make df
            edit_df = []
            for i, diff_list in enumerate(diffs):
                id_ = identifiers[i]
                for pos, ref, alt in diff_list:
                    edit_df.append({
                        'id': id_,
                        'relative_pos': pos,
                        'ref': ref,
                        'alt': alt
                    })

            edit_df = pd.DataFrame(edit_df)
            edit_df["chrom"] = chrom_
            edit_df["pos"] = edit_df["relative_pos"] + ref_seq_start + 1
            edit_df["variant_id"] = edit_df["chrom"] + ":" + edit_df["pos"].astype(str) + ":" + edit_df["ref"] + ":" + edit_df["alt"]
            edit_df["id"].value_counts()
            
            # Save edit_df
            edit_df[["id", "chrom", "pos", "ref", "alt", "variant_id", "relative_pos"]].to_csv(os.path.join(path_out, dir_name, subdir_name, f"ISL1_{chrom_}:{start_}-{end_}_screen_edits_pos-{midpoint}_flanks-{fl}_{run_hash}.tsv"), sep="\t", index=False)
            
            # Predict on MPRA integrated sequence context
            insert_len = 270
            barcode_len = 15
            insert_seqs = [optim_seq[seq_center - insert_len // 2: seq_center + insert_len // 2] for optim_seq in optim_seqs]
            results = []

            # For each insert sequence, create the full MPRA sequence by inserting into MPRA_integrated_sequence_1, then chop at the middle 2114 and predict
            for i, insert_seq in enumerate(insert_seqs):

                # Get Id
                curr_id = identifiers[i]

                # Get a random 15bp barcode sequence
                barcode_seq = ''.join(random.choices('ACGT', k=barcode_len))

                # Find the start and end of the 270 'N's in a row
                insert_start = MPRA_integrated_sequence_1.find('N' * insert_len)
                insert_end = insert_start + insert_len
                MPRA_seq = MPRA_integrated_sequence_1[:insert_start] + insert_seq + MPRA_integrated_sequence_1[insert_end:]

                # Find the start and end of exactly 15 'N's in a row and no more than 15 N's
                barcode_start = MPRA_seq.find('N' * barcode_len)
                barcode_end = barcode_start + barcode_len
                MPRA_seq = MPRA_seq[:barcode_start] + barcode_seq + MPRA_seq[barcode_end:]

                # Create seq by inserting insert_seq and barcode_seq into MPRA_integrated_sequence_1 at proper positions
                MPRA_seq = MPRA_seq.upper()

                # Take the middle seq_len of MPRA_seq
                insert_center = insert_start + insert_len // 2
                mpra_seq_start = insert_center - seq_len // 2
                mpra_seq_end = insert_center + seq_len // 2
                mpra_seq = MPRA_seq[mpra_seq_start:mpra_seq_end]

                # One-hot encode and predict
                X_mpra = torch.tensor(sp.ohe([mpra_seq], alphabet=sp.DNA).transpose(0, 2, 1), dtype=torch.uint8)
                y_mpra = predict(
                    model=count_wrapper.type(torch.float32),
                    X=X_mpra,
                    verbose=False
                )

                # Save results
                results.append({
                    "id": curr_id,
                    "mpra_sequence": mpra_seq,
                    "mpra_prediction": float(y_mpra.item())
                })

            # Collect into dataframe
            mpra_df = pd.DataFrame(results)
            
            # Save overall results
            design_df = pd.DataFrame({
                "id": identifiers,
                "reference_locus": design_res["reference_locus"],
                "reference_sequence": design_res["reference_sequence"],
                "reference_prediction": design_res["reference_prediction"],
                "design_method": "screen",
                "design_target": np.concatenate([np.repeat(y_bar.item(), k) for y_bar in y_bars]),
                "design_sequence": optim_seqs,
                "design_prediction": y_hat.numpy(force=True),
            })

            # merge designed_df with edited_df on id
            designed_merged_df = design_df.merge(mpra_df, on="id")
            
            # Overwrite design_sequence and design_prediction with edited versions
            header = {
                "method": "screen",
            }
            with open(os.path.join(path_out, dir_name, subdir_name, f"ISL1_{chrom_}:{start_}-{end_}_screen_design_{run_hash}.tsv"), 'w', newline='') as f:
                writer = csv.writer(f, delimiter='\t')
                for key, value in header.items():
                    writer.writerow([f"# {key}: {value}"])
                designed_merged_df.to_csv(f, sep="\t", index=False, float_format='%.4f')
            
            # Save X_optim to npz
            np.savez_compressed(os.path.join(path_out, dir_name, subdir_name, f"ISL1_{chrom_}:{start_}-{end_}_screen_X_optim_pos-{midpoint}_flanks-{fl}_{run_hash}.npz"), X_optim=X_optim.numpy(force=True))
            
            # Plot them all
            plots_dir = os.path.join(path_out, dir_name, subdir_name, "plots")
            os.makedirs(plots_dir, exist_ok=True)
            for seq_num in range(len(optim_seqs)):
                print(f"Plotting sequence {seq_num+1}/{len(optim_seqs)}")

                # --- Get sequences ---
                original_seq = sp.decode_ohe(
                    X_ref.cpu().numpy().astype(np.uint8), 
                    alphabet=sp.DNA, ohe_axis=1
                )
                original_seq = b''.join(original_seq[0]).decode('utf-8')

                mutated_seq = sp.decode_ohe(
                    np.expand_dims(X_optim[seq_num].detach().cpu().numpy().astype(np.uint8), axis=0),
                    alphabet=sp.DNA, ohe_axis=1
                )
                mutated_seq = b''.join(mutated_seq[0]).decode('utf-8')

                diffs = [
                    (i, original_seq[i], mutated_seq[i]) 
                    for i in range(len(original_seq)) if original_seq[i] != mutated_seq[i]
                ]

                # --- Get predictions ---
                original_profile, original_counts = predict(model, X_ref)
                original_counts = np.squeeze(np.exp(original_counts.cpu().detach().numpy()))
                original_profile = np.squeeze(torch.exp(log_softmax_profile(original_profile)).cpu().detach().numpy())

                mutated_profile, mutated_counts = predict(model, X_optim[seq_num:seq_num+1])
                mutated_counts = np.squeeze(np.exp(mutated_counts.cpu().detach().numpy()))
                mutated_profile = np.squeeze(torch.exp(log_softmax_profile(mutated_profile)).cpu().detach().numpy())

                original_pred = original_counts * original_profile
                mutated_pred = mutated_counts * mutated_profile

                # --- Get seqlets ---
                plot_design_seqlets = seqlets[seqlets["id"] == identifiers[seq_num]] \
                    .rename({"annotation": "motif_name", "attribution": "score"}, axis=1)
                plot_ref_seqlets = seqlets[seqlets["id"] == "ISL1_chr5:51382751-51383709_ref_sequence"] \
                    .rename({"annotation": "motif_name", "attribution": "score"}, axis=1)

                plot_design_seqlets = plot_design_seqlets[["motif_name", "start", "end", "score"]]
                plot_ref_seqlets = plot_ref_seqlets[["motif_name", "start", "end", "score"]]

                # --- Plot ---
                with sns.plotting_context("paper", font_scale=1.5):
                    fig, axes = plt.subplots(
                        nrows=3, ncols=1, figsize=(18, 8), sharex=True,
                        gridspec_kw={'height_ratios': [2, 1, 1]}
                    )

                    # Predictions
                    axes[0].plot(
                        range(0, plot_window),
                        original_pred[plot_target_start:plot_target_end],
                        color='r', linewidth=1, label=f"Original ({round(float(np.log(original_counts)), 2)})"
                    )
                    axes[0].plot(
                        range(0, plot_window),
                        mutated_pred[plot_target_start:plot_target_end],
                        color='b', linewidth=1, label=f"Mutated ({round(float(np.log(mutated_counts)), 2)})"
                    )
                    axes[0].set_ylabel("Predicted Accessibility")
                    axes[0].legend(fontsize=12, loc="upper left")
                    axes[0].set_title(f"{chrom_}:{plot_coordinate_start}-{plot_coordinate_end}")

                    # Attributions baseline
                    plot_logo(
                        X_ref_attr_[0], ax=axes[1],
                        start=plot_seq_start, end=plot_seq_end,
                        annotations=plot_ref_seqlets
                    )
                    axes[1].set_ylabel("Attributions")
                    axes[1].set_ylim(-0.05, 0.15)

                    # Attributions mutated
                    plot_logo(
                        X_optim_attr_[seq_num].detach().cpu(), ax=axes[2],
                        start=plot_seq_start, end=plot_seq_end,
                        annotations=plot_design_seqlets
                    )
                    axes[2].set_ylabel("Attributions")
                    axes[2].set_ylim(-0.05, 0.15)

                    # Mark edits
                    for diff in diffs:
                        pos = diff[0]
                        plot_pos = pos - plot_seq_start
                        axes[1].axvline(plot_pos, color='r', linestyle='--', linewidth=1)
                        axes[2].axvline(plot_pos, color='r', linestyle='--', linewidth=1)

                    #plt.tight_layout()

                    # Save figure
                    out_path = os.path.join(plots_dir, f"{identifiers[seq_num]}_plot.pdf")
                    plt.savefig(out_path, format="pdf")
                    plt.close(fig)

                    print(f"Saved plot → {out_path}")

            count += 1

    print("\nAll designs completed.")

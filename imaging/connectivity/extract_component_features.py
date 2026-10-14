#!/usr/bin/env python3
# Role: Per-participant component features for the Figure 4 chain: for every FWE-significant
#     component of an NBS output folder (imaging/connectivity/nbs_group_analysis.py), the mean Fisher-z FC over the
#     component's edge list (<contrast>_<sign>_comp<k>_suprathresh_fc; missing relmat cells
#     zero-filled as in the NBS engine), its PC1, and the covariate columns merged from
#     --covariates. Writes <nbs-dir>/nbs_suprathreshold_fc.csv.
# Manuscript item: connectivity features of Figure 4 and of the Results on connectivity and
#     outcomes.
"""
extract_component_features.py — Extract per-subject FC for supra-threshold edges
                                    across ALL NBS contrasts (not just main effects).

For each NBS contrast JSON in `nbs-dir`, this script:
1. Reconstructs the contrast vector used by `build_continuous_design`
2. Builds the projection-based reduced model (same as nbs_group_analysis._nbs_continuous)
3. Residualises edges on the reduced model
4. Computes the partial Freedman-Lane t-statistic
5. Thresholds at |t| > t_thresh (default 3.1, same as NBS)
6. Intersects with each significant component's node set to get the exact edges
   that drove the NBS result
7. Extracts per-subject mean FC and PC1 across those supra-threshold edges

Supports: main_apathy, main_quip, interaction_apathyxquip,
          diff_apathy_gt_quip, diff_quip_gt_apathy,
          shared_pos, shared_neg.

Output: nbs_suprathreshold_fc.csv with columns:
  subject_id, group,
  <contrast>_<sign>_comp<N>_suprathresh_fc (mean FC over supra-thresh edges),
  <contrast>_<sign>_comp<N>_n_edges (# supra-thresh edges actually found),
  <contrast>_<sign>_comp<N>_pc1 (subject score on first PC of edges),
  <all covariates merged from covariates_functional.csv>

Usage (on HPC):
    python3 extract_component_features.py
    python3 extract_component_features.py --nbs-dir <path>
    python3 extract_component_features.py --t-thresh 3.1

Notes
-----
• "Mirror" contrasts (diff_quip_gt_apathy = -diff_apathy_gt_quip; shared_neg =
  -shared_pos) are still processed. They identify the SAME physical edges as
  their counterparts and produce identical per-subject FC values under a sign
  flip; columns are emitted for completeness but are informationally redundant
  with the non-mirror contrast.
• The projection-based reduced model matches what _nbs_continuous does in
  nbs_group_analysis.py exactly, so suprathresh edge sets here coincide with the NBS
  edge sets (previous version used a simpler all-nuisance-columns model which
  diverged from NBS t-stats for main effects).
"""
from __future__ import annotations

import argparse
import csv
import json
import os
import sys
from pathlib import Path

import numpy as np
import pandas as pd

HOME = Path(os.environ.get("HOME", str(Path.home())))
GROUP_DIR = Path(os.environ.get("GROUP_DIR",
                                str(HOME / "PPMI/MRI_Data/derivatives/group")))
ANALYSIS_DIR = HOME / "PPMI/MRI_Analysis"
NBS_DIR = GROUP_DIR / "fc_connectome_continuous"
COV_FILE = ANALYSIS_DIR / "covariates_functional.csv"
EPHEMERAL = Path(f"/rds/general/user/{os.environ.get('USER','')}/ephemeral/PPMI/derivatives")


# Canonical list of contrast names that build_continuous_design emits.
# Any JSON in `nbs-dir` matching one of these will be processed; others
# are skipped with a warning.
KNOWN_CONTRASTS = (
    "main_apathy",
    "main_quip",
    "interaction_apathyxquip",
    "diff_apathy_gt_quip",
    "diff_quip_gt_apathy",
    "shared_pos",
    "shared_neg",
)


def _xcpd_subject_dir(subj: str, session: str = "ses-BL") -> Path | None:
    """Return the per-subject xcp-d root that contains a relmat in the
    requested session, preferring ephemeral over home."""
    for base in [EPHEMERAL / f"xcpd_{subj}" / f"sub-{subj}",
                 HOME / "PPMI/MRI_Data/derivatives/xcp_d" / f"sub-{subj}"]:
        if not base.is_dir():
            continue
        sess_dir = base / session
        if sess_dir.is_dir() and list(sess_dir.rglob("*atlas-4S256Parcels*_relmat.tsv")):
            return base
    return None


def load_subject_matrix(subj: str, session: str = "ses-BL") -> np.ndarray | None:
    """Load the per-subject FC matrix for the requested session ONLY.

    Older versions of this function used `sdir.rglob(...)` which crawls every
    session sub-directory and then picks the largest file. That silently mixed
    ses-BL with ses-24M whenever the 24M relmat happened to be larger,
    contaminating baseline analyses. We now restrict the search to the
    requested session directory.
    """
    sdir = _xcpd_subject_dir(subj, session=session)
    if sdir is None:
        return None
    sess_dir = sdir / session
    relmat_files = sorted(sess_dir.rglob("*atlas-4S256Parcels*_relmat.tsv"))
    if not relmat_files:
        return None
    # If multiple runs exist within the same session, pick the largest
    best = max(relmat_files, key=lambda p: p.stat().st_size)
    try:
        df = pd.read_csv(best, sep="\t", index_col=0)
        return df.values.astype(np.float64)
    except Exception:
        return None


def build_contrast_vec(name: str, ev_names: list[str]) -> np.ndarray | None:
    """Reconstruct the contrast vector for a named NBS contrast.

    Must match build_continuous_design in nbs_group_analysis.py.
    Returns None if a required EV column is missing.
    """
    def _idx(ev): return ev_names.index(ev) if ev in ev_names else -1
    n_ev = len(ev_names)
    idx_ap = _idx("apathy_composite_c")
    idx_qu = _idx("quip_c")
    idx_int = _idx("apathy_x_quip")
    if idx_ap < 0 or idx_qu < 0:
        return None
    cvec = np.zeros(n_ev)
    if name == "main_apathy":
        cvec[idx_ap] = 1.0
    elif name == "main_quip":
        cvec[idx_qu] = 1.0
    elif name == "interaction_apathyxquip":
        if idx_int < 0:
            return None
        cvec[idx_int] = 1.0
    elif name == "diff_apathy_gt_quip":
        cvec[idx_ap] = 1.0
        cvec[idx_qu] = -1.0
    elif name == "diff_quip_gt_apathy":
        cvec[idx_ap] = -1.0
        cvec[idx_qu] = 1.0
    elif name == "shared_pos":
        cvec[idx_ap] = 1.0
        cvec[idx_qu] = 1.0
    elif name == "shared_neg":
        cvec[idx_ap] = -1.0
        cvec[idx_qu] = -1.0
    else:
        return None
    return cvec


def _residualise(y: np.ndarray, X_nui: np.ndarray) -> np.ndarray:
    if X_nui.shape[1] == 0:
        return y - y.mean()
    X_int = np.column_stack([np.ones(len(y)), X_nui])
    beta, *_ = np.linalg.lstsq(X_int, y, rcond=None)
    return y - X_int @ beta


def compute_contrast_t_stats(design: np.ndarray, edges: np.ndarray,
                              cvec: np.ndarray) -> np.ndarray:
    """Compute the Freedman-Lane partial t-statistic for each edge under
    contrast `cvec`, using the SAME projection-based reduced model as
    _nbs_continuous in nbs_group_analysis.py.
    """
    n_subj, n_ev = design.shape
    n_edges = edges.shape[1]
    # Reduced model: project the contrast direction out of the full design
    cvec_norm = cvec / np.sqrt((cvec ** 2).sum())
    X_reduced = design - np.outer(design @ cvec_norm, cvec_norm)
    col_var = X_reduced.var(axis=0)
    keep_cols = np.where(col_var > 1e-12)[0]
    X_nui = X_reduced[:, keep_cols]

    # Residualise the effective regressor and each edge on the reduced model
    x = design @ cvec
    x_res = _residualise(x, X_nui)
    edges_res = np.empty_like(edges)
    for e in range(n_edges):
        edges_res[:, e] = _residualise(edges[:, e], X_nui)

    # Partial t-statistic of edge ~ x_res
    x_c = x_res - x_res.mean()
    ss_x = (x_c ** 2).sum()
    if ss_x < 1e-12:
        return np.zeros(n_edges)
    beta = (edges_res * x_c[:, None]).sum(axis=0) / ss_x
    y_hat = beta[None, :] * x_c[:, None]
    resid = edges_res - y_hat
    df = n_subj - max(1, X_nui.shape[1]) - 2
    if df < 1:
        return np.zeros(n_edges)
    sigma2 = (resid ** 2).sum(axis=0) / df
    se = np.sqrt(sigma2 / ss_x + 1e-12)
    return beta / se


def main():
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--t-thresh", type=float, default=3.1,
                   help="Primary edge t-threshold (must match NBS; default 3.1)")
    p.add_argument("--nbs-dir", default="",
                   help="NBS output directory; default: $GROUP_DIR/fc_connectome_continuous")
    p.add_argument("--out", default="",
                   help="Output CSV path; default: <nbs-dir>/nbs_suprathreshold_fc.csv")
    p.add_argument("--contrasts", nargs="*", default=None,
                   help="Optional subset of contrast names to process "
                        "(default: all NBS JSONs found). Use to skip mirror "
                        "contrasts, e.g. '--contrasts main_apathy main_quip "
                        "interaction_apathyxquip diff_apathy_gt_quip shared_pos'.")
    p.add_argument("--covariates", default=str(COV_FILE),
                   help="Covariates CSV to merge into the output")
    p.add_argument("--structural-csv",
                   default=str(HOME / "PPMI/MRI_Data/derivatives/all_structural_measures_ses-BL.csv"),
                   help="Optional all_structural_measures_ses-BL.csv to merge "
                        "FreeSurfer-derived atrophy covariates "
                        "(brain_atrophy_index, fs_eTIV, fs_TotalGray, "
                        "mean_cortical_thickness)")
    p.add_argument("--session", default="ses-BL",
                   help="Which session subdir of xcp-d to read from "
                        "(default ses-BL). CRITICAL: must match the NBS run's "
                        "session — mixing BL and 24M relmats silently "
                        "contaminates results.")
    args = p.parse_args()

    nbs_dir = Path(args.nbs_dir) if args.nbs_dir else NBS_DIR
    t_thresh = args.t_thresh

    # ── Load subject list and design matrix ──
    subj_file = nbs_dir / "design_subjects.tsv"
    design_file = nbs_dir / "design_matrix.tsv"
    for f in (subj_file, design_file):
        if not f.is_file():
            print(f"ERROR: {f} not found", file=sys.stderr)
            sys.exit(1)

    subj_df = pd.read_csv(subj_file, sep="\t")
    subj_df["sid"] = subj_df["subject_id"].str.replace("sub-", "")
    subjects = list(zip(subj_df["sid"], subj_df["group"]))
    n_subj = len(subjects)
    print(f"Subjects: {n_subj}")

    design_df = pd.read_csv(design_file, sep="\t")
    ev_names = list(design_df.columns)
    design = design_df.values
    print(f"Design:   {design.shape[0]} × {design.shape[1]} ({len(ev_names)} EVs)")

    # ── Discover contrasts from NBS JSONs ──
    nbs_jsons = sorted(nbs_dir.glob("nbs_*.json"))
    found_contrasts: list[str] = []
    for f in nbs_jsons:
        name = f.stem.replace("nbs_", "")
        if name in KNOWN_CONTRASTS:
            found_contrasts.append(name)
        else:
            print(f"  skipping unknown contrast: {name}")
    if args.contrasts:
        found_contrasts = [c for c in found_contrasts if c in args.contrasts]
    if not found_contrasts:
        print("ERROR: no known NBS contrasts to process", file=sys.stderr)
        sys.exit(1)
    print(f"Contrasts: {found_contrasts}")

    # ── Load all connectivity matrices ──
    print(f"Loading connectivity matrices (session: {args.session})...")
    all_mats = []
    for i, (sid, grp) in enumerate(subjects):
        if (i + 1) % 50 == 0:
            print(f"  {i+1}/{n_subj}...")
        all_mats.append(load_subject_matrix(sid, session=args.session))

    ref_shape = next((m.shape for m in all_mats if m is not None), None)
    if ref_shape is None:
        print("ERROR: no matrices loaded", file=sys.stderr)
        sys.exit(1)
    n_rois = ref_shape[0]

    mats_stacked = np.zeros((n_subj, n_rois, n_rois))
    for i, m in enumerate(all_mats):
        if m is not None and m.shape == ref_shape:
            mats_stacked[i] = m

    # Fisher-z, zero diagonal, nan→0 (same as NBS)
    mats_z = np.arctanh(np.clip(mats_stacked, -0.999, 0.999))
    for i in range(n_subj):
        np.fill_diagonal(mats_z[i], 0)
    mats_z = np.nan_to_num(mats_z, nan=0.0, posinf=0.0, neginf=0.0)

    triu = np.triu_indices(n_rois, k=1)
    n_edges = len(triu[0])
    edges = mats_z[:, triu[0], triu[1]]
    print(f"  {n_subj} × {n_rois} rois → {n_edges} upper-triangle edges")

    # ── Process each contrast ──
    component_edges: dict[str, dict] = {}

    for contrast_name in found_contrasts:
        cvec = build_contrast_vec(contrast_name, ev_names)
        if cvec is None:
            print(f"  {contrast_name}: required EV columns missing — skipping")
            continue

        t_obs = compute_contrast_t_stats(design, edges, cvec)
        n_pos = int((t_obs > t_thresh).sum())
        n_neg = int((t_obs < -t_thresh).sum())
        print(f"\n  {contrast_name}: {n_pos} pos edges, {n_neg} neg edges "
              f"supra-threshold at |t|>{t_thresh}")

        json_file = nbs_dir / f"nbs_{contrast_name}.json"
        if not json_file.is_file():
            print(f"    WARNING: {json_file} missing — skipping")
            continue
        nbs_data = json.load(open(json_file))

        # Build a fast (i,j) -> upper-triangle linear index map for this n_rois
        # so we can convert each JSON edge pair to an edge-vector index.
        ij_to_eidx = {(int(triu[0][e]), int(triu[1][e])): e
                      for e in range(n_edges)}

        for sign, key in [("pos", "obs_components_pos"),
                          ("neg", "obs_components_neg")]:
            sig_comps = [c for c in nbs_data.get(key, []) if c["p"] < 0.05]
            for ci, comp in enumerate(sig_comps):
                nodes = set(comp["nodes"])
                n_edges_nbs = comp["n_edges"]
                # Use the JSON's explicit edge list — these are the canonical
                # NBS suprathreshold edges that defined the component. (The
                # previous strategy of re-computing t-stats and re-thresholding
                # produced a non-identical subset because the script's
                # reduced-model t-stats can diverge by ~1e-3 from the NBS run.)
                comp_edge_idx = []
                json_edges = comp.get("edges", [])
                for e_ij in json_edges:
                    a, b = int(e_ij[0]), int(e_ij[1])
                    if a > b:
                        a, b = b, a
                    idx = ij_to_eidx.get((a, b))
                    if idx is not None:
                        comp_edge_idx.append(idx)
                label = f"{contrast_name}_{sign}_comp{ci+1}"
                component_edges[label] = {
                    "edge_indices": comp_edge_idx,
                    "n_edges_nbs": n_edges_nbs,
                    "n_edges_suprathresh": len(comp_edge_idx),
                    "t_values": t_obs[comp_edge_idx].tolist() if comp_edge_idx else [],
                    "nodes": sorted(nodes),
                    "p": comp["p"],
                }
                if comp_edge_idx:
                    mean_abs_t = float(np.abs(t_obs[comp_edge_idx]).mean())
                    match_flag = "" if len(comp_edge_idx) == n_edges_nbs else \
                                 f"  [WARN: edge-list/{n_edges_nbs} mismatch]"
                    print(f"    {label}: {len(comp_edge_idx)} edges from JSON "
                          f"(NBS reported {n_edges_nbs}), mean |t|={mean_abs_t:.2f}, "
                          f"p={comp['p']:.4f}{match_flag}")
                else:
                    print(f"    {label}: 0 edges loaded from JSON (NBS reported "
                          f"{n_edges_nbs}) — missing 'edges' field?")

    if not component_edges:
        print("\nNo significant components found across any contrast.", file=sys.stderr)
        sys.exit(0)

    # ── Build whole-FC specificity-control mask ──
    # Union of all apathy-circuit suprathreshold edges (pos + neg components of
    # main_apathy). whole_FC_mean_excl_apathy = per-subject mean FC over every
    # upper-triangle edge NOT in this union — the FC analogue of
    # whole_connectome_mean_excl_apathy in the SC pipeline. Used as a
    # circuit-specificity covariate in downstream cognition models.
    apathy_edge_indices = set()
    for label, info in component_edges.items():
        if label.startswith("main_apathy_"):
            apathy_edge_indices.update(info["edge_indices"])
    print(f"\nUnion of main_apathy suprathresh edges: {len(apathy_edge_indices)} "
          f"(out of {n_edges} upper-triangle edges)")
    non_apathy_mask = np.ones(n_edges, dtype=bool)
    for ei in apathy_edge_indices:
        non_apathy_mask[ei] = False

    # ── Extract per-subject FC for supra-threshold edges ──
    print("\nExtracting per-subject FC + PC1...")
    results = []
    for i, (sid, grp) in enumerate(subjects):
        row = {"subject_id": sid, "group": grp}
        # Whole-FC specificity covariate (mean across non-apathy edges)
        row["whole_FC_mean_excl_apathy"] = float(
            np.mean(mats_z[i, triu[0][non_apathy_mask], triu[1][non_apathy_mask]])
        )
        # Also a global all-edges mean, parallel to the SC variable
        row["whole_FC_mean"] = float(
            np.mean(mats_z[i, triu[0], triu[1]])
        )
        for label, info in component_edges.items():
            edge_idx = info["edge_indices"]
            if not edge_idx:
                row[f"{label}_suprathresh_fc"] = np.nan
                row[f"{label}_n_edges"] = 0
                continue
            vals = mats_z[i, triu[0][edge_idx], triu[1][edge_idx]]
            row[f"{label}_suprathresh_fc"] = float(np.mean(vals))
            row[f"{label}_n_edges"] = len(edge_idx)
        results.append(row)

    df = pd.DataFrame(results)

    # PC1 per component (sign-flipped so that PC1 > 0 ⇔ higher mean FC)
    for label, info in component_edges.items():
        edge_idx = info["edge_indices"]
        if len(edge_idx) < 2:
            continue
        edge_data = mats_z[:, triu[0][edge_idx], triu[1][edge_idx]]
        edge_std = (edge_data - edge_data.mean(axis=0)) / (edge_data.std(axis=0) + 1e-12)
        U, S, Vt = np.linalg.svd(edge_std, full_matrices=False)
        pc1 = U[:, 0] * S[0]
        if np.corrcoef(pc1, edge_data.mean(axis=1))[0, 1] < 0:
            pc1 = -pc1
        df[f"{label}_pc1"] = pc1
        print(f"  {label} PC1: {(S[0]**2)/(S**2).sum()*100:.1f}% variance explained")

    # Merge covariates
    cov = pd.read_csv(args.covariates)
    cov["subject_id"] = cov["subject_id"].astype(str).str.replace("sub-", "")
    df["subject_id"] = df["subject_id"].astype(str)
    df = df.merge(cov, on="subject_id", how="left", suffixes=("", "_cov"))
    if "group_cov" in df.columns:
        df.drop(columns=["group_cov"], inplace=True)

    # Merge FreeSurfer-derived atrophy covariates from all_structural_measures
    struct_path = Path(args.structural_csv)
    if struct_path.is_file():
        struct = pd.read_csv(struct_path)
        struct["PATNO"] = struct["PATNO"].astype(str).str.replace("sub-", "")
        # Derive brain atrophy index = BrainSegNotVent / eTIV
        if {"fs_BrainSegNotVent", "fs_EstimatedTotalIntraCranialVol"}.issubset(struct.columns):
            struct["brain_atrophy_index"] = (
                struct["fs_BrainSegNotVent"]
                / struct["fs_EstimatedTotalIntraCranialVol"]
            )
        # Derive mean cortical thickness across lh+rh DK regions
        thk_cols = [c for c in struct.columns
                    if c.startswith(("fs_lh_", "fs_rh_"))
                    and c.endswith("_thickness")
                    and not c.endswith("_thickstd")]
        if thk_cols:
            struct["mean_cortical_thickness"] = struct[thk_cols].mean(axis=1)
        keep = ["PATNO", "brain_atrophy_index", "mean_cortical_thickness",
                "fs_EstimatedTotalIntraCranialVol", "fs_BrainSegNotVent",
                "fs_TotalGray", "fs_SubCortGray", "fs_Cortex"]
        keep = [c for c in keep if c in struct.columns]
        struct_keep = struct[keep].rename(
            columns={"PATNO": "subject_id",
                     "fs_EstimatedTotalIntraCranialVol": "fs_eTIV"})
        struct_keep["subject_id"] = struct_keep["subject_id"].astype(str)
        df = df.merge(struct_keep, on="subject_id", how="left",
                      suffixes=("", "_struct"))
        n_matched = df["brain_atrophy_index"].notna().sum() if "brain_atrophy_index" in df.columns else 0
        print(f"  Merged structural atrophy covariates: "
              f"{n_matched}/{len(df)} subjects matched")
    else:
        print(f"  WARN: structural CSV not found at {struct_path} — "
              f"skipping brain_atrophy_index merge")

    out_path = Path(args.out) if args.out else nbs_dir / "nbs_suprathreshold_fc.csv"
    df.to_csv(out_path, index=False)
    print(f"\nSaved: {out_path}")
    print(f"  {len(df)} subjects, {len(df.columns)} columns, "
          f"{len(component_edges)} components across {len(found_contrasts)} contrasts")

    # Quick correlation sanity checks
    for predictor in ("apathy_composite", "quip"):
        if predictor not in df.columns:
            continue
        pred = pd.to_numeric(df[predictor], errors="coerce")
        print(f"\nCorrelations with {predictor}:")
        for label in component_edges:
            for suffix in ("_suprathresh_fc", "_pc1"):
                col = f"{label}{suffix}"
                if col not in df.columns:
                    continue
                vals = pd.to_numeric(df[col], errors="coerce")
                mask = pred.notna() & vals.notna()
                if mask.sum() <= 10:
                    continue
                from scipy import stats
                r, p = stats.pearsonr(pred[mask], vals[mask])
                sig = " *" if p < 0.05 else ""
                print(f"  {col:55s} r={r:+.3f}  p={p:.4f}{sig}")

    # Edge info JSON (for the mediation / figure scripts)
    edge_info = {}
    for label, info in component_edges.items():
        tvals = np.array(info["t_values"]) if info["t_values"] else np.array([0.0])
        edge_info[label] = {
            "n_edges_nbs": info["n_edges_nbs"],
            "n_edges_suprathresh": info["n_edges_suprathresh"],
            "mean_abs_t": float(np.abs(tvals).mean()),
            "max_abs_t": float(np.abs(tvals).max()),
            "p_nbs": info["p"],
        }
    with open(nbs_dir / "nbs_suprathreshold_edge_info.json", "w") as f:
        json.dump(edge_info, f, indent=2)
    print(f"  Edge info: {nbs_dir / 'nbs_suprathreshold_edge_info.json'}")


if __name__ == "__main__":
    main()

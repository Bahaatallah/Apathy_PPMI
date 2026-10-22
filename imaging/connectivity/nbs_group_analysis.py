#!/usr/bin/env python3
# Role: NBS engine (connectome_continuous analysis): edge-wise GLM, residualised-regressor
#     permutation, component FWE, 10,000 permutations, |t| > 3.1 by default. Reads the largest
#     4S256Parcels relmat per participant.
# Manuscript item: Table 2, Figure 3 and the NBS sensitivity analyses of the Supplement.
"""
Dimensional Network-Based Statistic (NBS) analysis of resting-state functional
connectivity from XCP-D connectivity matrices.

Input
  * One connectivity matrix per participant for the selected atlas
    (XCPD_ATLAS or --atlas) and session (SESSION, default ses-BL), looked up
    in DERIV_DIR/xcp_d and then in --fmriprep-root. When several runs exist
    the largest relmat file is used. Matrices are Fisher-z transformed, with
    the diagonal and non-finite entries set to zero.
  * Participants are read from --subjects-file (one subject|group pair per
    line) or, without it, discovered from the BIDS group folders.

Design
  * Edge-wise general linear model. The regressors of interest are the
    mean-centred apathy composite, the mean-centred QUIP score and their
    product. The nuisance covariates are the depression factor, age, sex,
    site, TIV (z-scored) and mean framewise displacement, read from
    --covariates. Participants missing any required covariate are excluded.
  * Contrasts are the main effects of apathy and QUIP, the apathy-by-QUIP
    interaction, the differential contrasts (apathy > QUIP, QUIP > apathy)
    and the shared contrasts (both positive, both negative).

Inference
  * For each contrast, the edges and the contrast regressor are residualised
    on the reduced model (the design with the contrast direction projected
    out), and the residualised regressor is permuted with a fixed seed.
  * Edges with t above --t-thresh (default 3.1), or below its negative, form
    connected components whose size is their number of edges. Each component
    is compared with the permutation distribution of the largest component of
    the same sign, which gives family-wise error corrected p-values
    (--nperm permutations, default 10,000).

Sensitivity options
  --extended-covs LIST   adds any of disease_duration_months, updrs_iii_off,
                         ledd and moca as nuisance covariates
  --include-extended     adds all four
  --quip-binary          replaces the QUIP score by its presence (QUIP > 0)

Outputs, in GROUP_DIR/fc_connectome_continuous
  nbs_<contrast>.json, design_matrix.tsv, design_subjects.tsv,
  design_diagnostics.json and fc_connectome_continuous_summary.json

Usage
  python3 nbs_group_analysis.py --analysis connectome_continuous
      --nperm 10000 --t-thresh 3.1 --subjects-file <roster>
      --covariates <table> --fmriprep-root <path>
      [--extended-covs <list> | --include-extended] [--quip-binary]

Requires numpy, and pandas to read relmat files with row and column labels.
"""

from __future__ import annotations

import argparse
import csv
import json
import os
import sys
from pathlib import Path

import numpy as np

# ── Configuration ────────────────────────────────────────────────────
RDS_HOME = os.environ.get("RDS_HOME",
                           f"/rds/general/user/{os.environ.get('USER', 'user')}/home")
RDS_EPHEMERAL = os.environ.get("RDS_EPHEMERAL",
                                f"/rds/general/user/{os.environ.get('USER', 'user')}/ephemeral")
DERIV_DIR = Path(os.environ.get("DERIV_DIR", f"{RDS_HOME}/PPMI/MRI_Data/derivatives"))
GROUP_DIR = Path(os.environ.get("GROUP_DIR", str(DERIV_DIR / "group")))
BIDS_ROOT = Path(os.environ.get("BIDS_ROOT", f"{RDS_HOME}/PPMI/MRI_Data/BIDS"))
SCRIPTS_DIR = Path(os.path.dirname(os.path.abspath(__file__)))

# Root under which the XCP-D outputs are searched when they are not under DERIV_DIR
# (xcpd_<id>/sub-<id>/<session> or xcp_d/sub-<id>/<session>); overridable with --fmriprep-root.
FMRIPREP_ROOT = Path(os.environ.get(
    "FMRIPREP_ROOT",
    os.environ.get("EPHEMERAL_DERIV", f"{RDS_EPHEMERAL}/PPMI/derivatives")))

GROUPS = ["Apathy_BL", "Both_BL", "ICB_BL", "Neither_BL"]
# Session timepoint. Default BL; for 24M analysis set SESSION=ses-24M in env.
SESSION = os.environ.get("SESSION", "ses-BL")
XCPD_ATLAS = os.environ.get("XCPD_ATLAS", "4S256Parcels")


def _xcpd_subject_dir(subj: str) -> Path:
    """Return the XCP-D directory for the CURRENT SESSION that contains a
    denoised BOLD.

    Session-aware: candidate paths include `ses-{SESSION}` so that for a
    subject processed at both BL and 24M, we only return the session
    requested by the SESSION env var. When the requested session is
    missing for a subject, returns a non-existent path so downstream
    rglob() yields nothing and the subject is correctly skipped.

    Critical: the 05b pipeline rsyncs XCP-D outputs HOME after the job
    finishes. When HOME quota is exceeded, rsync leaves behind an empty
    sub-<id>/ stub on HOME while the real data lives on EPHEMERAL. We
    must prefer whichever dir actually has data, not the first that
    exists.

    Search order (first with *desc-denoised*bold.nii.gz wins):
      1. $DERIV_DIR/xcp_d/sub-<subj>/<SESSION>          (home, primary)
      2. $FMRIPREP_ROOT/xcpd_<subj>/sub-<subj>/<SESSION> (05b on ephemeral)
      3. $FMRIPREP_ROOT/xcp_d/sub-<subj>/<SESSION>      (legacy layout)
    """
    candidates = [
        DERIV_DIR / "xcp_d" / f"sub-{subj}" / SESSION,
        FMRIPREP_ROOT / f"xcpd_{subj}" / f"sub-{subj}" / SESSION,
        FMRIPREP_ROOT / "xcp_d" / f"sub-{subj}" / SESSION,
    ]
    for c in candidates:
        if c.is_dir():
            try:
                if next(c.rglob("*desc-denoised*bold.nii.gz"), None) is not None:
                    return c
            except OSError:
                continue
    # No populated dir found for this SESSION — return first candidate so
    # downstream rglob() returns empty and the subject is skipped, NOT
    # silently substituted with another session's data.
    return candidates[0]

# Required covariate columns (group analyses fail without these)


# ======================================================================
# Subject discovery
# ======================================================================
def find_xcp_subjects() -> list[tuple[str, str]]:
    """Find subjects with XCP-D denoised rest BOLD + connectome for the chosen atlas."""
    subjects = []
    for grp in GROUPS:
        grp_dir = BIDS_ROOT / grp
        if not grp_dir.is_dir():
            continue
        for sd in sorted(grp_dir.glob("sub-*")):
            subj = sd.name.replace("sub-", "")
            xcp_dir = _xcpd_subject_dir(subj)
            if not xcp_dir.is_dir():
                continue
            has_rest = bool(list(xcp_dir.rglob("*task-rest*desc-denoised*bold.nii.gz")))
            if has_rest:
                subjects.append((subj, grp))
    return subjects


# ======================================================================
# Continuous (dimensional) design matrix
# ======================================================================
# Uses raw clinical scores as continuous regressors rather than dichotomised
# group labels. More power + tests apathy × ICB interaction directly.

CONTINUOUS_PRIMARY = ["apathy_composite", "quip"]          # main effects


def _read_cov_csv(cov_file: Path) -> tuple[list[str], dict[str, dict]]:
    """Minimal CSV reader that returns (fieldnames, dict keyed by subject_id)."""
    rows: dict[str, dict] = {}
    fieldnames: list[str] = []
    with open(cov_file) as f:
        reader = csv.DictReader(f)
        fieldnames = list(reader.fieldnames or [])
        for row in reader:
            sid = (row.get("subject_id")
                   or row.get("participant_id", "")).replace("sub-", "").strip()
            if sid:
                rows[sid] = row
    return fieldnames, rows


EXTENDED_COVS_VALID = ["disease_duration_months", "updrs_iii_off", "ledd", "moca"]


def build_continuous_design(cov_file: Path,
                             subjects: list[tuple[str, str]],
                             include_depression: bool = True,
                             extended_covs: list[str] | None = None,
                             quip_binary: bool = False):
    """Build a dimensional design matrix for apathy × ICB analysis.

    Returns (kept, design, ev_names, contrasts, diag) where:
        kept        list[(subject_id, group)]  subjects surviving completeness
        design      (n × n_ev) np.ndarray  mean-centred, with sex/site dummies
        ev_names    list[str]  one name per column
        contrasts   list[(name, np.ndarray)]  primary contrasts + F spec
        diag        dict  diagnostic info (VIFs, correlations, dropped subjects)

    EV order:
        [0] apathy_composite_c     (mean-centred)
        [1] quip_c                 (mean-centred)
        [2] apathy_x_quip          (product of demeaned apathy and quip)
        [3] depression_factor_c    (optional; primary nuisance)
        [4] age_c
        [5..] any extended nuisance covariates present in `extended_covs`
              (disease_duration_c, updrs_iii_off_c, ledd_c, moca_c — in that
              canonical order, regardless of the order supplied)
        [..] sex dummy, site dummies
        [n-2] TIV_c
        [n-1] mean_FD_c

    `extended_covs` is any subset of EXTENDED_COVS_VALID. Pass None/[] for the
    lean model. Subjects missing ANY required column are dropped.
    """
    extended_covs = list(extended_covs or [])
    invalid = [c for c in extended_covs if c not in EXTENDED_COVS_VALID]
    if invalid:
        raise ValueError(
            f"unknown extended covariates: {invalid}. "
            f"valid names: {EXTENDED_COVS_VALID}")
    # Normalise to canonical order so EV naming is stable
    extended_covs = [c for c in EXTENDED_COVS_VALID if c in extended_covs]

    fieldnames, rows = _read_cov_csv(cov_file)

    required = list(CONTINUOUS_PRIMARY) + ["age", "sex", "site", "TIV", "mean_FD"]
    if include_depression:
        required.append("depression_factor")
    required.extend(extended_covs)

    missing_cols = [c for c in required if c not in fieldnames]
    if missing_cols:
        raise RuntimeError(
            f"covariates file missing columns required for continuous "
            f"analysis: {missing_cols}\n  header: {fieldnames}")

    # ── Filter subjects ──────────────────────────────────────────────
    # Recognise the common missing-value tokens pandas / numpy / R emit when
    # dumping to CSV. Previously only empty and "NA" were caught, which let
    # literal "NaN"/"nan" values (what pandas writes for float NaN by default)
    # slip through and poison lstsq with DLASCL errors downstream.
    _MISSING_TOKENS = {"", "na", "nan", "null", "none", "n/a", "#n/a"}
    kept: list[tuple[str, str]] = []
    dropped: list[tuple[str, str]] = []
    for sid, grp in subjects:
        if sid not in rows:
            dropped.append((sid, "not_in_covariates"))
            continue
        row = rows[sid]
        missing = []
        for c in required:
            v = (row.get(c) or "").strip()
            if v.lower() in _MISSING_TOKENS:
                missing.append(c)
        if missing:
            dropped.append((sid, f"missing:{','.join(missing)}"))
            continue
        kept.append((sid, grp))

    if len(kept) < 10:
        raise RuntimeError(
            f"only {len(kept)} subjects have all continuous regressors "
            f"populated (need >=10). First few dropped:\n  " +
            "\n  ".join(f"sub-{s}: {r}" for s, r in dropped[:15]))

    n = len(kept)

    def _col(name: str) -> np.ndarray:
        return np.array([float(rows[sid][name]) for sid, _ in kept],
                        dtype=np.float64)

    # ── Continuous regressors (raw, then demean) ─────────────────────
    apathy = _col("apathy_composite")
    quip = _col("quip")
    if quip_binary:
        # Sensitivity: dichotomize zero-inflated QUIP (see PPMI cohort —
        # 77% of subjects have QUIP=0, only 70/310 carry any signal, and
        # within the non-zero subset SD≈0.6 over range 1–4).
        n_nonzero = int((quip > 0).sum())
        quip = (quip > 0).astype(np.float64)
        print(f"  [quip_binary] QUIP dichotomised to presence/absence: "
              f"{n_nonzero}/{len(quip)} subjects = 1 (QUIP>0), "
              f"{len(quip)-n_nonzero} = 0",
              file=sys.stderr)
    age = _col("age")
    tiv = _col("TIV")
    mfd = _col("mean_FD")
    dep = _col("depression_factor") if include_depression else None
    disease_dur = _col("disease_duration_months") if "disease_duration_months" in extended_covs else None
    updrs = _col("updrs_iii_off") if "updrs_iii_off" in extended_covs else None
    ledd = _col("ledd") if "ledd" in extended_covs else None
    moca = _col("moca") if "moca" in extended_covs else None

    apathy_c = apathy - apathy.mean()
    quip_c = quip - quip.mean()
    age_c = age - age.mean()
    # Z-score TIV to avoid ill-conditioned design (raw TIV has std~200k)
    tiv_c = (tiv - tiv.mean()) / (tiv.std() + 1e-12)
    mfd_c = mfd - mfd.mean()
    dep_c = (dep - dep.mean()) if dep is not None else None
    disease_dur_c = (disease_dur - disease_dur.mean()) if disease_dur is not None else None
    updrs_c = (updrs - updrs.mean()) if updrs is not None else None
    ledd_c = (ledd - ledd.mean()) if ledd is not None else None
    moca_c = (moca - moca.mean()) if moca is not None else None

    # Interaction = product of centred variables (standard practice, reduces
    # collinearity between interaction and main effects)
    apathy_x_quip = apathy_c * quip_c

    # ── Categorical dummies (sex, site) ──────────────────────────────
    def _dummies(name: str, drop_first: bool = True) -> tuple[list[str], np.ndarray]:
        raw = [(rows[sid].get(name) or "").strip() for sid, _ in kept]
        levels = sorted(set(v for v in raw if v))
        if len(levels) <= 1:
            return [], np.zeros((n, 0))
        if drop_first:
            ref = levels[0]
            keep_levels = levels[1:]
        else:
            keep_levels = levels
            ref = None
        cols = []
        labels = []
        for lv in keep_levels:
            c = np.array([1.0 if r == lv else 0.0 for r in raw])
            c = c - c.mean()   # demean so contrast weights of 0 are true
            cols.append(c)
            labels.append(f"{name}_{lv}" + (f"_vs_{ref}" if ref else ""))
        return labels, np.column_stack(cols) if cols else np.zeros((n, 0))

    sex_labels, sex_mat = _dummies("sex")
    site_labels, site_mat = _dummies("site")

    # ── Assemble design matrix ───────────────────────────────────────
    ev_names: list[str] = []
    ev_cols: list[np.ndarray] = []

    ev_names.append("apathy_composite_c");   ev_cols.append(apathy_c)
    ev_names.append("quip_c");               ev_cols.append(quip_c)
    ev_names.append("apathy_x_quip");        ev_cols.append(apathy_x_quip)
    if dep_c is not None:
        ev_names.append("depression_factor_c"); ev_cols.append(dep_c)
    ev_names.append("age_c");                ev_cols.append(age_c)
    if disease_dur_c is not None:
        ev_names.append("disease_duration_c");   ev_cols.append(disease_dur_c)
    if updrs_c is not None:
        ev_names.append("updrs_iii_off_c");      ev_cols.append(updrs_c)
    if ledd_c is not None:
        ev_names.append("ledd_c");               ev_cols.append(ledd_c)
    if moca_c is not None:
        ev_names.append("moca_c");               ev_cols.append(moca_c)
    ev_names.extend(sex_labels)
    for c in sex_mat.T:
        ev_cols.append(c)
    ev_names.extend(site_labels)
    for c in site_mat.T:
        ev_cols.append(c)
    ev_names.append("TIV_c");                ev_cols.append(tiv_c)
    ev_names.append("mean_FD_c");            ev_cols.append(mfd_c)

    design = np.column_stack(ev_cols)
    n_ev = design.shape[1]

    # ── Sanity: no NaN/Inf allowed in the design (cheap, early, explicit) ─
    # A single NaN in any column detonates np.linalg.lstsq / SVD with the
    # utterly opaque "DLASCL parameter 4 illegal value" LAPACK message. Fail
    # loudly with a column-level report instead.
    if not np.all(np.isfinite(design)):
        bad_cols = [ev_names[i] for i in range(n_ev)
                    if not np.all(np.isfinite(design[:, i]))]
        n_bad_rows = int((~np.isfinite(design).all(axis=1)).sum())
        raise RuntimeError(
            f"design matrix contains NaN/Inf in columns: {bad_cols} "
            f"({n_bad_rows} subjects affected). "
            "Upstream filter likely failed to recognise a missing-value "
            "token in the covariates file.")

    # ── Diagnostics: correlations + VIF for primary regressors ───────
    primary_idx = [ev_names.index(n) for n in
                   ["apathy_composite_c", "quip_c", "apathy_x_quip"]]
    primary_corr = np.corrcoef(design[:, primary_idx].T)

    def _vif(X: np.ndarray, col: int) -> float:
        others = np.delete(X, col, axis=1)
        if others.shape[1] == 0:
            return 1.0
        # R² of col regressed on the others
        y = X[:, col]
        X_int = np.column_stack([np.ones(len(y)), others])
        beta, *_ = np.linalg.lstsq(X_int, y, rcond=None)
        y_hat = X_int @ beta
        ss_res = ((y - y_hat) ** 2).sum()
        ss_tot = ((y - y.mean()) ** 2).sum()
        if ss_tot < 1e-12:
            return float("inf")
        r2 = 1.0 - ss_res / ss_tot
        if r2 >= 0.9999:
            return float("inf")
        return 1.0 / (1.0 - r2)

    vifs = {ev_names[i]: _vif(design, i) for i in primary_idx}

    diag = {
        "n_kept": n,
        "n_dropped": len(dropped),
        "dropped": dropped,
        "primary_correlation": primary_corr.tolist(),
        "primary_ev_order": [ev_names[i] for i in primary_idx],
        "vifs_primary": vifs,
        "score_ranges": {
            "apathy_composite": [float(apathy.min()), float(apathy.max()),
                                  float(apathy.mean()), float(apathy.std())],
            "quip": [float(quip.min()), float(quip.max()),
                     float(quip.mean()), float(quip.std())],
        },
    }

    # ── Contrasts: main effects + interaction + differential + shared ─
    contrasts: list[tuple[str, np.ndarray]] = []
    for name in ("apathy_composite_c", "quip_c", "apathy_x_quip"):
        c = np.zeros(n_ev)
        c[ev_names.index(name)] = 1.0
        # Strip the "_c" / keep friendly name
        pretty = {"apathy_composite_c": "main_apathy",
                  "quip_c": "main_quip",
                  "apathy_x_quip": "interaction_apathyxquip"}[name]
        contrasts.append((pretty, c))

    # Differential contrasts: apathy vs QUIP (tests specificity)
    idx_ap = ev_names.index("apathy_composite_c")
    idx_qu = ev_names.index("quip_c")

    c_ap_gt_qu = np.zeros(n_ev)
    c_ap_gt_qu[idx_ap] = 1.0
    c_ap_gt_qu[idx_qu] = -1.0
    contrasts.append(("diff_apathy_gt_quip", c_ap_gt_qu))

    c_qu_gt_ap = np.zeros(n_ev)
    c_qu_gt_ap[idx_qu] = 1.0
    c_qu_gt_ap[idx_ap] = -1.0
    contrasts.append(("diff_quip_gt_apathy", c_qu_gt_ap))

    # Shared contrasts: common motivational dysregulation
    c_shared_pos = np.zeros(n_ev)
    c_shared_pos[idx_ap] = 1.0
    c_shared_pos[idx_qu] = 1.0
    contrasts.append(("shared_pos", c_shared_pos))

    c_shared_neg = np.zeros(n_ev)
    c_shared_neg[idx_ap] = -1.0
    c_shared_neg[idx_qu] = -1.0
    contrasts.append(("shared_neg", c_shared_neg))

    return kept, design, ev_names, contrasts, diag


# ======================================================================
# Functional connectome analysis (XCP-D)
# ======================================================================
def _load_xcpd_relmat(path: Path) -> np.ndarray:
    """Robustly load an XCP-D LINC relmat TSV/CSV.

    Handles three formats in priority order:
      1. Pure numeric (legacy XCP-D, square matrix, no labels)
      2. LINC mode: header row of node names + first-column row labels
         (np.loadtxt cannot parse this — uses pandas)
      3. Header row only, no row labels
    """
    sep = "\t" if str(path).endswith(".tsv") else ","
    # Format 1: pure numeric — fastest, works for legacy outputs
    try:
        mat = np.loadtxt(path, delimiter=sep)
        if mat.ndim == 2 and mat.shape[0] == mat.shape[1] and mat.shape[0] > 0:
            return mat
    except Exception:
        pass
    # Format 2: LINC — header row + first-col row labels (pandas)
    try:
        import pandas as pd
        df = pd.read_csv(path, sep=sep, index_col=0)
        mat = df.to_numpy(dtype=np.float64)
        if mat.ndim == 2 and mat.shape[0] == mat.shape[1] and mat.shape[0] > 0:
            return mat
    except Exception:
        pass
    # Format 3: header row only, no row labels
    mat = np.loadtxt(path, delimiter=sep, skiprows=1)
    return mat


def load_xcpd_connectomes(subjects, atlas):
    """Load XCP-D connectomes for the EXACT atlas name, rejecting any subject
    whose matrix name doesn't match. Returns dict subj -> (group, matrix)."""
    matrices = {}
    skipped = 0
    load_errors = 0
    for subj, grp in subjects:
        xcp_dir = _xcpd_subject_dir(subj)
        # XCP-D naming: *task-rest*_atlas-<atlas>_*relmat.tsv (or conmat)
        pattern_candidates = [
            f"*task-rest*_atlas-{atlas}*_*relmat.tsv",
            f"*task-rest*_atlas-{atlas}*_*relmat.csv",
            f"*task-rest*_atlas-{atlas}*_*measure-pearsoncorrelation*.tsv",
        ]
        found = None
        for pat in pattern_candidates:
            matches = list(xcp_dir.rglob(pat))
            if matches:
                # Prefer the longest-run matrix when multiple exist (the
                # real ~240-vol rest scan, not a 10-vol reverse-PE run).
                # Sort by file size (larger = longer run = more edges).
                matches.sort(key=lambda p: p.stat().st_size, reverse=True)
                found = matches[0]
                break
        if not found:
            skipped += 1
            continue
        try:
            mat = _load_xcpd_relmat(found)
            if mat.ndim != 2 or mat.shape[0] != mat.shape[1]:
                raise ValueError(f"non-square {mat.shape}")
            matrices[subj] = (grp, mat)
        except Exception as e:
            load_errors += 1
            if load_errors <= 5:
                print(f"  WARNING: could not load {found.name}: {e}",
                      file=sys.stderr)
            skipped += 1
    if load_errors > 5:
        print(f"  ... and {load_errors - 5} more load errors", file=sys.stderr)
    if skipped:
        print(f"  {skipped} subjects had no loadable {atlas} connectome",
              file=sys.stderr)
    return matrices


# ----------------------------------------------------------------------
# Continuous connectome analysis (apathy × ICB)
# ----------------------------------------------------------------------
def _nbs_continuous(all_mats: np.ndarray, design: np.ndarray,
                    contrast: "int | np.ndarray",
                    cov_cols: "list[int] | None" = None,
                    t_thresh: float = 3.1, nperm: int = 5000,
                    seed: int = 42):
    """Permutation-based NBS for a CONTINUOUS regressor of interest.

    Uses residualised-regressor permutation (Freedman-Lane type; Winkler et al. 2014): for each edge,
    regress out ALL other regressors (the reduced model), then permute
    the residualised contrast projection and recompute the partial
    t-statistic.  This is equivalent to the partial t from the full GLM.

    Parameters
    ----------
    all_mats : (n_subj, n_rois, n_rois)   Fisher-z FC matrices
    design   : (n_subj, n_ev)             full design matrix
    contrast : int or np.ndarray          either a column index (legacy)
               or a contrast vector of length n_ev (e.g. [1,-1,0,...])
    cov_cols : list[int] or None          DEPRECATED — ignored.  The
               reduced model is now derived automatically from the
               contrast and the full design (= all columns not in the
               contrast direction).
    """
    rng = np.random.default_rng(seed)
    n_subj, n_rois, _ = all_mats.shape
    n_ev = design.shape[1]
    triu = np.triu_indices(n_rois, k=1)
    n_edges = len(triu[0])
    edges = all_mats[:, triu[0], triu[1]]   # (n_subj, n_edges)

    # Build contrast vector
    if isinstance(contrast, (int, np.integer)):
        cvec = np.zeros(n_ev)
        cvec[contrast] = 1.0
    else:
        cvec = np.asarray(contrast, dtype=np.float64)

    # Contrast projection (the regressor of interest)
    x = design @ cvec

    # Reduced model = full design projected out of the contrast direction.
    # This is ALL columns that are orthogonal to the contrast, which is
    # equivalent to including all other regressors (primary + nuisance)
    # as covariates — matching what FSL randomise does.
    cvec_norm = cvec / np.sqrt((cvec ** 2).sum())
    # Project each design column onto the contrast direction and subtract
    X_reduced = design - np.outer(design @ cvec_norm, cvec_norm)
    # Keep only columns with non-zero variance (some may collapse)
    col_var = X_reduced.var(axis=0)
    keep_cols = np.where(col_var > 1e-12)[0]
    X_nui = X_reduced[:, keep_cols]

    # Residualise each edge AND the regressor on the reduced model
    # (permute the residualised regressor, not the raw data)
    def _residualise(y, X_nui):
        if X_nui.shape[1] == 0:
            return y - y.mean()
        X_int = np.column_stack([np.ones(len(y)), X_nui])
        beta, *_ = np.linalg.lstsq(X_int, y, rcond=None)
        return y - X_int @ beta

    x_res = _residualise(x, X_nui)

    edges_res = np.zeros_like(edges)
    for e in range(n_edges):
        edges_res[:, e] = _residualise(edges[:, e], X_nui)

    def _partial_tstat(x_vec: np.ndarray) -> np.ndarray:
        """t statistic for each edge ~ x_vec (after nuisance residualisation)."""
        x_c = x_vec - x_vec.mean()
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

    def _max_component_size(tvals: np.ndarray, sign: str):
        mask = tvals > t_thresh if sign == "pos" else tvals < -t_thresh
        if not np.any(mask):
            return 0, None
        adj = np.zeros((n_rois, n_rois), dtype=bool)
        sel_edges = np.where(mask)[0]
        for e in sel_edges:
            r, c = triu[0][e], triu[1][e]
            adj[r, c] = True
            adj[c, r] = True
        visited = np.zeros(n_rois, dtype=bool)
        components = []
        for start in range(n_rois):
            if visited[start] or not adj[start].any():
                continue
            stack = [start]
            comp = []
            while stack:
                node = stack.pop()
                if visited[node]:
                    continue
                visited[node] = True
                comp.append(node)
                stack.extend(np.where(adj[node] & ~visited)[0].tolist())
            comp_set = set(comp)
            n_comp_edges = sum(
                1 for e in sel_edges
                if triu[0][e] in comp_set and triu[1][e] in comp_set)
            if n_comp_edges > 0:
                comp_edges = [[int(triu[0][e]), int(triu[1][e])]
                              for e in sel_edges
                              if triu[0][e] in comp_set and triu[1][e] in comp_set]
                components.append((n_comp_edges, sorted(comp), comp_edges))
        if not components:
            return 0, None
        components.sort(key=lambda x: -x[0])
        return components[0][0], components

    # Observed
    obs_t = _partial_tstat(x_res)
    obs_pos_size, obs_pos_comps = _max_component_size(obs_t, "pos")
    obs_neg_size, obs_neg_comps = _max_component_size(obs_t, "neg")

    # Permutation null (shuffle the residualised regressor)
    null_pos = np.zeros(nperm, dtype=np.int32)
    null_neg = np.zeros(nperm, dtype=np.int32)
    for p in range(nperm):
        x_perm = rng.permutation(x_res)
        t_perm = _partial_tstat(x_perm)
        null_pos[p], _ = _max_component_size(t_perm, "pos")
        null_neg[p], _ = _max_component_size(t_perm, "neg")

    def _pvalue(size, null):
        return float((null >= size).sum() + 1) / (nperm + 1)

    return {
        "nperm": nperm,
        "t_thresh": t_thresh,
        "obs_components_pos": [
            (size, nodes, _pvalue(size, null_pos), edges)
            for size, nodes, edges in (obs_pos_comps or [])
        ],
        "obs_components_neg": [
            (size, nodes, _pvalue(size, null_neg), edges)
            for size, nodes, edges in (obs_neg_comps or [])
        ],
        "null_max_pos_p95": float(np.percentile(null_pos, 95)),
        "null_max_neg_p95": float(np.percentile(null_neg, 95)),
    }


def run_connectome_continuous(subjects, cov_file, nperm, atlas,
                              extended_covs: list[str] | None = None,
                              quip_binary: bool = False,
                              t_thresh: float = 3.1):
    """Continuous (dimensional) FC connectome analysis.

    Tests main effects of apathy_composite, quip, and their interaction
    with nuisance covariates (depression_factor, age, sex, site, TIV, mean_FD).
    Outputs go to GROUP_DIR/fc_connectome_continuous/ so the categorical
    results are preserved.
    """
    print(f"\n=== Continuous Connectome Analysis "
          f"(apathy × ICB, atlas={atlas}) ===")
    out_dir = GROUP_DIR / "fc_connectome_continuous"
    out_dir.mkdir(parents=True, exist_ok=True)

    matrices = load_xcpd_connectomes(subjects, atlas)
    if len(matrices) < 10:
        print(f"  Only {len(matrices)} subjects with {atlas} connectome")
        return

    ref_shape = list(matrices.values())[0][1].shape
    aligned = [(s, g, m) for s, (g, m) in matrices.items()
               if m.shape == ref_shape]
    print(f"  {len(aligned)} subjects × {ref_shape[0]} ROIs")

    subjs_only = [(s, g) for s, g, _ in aligned]

    try:
        kept, design, ev_names, contrasts, diag = build_continuous_design(
            Path(cov_file), subjs_only, include_depression=True,
            extended_covs=extended_covs, quip_binary=quip_binary)
    except Exception as e:
        print(f"  ERROR: {e}")
        return

    # Print VIFs + score ranges before running expensive permutations
    print(f"\n  Subjects with complete continuous regressors: {diag['n_kept']}"
          f"  (dropped {diag['n_dropped']})")
    print(f"  apathy_composite: min={diag['score_ranges']['apathy_composite'][0]:.2f} "
          f"max={diag['score_ranges']['apathy_composite'][1]:.2f} "
          f"mean={diag['score_ranges']['apathy_composite'][2]:.2f} "
          f"sd={diag['score_ranges']['apathy_composite'][3]:.2f}")
    print(f"  quip:             min={diag['score_ranges']['quip'][0]:.2f} "
          f"max={diag['score_ranges']['quip'][1]:.2f} "
          f"mean={diag['score_ranges']['quip'][2]:.2f} "
          f"sd={diag['score_ranges']['quip'][3]:.2f}")
    print(f"  VIFs (primary regressors): {diag['vifs_primary']}")
    any_high = any(v > 5 for v in diag['vifs_primary'].values()
                   if np.isfinite(v))
    if any_high:
        print(f"  WARNING: at least one VIF > 5 — interpret cautiously "
              f"(high collinearity between main effects and interaction)",
              file=sys.stderr)

    kept_set = {s for s, _ in kept}
    aligned = [(s, g, m) for s, g, m in aligned if s in kept_set]
    order = {s: i for i, (s, _) in enumerate(kept)}
    aligned.sort(key=lambda x: order[x[0]])

    n_subj = len(aligned)
    n_rois = ref_shape[0]
    all_mats = np.stack([m for _, _, m in aligned])

    # NaN/inf hardening
    n_bad = int((~np.isfinite(all_mats)).sum())
    if n_bad:
        print(f"  NOTE: {n_bad} non-finite entries in raw matrices "
              f"(replaced with 0)", file=sys.stderr)
        all_mats = np.nan_to_num(all_mats, nan=0.0, posinf=0.0, neginf=0.0)

    all_mats_z = np.arctanh(np.clip(all_mats, -0.999, 0.999))
    for i in range(n_subj):
        np.fill_diagonal(all_mats_z[i], 0)
    all_mats_z = np.nan_to_num(all_mats_z, nan=0.0, posinf=0.0, neginf=0.0)

    summary = {}
    for cname, cvec in contrasts:
        print(f"\n  NBS: {cname} (continuous, {nperm} perms)")
        res = _nbs_continuous(all_mats_z, design, cvec,
                              t_thresh=t_thresh,
                              nperm=nperm, seed=42)
        for label, key in [("pos", "obs_components_pos"),
                           ("neg", "obs_components_neg")]:
            comps = res[key][:5]
            if comps:
                print(f"    Top {label} components:")
                for size, nodes, p, edges in comps:
                    star = " *" if p < 0.05 else ""
                    print(f"      {size} edges, {len(nodes)} nodes, "
                          f"p={p:.4f}{star}")
        # Save JSON
        with open(out_dir / f"nbs_{cname}.json", "w") as f:
            json.dump({
                "contrast": cname,
                "nperm": res["nperm"],
                "t_thresh": res["t_thresh"],
                "null_max_pos_p95": res["null_max_pos_p95"],
                "null_max_neg_p95": res["null_max_neg_p95"],
                "obs_components_pos": [
                    {"n_edges": s, "nodes": n, "p": p, "edges": e}
                    for s, n, p, e in res["obs_components_pos"]],
                "obs_components_neg": [
                    {"n_edges": s, "nodes": n, "p": p, "edges": e}
                    for s, n, p, e in res["obs_components_neg"]],
            }, f, indent=2)
        summary[cname] = {
            "n_sig_pos": sum(1 for _, _, p, _ in res["obs_components_pos"]
                             if p < 0.05),
            "n_sig_neg": sum(1 for _, _, p, _ in res["obs_components_neg"]
                             if p < 0.05),
        }

    # Save the full design + diagnostics
    np.savetxt(out_dir / "design_matrix.tsv", design,
               delimiter="\t", header="\t".join(ev_names), comments="")
    with open(out_dir / "design_subjects.tsv", "w") as f:
        f.write("subject_id\tgroup\n")
        for sid, grp in kept:
            f.write(f"sub-{sid}\t{grp}\n")
    with open(out_dir / "design_diagnostics.json", "w") as f:
        json.dump(diag, f, indent=2, default=str)

    primary_names = {"apathy_composite_c", "quip_c", "apathy_x_quip"}
    with open(out_dir / "fc_connectome_continuous_summary.json", "w") as f:
        json.dump({
            "atlas": atlas,
            "n_subjects": n_subj,
            "n_rois": int(n_rois),
            "ev_names": ev_names,
            "primary_regressors": list(primary_names),
            "nbs_summary": summary,
            "vifs": diag["vifs_primary"],
            "score_ranges": diag["score_ranges"],
        }, f, indent=2)
    print(f"\n  Results: {out_dir}/")


# ======================================================================
# Main
# ======================================================================
def main():
    # Late-rebind: allow overriding FMRIPREP_ROOT from the CLI.
    # Must be declared before any read of the global inside this function.
    global FMRIPREP_ROOT

    parser = argparse.ArgumentParser(description="PPMI Group FC & Connectomics")
    parser.add_argument("--analysis", default="connectome_continuous",
                        choices=["connectome_continuous"])
    parser.add_argument("--include-extended", action="store_true",
                        help="Shortcut: include ALL four extended nuisance "
                             "covariates (disease_duration, UPDRS-III, LEDD, MoCA). "
                             "Equivalent to --extended-covs "
                             "disease_duration_months,updrs_iii_off,ledd,moca.")
    parser.add_argument("--extended-covs", default="",
                        help="Comma-separated subset of extended nuisance "
                             "covariates to add to the lean model. "
                             "Valid names: "
                             f"{','.join(EXTENDED_COVS_VALID)}. "
                             "Empty (default) = lean model "
                             "(age/sex/site/TIV/FD/depression only). "
                             "Drops subjects missing any requested covariate.")
    parser.add_argument("--quip-binary", action="store_true",
                        help="Dichotomise QUIP to presence/absence (>0 vs =0) "
                             "before use as a regressor. Sensitivity for the "
                             "PPMI cohort where QUIP is severely zero-inflated "
                             "(77%% of subjects have QUIP=0).")
    parser.add_argument("--nperm", type=int, default=10000)
    parser.add_argument("--t-thresh", type=float, default=3.1,
                        help="NBS primary threshold (t-statistic, default 3.1)")
    parser.add_argument("--atlas", default=XCPD_ATLAS,
                        help="XCP-D atlas name (e.g. 4S256Parcels)")
    parser.add_argument("--covariates",
                        default=str(SCRIPTS_DIR / "covariates.csv"))
    parser.add_argument("--fmriprep-root", default=str(FMRIPREP_ROOT),
                        help="Root searched for XCP-D outputs after DERIV_DIR "
                             "(xcpd_<subj>/sub-<subj>/ or xcp_d/sub-<subj>/)")
    parser.add_argument("--subjects-file", default="",
                        help="Optional subject|group list (e.g. from "
                             "validate_functional_outputs.py). When given, "
                             "overrides automatic XCP-D subject discovery.")
    args = parser.parse_args()

    # Resolve extended covariates (CLI mutex: --include-extended wins)
    if args.include_extended:
        extended_covs = list(EXTENDED_COVS_VALID)
    else:
        extended_covs = [c.strip() for c in args.extended_covs.split(",")
                         if c.strip()]
        invalid = [c for c in extended_covs if c not in EXTENDED_COVS_VALID]
        if invalid:
            print(f"ERROR: --extended-covs contains unknown names: {invalid}\n"
                  f"  valid: {EXTENDED_COVS_VALID}", file=sys.stderr)
            sys.exit(1)

    FMRIPREP_ROOT = Path(args.fmriprep_root)

    print("============================================")
    print("PPMI Group FC & Connectomics Analysis")
    print("============================================")
    print(f"  DERIV_DIR:     {DERIV_DIR}")
    print(f"  FMRIPREP_ROOT: {FMRIPREP_ROOT}")
    print(f"  GROUP_DIR:     {GROUP_DIR}")
    print(f"  Atlas:         {args.atlas}")
    print(f"  Nperm:         {args.nperm}")
    print(f"  Covariates:    {args.covariates}")
    print(f"  Extended covs: "
          f"{','.join(extended_covs) if extended_covs else '(none — lean model)'}")
    print(f"  QUIP binary:   {args.quip_binary}")

    if args.subjects_file:
        path = Path(args.subjects_file)
        if not path.is_file():
            print(f"ERROR: --subjects-file not found: {path}", file=sys.stderr)
            sys.exit(1)
        xcp_subjects = []
        with open(path) as f:
            for line in f:
                line = line.strip()
                if not line or "|" not in line:
                    continue
                sid, grp = line.split("|", 1)
                sid = sid.replace("sub-", "").strip()
                xcp_subjects.append((sid, grp.strip()))
        print(f"  Using subjects from: {path}  ({len(xcp_subjects)} subjects)")
    else:
        xcp_subjects = find_xcp_subjects()
    print(f"\nSubjects with XCP-D rest:    {len(xcp_subjects)}")

    # Continuous (dimensional) analyses — apathy × ICB interaction
    if args.analysis == "connectome_continuous":
        run_connectome_continuous(xcp_subjects, args.covariates,
                                   args.nperm, args.atlas,
                                   extended_covs=extended_covs,
                                   quip_binary=args.quip_binary,
                                   t_thresh=args.t_thresh)

    print("\n============================================")
    print("Analysis complete.")
    print(f"Results in: {GROUP_DIR}/")


if __name__ == "__main__":
    main()

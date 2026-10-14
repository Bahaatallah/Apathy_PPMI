#!/usr/bin/env python3
# Role: NBS engine (connectome_continuous analysis): edge-wise GLM, residualised-regressor
#     permutation, component FWE, 10,000 permutations, |t| > 3.1 by default. Reads the largest
#     4S256Parcels relmat per participant.
# Manuscript item: Table 2, Figure 3 and the NBS sensitivity analyses of the Supplement.
"""
nbs_group_analysis.py — Group functional connectivity & connectomics analysis

Analyses XCP-D denoised outputs and structural-connectome outputs for group
comparisons between Apathy_BL, Both_BL, and ICB_BL.

Key design decisions:
  * Atlas selection is explicit (XCPD_ATLAS env var or --atlas). The script
    refuses to mix matrices from different atlases.
  * Seed-based FC uses XCP-D denoised residuals (*desc-denoised*bold.nii.gz),
    NOT raw fMRIPrep preproc_bold. Resting runs only (task-rest).
  * Seed ROI masks are built from MNI coordinates using the image affine
    (nibabel), not voxel-grid arithmetic.
  * Seed FC is a true Pearson correlation map (z-scored seed × z-scored voxel),
    then Fisher-z transformed for group statistics.
  * Connectome-level inference uses real permutation-based NBS
    (label shuffling -> max component size null).
  * TR is read from each subject's JSON sidecar, never hard-coded.
  * Nuisance covariates (age, sex, TIV, mean_FD, site) are applied as GLM
    regressors in both voxelwise and edge-wise analyses.
  * Subjects must be present in covariates.csv or they are dropped (no
    silent imputation).

Analyses:
  1. Functional connectome comparison (XCP-D connectivity matrices)
       - Edge-wise GLM with nuisance covariates
       - Permutation-based NBS
       - Graph theory summary (bctpy)
  2. Seed-based FC from XCP-D denoised BOLD
       - Pearson correlation -> Fisher z
       - Voxelwise randomise with TFCE
  3. Structural connectome group analysis (from 09_tractography_connectome.sh)

Usage:
  python3 nbs_group_analysis.py --analysis all|connectome|seed|structural
                         [--nperm 10000]
                         [--atlas 4S256Parcels]
                         [--covariates covariates.csv]

Requires: numpy, scipy, pandas, nibabel, bctpy (optional for graph metrics)
"""

from __future__ import annotations

import argparse
import csv
import json
import os
import subprocess
import sys
from glob import glob
from pathlib import Path
from typing import Iterable

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
FSLDIR = Path(os.environ.get("FSLDIR", "/opt/fsl"))

# fMRIPrep outputs may live on ephemeral (05b keeps per-subject dirs
# `fmriprep_<id>/sub-<id>/` on EPHEMERAL to save home quota). FMRIPREP_ROOT
# is where we look first; subject dirs are then searched both under
# FMRIPREP_ROOT/fmriprep/sub-<id> (legacy) and
# FMRIPREP_ROOT/fmriprep_<id>/sub-<id> (05b).
FMRIPREP_ROOT = Path(os.environ.get(
    "FMRIPREP_ROOT",
    os.environ.get("EPHEMERAL_DERIV", f"{RDS_EPHEMERAL}/PPMI/derivatives")))

GROUPS = ["Apathy_BL", "Both_BL", "ICB_BL", "Neither_BL"]
# Session timepoint. Default BL; for 24M analysis set SESSION=ses-24M in env.
SESSION = os.environ.get("SESSION", "ses-BL")
XCPD_ATLAS = os.environ.get("XCPD_ATLAS", "4S256Parcels")


def _fmriprep_subject_dir(subj: str) -> Path:
    """Return the first fMRIPrep directory that actually exists for `subj`.

    Search order:
      1. $FMRIPREP_ROOT/fmriprep_<subj>/sub-<subj>       (05b on ephemeral)
      2. $FMRIPREP_ROOT/fmriprep/sub-<subj>              (legacy / non-split)
      3. $DERIV_DIR/fmriprep/sub-<subj>                  (home fallback)
    Falls back to the legacy home path if nothing exists so rglob still works.
    """
    candidates = [
        FMRIPREP_ROOT / f"fmriprep_{subj}" / f"sub-{subj}",
        FMRIPREP_ROOT / "fmriprep" / f"sub-{subj}",
        DERIV_DIR / "fmriprep" / f"sub-{subj}",
    ]
    for c in candidates:
        if c.is_dir():
            return c
    return DERIV_DIR / "fmriprep" / f"sub-{subj}"


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
REQUIRED_COVS = ["age", "sex", "TIV", "mean_FD", "site"]
CATEGORICAL = {"sex", "site", "scanner_model", "handedness"}

# MNI-space seed coordinates (mm) for PD apathy / ICB circuits
SEED_ROIS = {
    "ventral_striatum_L": (-10, 12, -6),
    "ventral_striatum_R": (10, 12, -6),
    "dorsal_caudate_L":   (-13, 15, 9),
    "dorsal_caudate_R":   (13, 15, 9),
    "putamen_L":          (-24, 4, -2),
    "putamen_R":          (24, 4, -2),
    "NAcc_L":             (-9, 9, -8),
    "NAcc_R":             (9, 9, -8),
    "VTA":                (0, -16, -8),
    "SN":                 (0, -20, -10),
    "ACC":                (0, 35, 14),
    "OFC_L":              (-24, 36, -12),
    "OFC_R":              (24, 36, -12),
}
SEED_RADIUS_MM = 5.0


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


def find_tractography_subjects() -> list[tuple[str, str]]:
    subjects = []
    for grp in GROUPS:
        grp_dir = BIDS_ROOT / grp
        if not grp_dir.is_dir():
            continue
        for sd in sorted(grp_dir.glob("sub-*")):
            subj = sd.name.replace("sub-", "")
            if (DERIV_DIR / "tractography" / f"sub-{subj}" / SESSION
                / "connectome_sift2.csv").is_file():
                subjects.append((subj, grp))
    return subjects


# ======================================================================
# Covariates + design matrix
# ======================================================================
# Delegates to the shared helper so FC uses exactly the same required
# covariate list, provenance integration, categorical handling, and
# group-relabel-after-filter policy as TBSS/VBM/morphometry/surface.
sys.path.insert(0, str(SCRIPTS_DIR))
from group_covariates import load_covariates_with_provenance  # noqa: E402


def load_covariates_strict(cov_file: Path, subjects: list[tuple[str, str]],
                            cov_list: list[str] | None = None):
    """
    Thin wrapper over group_covariates.load_covariates_with_provenance
    that keeps the original 4-tuple return signature used by the rest
    of this script (kept, design, ev_names, group_labels).
    """
    kept, design, ev_names, group_labels, _nuisance_cols = \
        load_covariates_with_provenance(
            cov_file=cov_file,
            subjects=subjects,
            deriv_dir=DERIV_DIR,
            extra_cov_list=cov_list,
            verbose=True)
    return kept, design, ev_names, group_labels


def group_contrast_vectors(ev_names: list[str], group_labels: list[str]) -> list[tuple[str, np.ndarray]]:
    """Return list of (name, contrast_vector) for pairwise group comparisons,
    with zero weights on all covariate EVs. Contrast length = len(ev_names)."""
    n_ev = len(ev_names)
    idx = {name: i for i, name in enumerate(ev_names)}
    contrasts = []
    for a in group_labels:
        for b in group_labels:
            if a == b:
                continue
            c = np.zeros(n_ev)
            c[idx[f"group_{a}"]] = 1.0
            c[idx[f"group_{b}"]] = -1.0
            contrasts.append((f"{a}_gt_{b}", c))
    return contrasts


# ======================================================================
# Continuous (dimensional) design matrix
# ======================================================================
# Uses raw clinical scores as continuous regressors rather than dichotomised
# group labels. More power + tests apathy × ICB interaction directly.

CONTINUOUS_PRIMARY = ["apathy_composite", "quip"]          # main effects
CONTINUOUS_INTERACTION = ("apathy_composite", "quip")      # product term
CONTINUOUS_REQUIRED = CONTINUOUS_PRIMARY + ["depression_factor",
                                             "age", "sex", "site", "TIV", "mean_FD"]


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
# Utility: residualise data on covariates (for edge-wise + NBS)
# ======================================================================
def residualise(y: np.ndarray, X: np.ndarray, idx_keep: Iterable[int]) -> np.ndarray:
    """
    Residualise y (n,) on X[:, idx_keep] using least squares.
    Returns residuals of same shape as y.
    idx_keep are the covariate columns only (group columns should be excluded).
    """
    X_sub = X[:, list(idx_keep)]
    if X_sub.shape[1] == 0:
        return y - y.mean()
    # Append intercept
    X_int = np.column_stack([np.ones(len(y)), X_sub])
    beta, *_ = np.linalg.lstsq(X_int, y, rcond=None)
    return y - X_int @ beta


# ======================================================================
# 1. Functional connectome analysis (XCP-D)
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


def run_nbs(all_mats: np.ndarray, design: np.ndarray,
            contrast: np.ndarray, cov_cols: list[int],
            t_thresh: float = 3.1, nperm: int = 5000, seed: int = 0):
    """
    Permutation-based Network-Based Statistic (Zalesky 2010).

    Parameters
    ----------
    all_mats : (n_subj, n_rois, n_rois)   Fisher-z FC matrices
    design   : (n_subj, n_ev)             design matrix (group + covariates)
    contrast : (n_ev,)                    pairwise group contrast vector
    cov_cols : list[int]                  indices of covariate columns in design
    t_thresh : float                      primary edge-level t-statistic threshold
    nperm    : int                        number of permutations
    seed     : int                        RNG seed

    Returns
    -------
    dict with keys: 'obs_components_pos', 'obs_components_neg'
                    (lists of (n_edges, nodes, p, edges)),
                    'null_max_pos', 'null_max_neg', 'nperm', 't_thresh'
    """
    rng = np.random.default_rng(seed)
    n_subj, n_rois, _ = all_mats.shape
    triu = np.triu_indices(n_rois, k=1)
    n_edges = len(triu[0])

    # Residualise each edge on covariates
    edges = all_mats[:, triu[0], triu[1]]  # (n_subj, n_edges)
    if cov_cols:
        edges_res = np.zeros_like(edges)
        for e in range(n_edges):
            edges_res[:, e] = residualise(edges[:, e], design, cov_cols)
        edges = edges_res

    # Group columns: infer from the two non-zero contrast positions
    pos_col = int(np.argmax(contrast))
    neg_col = int(np.argmin(contrast))
    idx_a = np.where(design[:, pos_col] > 0.5)[0]
    idx_b = np.where(design[:, neg_col] > 0.5)[0]
    if len(idx_a) < 3 or len(idx_b) < 3:
        return None

    def edge_tstats(sel_a, sel_b):
        mean_a = edges[sel_a].mean(axis=0)
        mean_b = edges[sel_b].mean(axis=0)
        var_a = edges[sel_a].var(axis=0, ddof=1)
        var_b = edges[sel_b].var(axis=0, ddof=1)
        n_a = len(sel_a)
        n_b = len(sel_b)
        se = np.sqrt(var_a / n_a + var_b / n_b + 1e-12)
        return (mean_a - mean_b) / se

    def max_component_size(tvals, sign):
        """Find largest connected component among supra-threshold edges."""
        if sign == "pos":
            mask = tvals > t_thresh
        else:
            mask = tvals < -t_thresh
        if not np.any(mask):
            return 0, None
        # Build adjacency (n_rois x n_rois) from suprathreshold edges
        adj = np.zeros((n_rois, n_rois), dtype=bool)
        sel_edges = np.where(mask)[0]
        for e in sel_edges:
            r, c = triu[0][e], triu[1][e]
            adj[r, c] = True
            adj[c, r] = True
        # Connected components via BFS
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
            # Count edges within the component
            comp_set = set(comp)
            n_comp_edges = sum(
                1 for e in sel_edges
                if triu[0][e] in comp_set and triu[1][e] in comp_set
            )
            if n_comp_edges > 0:
                comp_edges = [[int(triu[0][e]), int(triu[1][e])]
                              for e in sel_edges
                              if triu[0][e] in comp_set and triu[1][e] in comp_set]
                components.append((n_comp_edges, sorted(comp), comp_edges))
        if not components:
            return 0, None
        components.sort(key=lambda x: -x[0])
        return components[0][0], components

    # Observed statistic
    obs_t = edge_tstats(idx_a, idx_b)
    obs_max_size_pos, obs_comps_pos = max_component_size(obs_t, "pos")
    obs_max_size_neg, obs_comps_neg = max_component_size(obs_t, "neg")

    # Permutation null (shuffle group assignment between the two groups only)
    combined = np.concatenate([idx_a, idx_b])
    n_a = len(idx_a)
    null_max_pos = np.zeros(nperm, dtype=np.int32)
    null_max_neg = np.zeros(nperm, dtype=np.int32)
    for p in range(nperm):
        perm = rng.permutation(combined)
        sel_a = perm[:n_a]
        sel_b = perm[n_a:]
        t_perm = edge_tstats(sel_a, sel_b)
        null_max_pos[p], _ = max_component_size(t_perm, "pos")
        null_max_neg[p], _ = max_component_size(t_perm, "neg")

    # P-values for observed components
    def p_for(size, null):
        return float((null >= size).sum() + 1) / (nperm + 1)

    out = {
        "nperm": nperm,
        "t_thresh": t_thresh,
        "obs_components_pos": [
            (size, nodes, p_for(size, null_max_pos), edges)
            for size, nodes, edges in (obs_comps_pos or [])
        ],
        "obs_components_neg": [
            (size, nodes, p_for(size, null_max_neg), edges)
            for size, nodes, edges in (obs_comps_neg or [])
        ],
        "null_max_pos": null_max_pos.tolist(),
        "null_max_neg": null_max_neg.tolist(),
    }
    return out


def run_connectome_analysis(subjects, cov_file, nperm, atlas, t_thresh: float = 3.1):
    print(f"\n=== Functional Connectome Analysis (atlas={atlas}) ===")
    out_dir = GROUP_DIR / "fc_connectome"
    out_dir.mkdir(parents=True, exist_ok=True)

    matrices = load_xcpd_connectomes(subjects, atlas)
    if len(matrices) < 10:
        print(f"  Only {len(matrices)} subjects with {atlas} connectome — need >= 10.")
        return

    # Align by shape
    ref_shape = list(matrices.values())[0][1].shape
    aligned = [(s, g, m) for s, (g, m) in matrices.items() if m.shape == ref_shape]
    print(f"  {len(aligned)} subjects × {ref_shape[0]} ROIs")

    subjs_only = [(s, g) for s, g, _ in aligned]
    try:
        kept, design, ev_names, group_labels = load_covariates_strict(Path(cov_file), subjs_only)
    except (FileNotFoundError, ValueError, RuntimeError) as e:
        print(f"  ERROR: {e}")
        return

    kept_set = {s for s, _ in kept}
    aligned = [(s, g, m) for s, g, m in aligned if s in kept_set]
    # Reorder to match `kept`
    order = {s: i for i, (s, _) in enumerate(kept)}
    aligned.sort(key=lambda x: order[x[0]])

    n_subj = len(aligned)
    n_rois = ref_shape[0]
    all_mats = np.stack([m for _, _, m in aligned])

    # Sanitise raw correlation matrices BEFORE Fisher-z. Partial-coverage
    # subjects can have NaN rows/columns where an atlas ROI fell entirely
    # outside the BOLD FOV, and arctanh(1.0) from the diagonal is inf.
    # Replace non-finite with 0 — safe because the diagonal is zeroed
    # cosmetically below and NaN-off-diagonal edges should not contribute.
    n_bad = int((~np.isfinite(all_mats)).sum())
    if n_bad:
        print(f"  NOTE: {n_bad} non-finite entries across raw matrices "
              f"(replaced with 0)", file=sys.stderr)
        all_mats = np.nan_to_num(all_mats, nan=0.0, posinf=0.0, neginf=0.0)

    # Fisher z with clipping to [-0.999, 0.999] so arctanh is finite
    all_mats_z = np.arctanh(np.clip(all_mats, -0.999, 0.999))
    for i in range(n_subj):
        np.fill_diagonal(all_mats_z[i], 0)
    # Final guarantee
    all_mats_z = np.nan_to_num(all_mats_z, nan=0.0, posinf=0.0, neginf=0.0)

    # Save group means
    for g in group_labels:
        idx = [i for i, (_, gg, _) in enumerate(aligned) if gg == g]
        if idx:
            np.savetxt(out_dir / f"mean_fc_{g}.csv",
                       all_mats_z[idx].mean(axis=0), delimiter=",", fmt="%.6f")

    # Covariate columns = everything except group_* EVs
    cov_cols = [i for i, n in enumerate(ev_names) if not n.startswith("group_")]

    # Run NBS for each pairwise contrast
    contrasts = group_contrast_vectors(ev_names, group_labels)
    summary = {}
    for cname, cvec in contrasts:
        # Only run unique pairs (A>B; skip B>A which is symmetric)
        a, b = cname.split("_gt_")
        if a > b:
            continue
        print(f"\n  NBS: {cname} ({nperm} perms)")
        res = run_nbs(all_mats_z, design, cvec, cov_cols,
                      t_thresh=t_thresh, nperm=nperm, seed=42)
        if res is None:
            continue
        # Log top components
        for label, key in [("pos", "obs_components_pos"), ("neg", "obs_components_neg")]:
            comps = res[key][:5]
            if comps:
                print(f"    Top {label} components:")
                for size, nodes, p, edges in comps:
                    star = " *" if p < 0.05 else ""
                    print(f"      {size} edges, {len(nodes)} nodes, p={p:.4f}{star}")
        # Save
        with open(out_dir / f"nbs_{cname}.json", "w") as f:
            # Don't dump the full null arrays; keep summary
            out = {k: v for k, v in res.items() if "null_" not in k}
            out["null_max_pos_mean"] = float(np.mean(res["null_max_pos"]))
            out["null_max_pos_p95"] = float(np.percentile(res["null_max_pos"], 95))
            out["null_max_neg_mean"] = float(np.mean(res["null_max_neg"]))
            out["null_max_neg_p95"] = float(np.percentile(res["null_max_neg"], 95))
            # Convert tuple components to lists for JSON
            out["obs_components_pos"] = [
                {"n_edges": s, "nodes": n, "p": p, "edges": e}
                for s, n, p, e in out["obs_components_pos"]
            ]
            out["obs_components_neg"] = [
                {"n_edges": s, "nodes": n, "p": p, "edges": e}
                for s, n, p, e in out["obs_components_neg"]
            ]
            json.dump(out, f, indent=2)
        summary[cname] = {
            "n_sig_pos": sum(1 for _, _, p, _ in res["obs_components_pos"] if p < 0.05),
            "n_sig_neg": sum(1 for _, _, p, _ in res["obs_components_neg"] if p < 0.05),
        }

    # Graph theory
    try:
        import bct
        print("\n  --- Graph theory metrics ---")
        import pandas as pd
        records = []
        skipped_graph: list[tuple[str, str, float, str]] = []
        # Sensitivity analysis across densities 10, 15, 20, 25%
        for density in [0.10, 0.15, 0.20, 0.25]:
            for i, (subj, grp, _) in enumerate(aligned):
                mat = all_mats_z[i].copy()
                # Sanitise: Fisher-z can leave +/- inf on the diagonal and
                # any numerical NaNs in the original correlation matrix
                # propagate to the graph metrics. zero them BEFORE
                # thresholding so bctpy never sees a non-finite value.
                if not np.all(np.isfinite(mat)):
                    n_bad = int((~np.isfinite(mat)).sum())
                    mat = np.nan_to_num(mat, nan=0.0, posinf=0.0, neginf=0.0)
                    # Only log the first occurrence per subject per density
                    if density == 0.10:
                        print(f"    sub-{subj}: replaced {n_bad} non-finite "
                              f"values before graph metrics", file=sys.stderr)
                mat[mat < 0] = 0
                np.fill_diagonal(mat, 0)
                try:
                    mat_thr = bct.threshold_proportional(mat, density)
                    mat_bin = (mat_thr > 0).astype(float)
                    if not (np.all(np.isfinite(mat_thr))
                            and np.all(np.isfinite(mat_bin))):
                        raise ValueError("non-finite after thresholding")
                    rec = {
                        "subject_id": f"sub-{subj}",
                        "group": grp,
                        "density": density,
                        "mean_degree": float(np.mean(bct.degrees_und(mat_bin))),
                        "mean_strength": float(np.mean(bct.strengths_und(mat_thr))),
                        "mean_clustering_bin": float(np.mean(bct.clustering_coef_bu(mat_bin))),
                        "mean_clustering_w": float(np.mean(bct.clustering_coef_wu(mat_thr))),
                        "char_path_len": float(bct.charpath(bct.distance_bin(mat_bin))[0]),
                        "modularity": float(bct.modularity_und(mat_thr)[1]),
                        "global_efficiency": float(bct.efficiency_bin(mat_bin)),
                    }
                except Exception as e:
                    skipped_graph.append((subj, grp, density, str(e)[:60]))
                    continue
                records.append(rec)
        df = pd.DataFrame(records)
        df.to_csv(out_dir / "graph_metrics_all_densities.csv", index=False)
        print(f"  Saved graph metrics (sensitivity sweep): {len(df)} rows "
              f"(skipped {len(skipped_graph)})")
        if skipped_graph:
            with open(out_dir / "graph_metrics_skipped.txt", "w") as f:
                for subj, grp, d, err in skipped_graph:
                    f.write(f"sub-{subj}\t{grp}\tdensity={d}\t{err}\n")
    except ImportError:
        print("  bctpy not installed — skipping graph metrics (pip install bctpy)")

    with open(out_dir / "fc_connectome_summary.json", "w") as f:
        json.dump({
            "atlas": atlas,
            "n_subjects": n_subj,
            "n_rois": int(n_rois),
            "ev_names": ev_names,
            "nbs_summary": summary,
        }, f, indent=2)
    print(f"\n  Results: {out_dir}/")


# ----------------------------------------------------------------------
# Continuous connectome analysis (apathy × ICB)
# ----------------------------------------------------------------------
def _nbs_continuous(all_mats: np.ndarray, design: np.ndarray,
                    contrast: "int | np.ndarray",
                    cov_cols: "list[int] | None" = None,
                    t_thresh: float = 3.1, nperm: int = 5000,
                    seed: int = 42):
    """Permutation-based NBS for a CONTINUOUS regressor of interest.

    Uses Freedman-Lane permutation on residualised data: for each edge,
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
    # (Freedman-Lane: permute the residualised regressor, not the raw data)
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

    # NaN/inf hardening (same as categorical path)
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
# 2. Seed-based FC (XCP-D denoised residuals)
# ======================================================================
def _build_sphere_mask_mni(ref_img, mni_xyz_mm, radius_mm):
    """Build a spherical ROI mask at MNI coordinate `mni_xyz_mm` within ref_img space."""
    import nibabel as nib
    from nibabel.affines import apply_affine

    affine = ref_img.affine
    inv = np.linalg.inv(affine)
    # Voxel coordinate of seed center
    center_vox = apply_affine(inv, np.array(mni_xyz_mm))
    shape = ref_img.shape[:3]

    # Build voxel grid and compute Euclidean distance (mm) to center
    ii, jj, kk = np.indices(shape)
    # Convert every voxel to mm via affine
    vox_coords = np.stack([ii.ravel(), jj.ravel(), kk.ravel(), np.ones(ii.size)], axis=0)
    mm_coords = affine @ vox_coords
    mm = mm_coords[:3].T  # (N, 3)
    dist = np.linalg.norm(mm - np.array(mni_xyz_mm), axis=1).reshape(shape)
    mask = (dist <= radius_mm).astype(np.int16)
    if mask.sum() == 0:
        return None
    return nib.Nifti1Image(mask, affine, ref_img.header)


def _read_tr_from_json(bold_nii: Path) -> float | None:
    json_file = bold_nii.with_name(bold_nii.name.replace(".nii.gz", ".json"))
    if not json_file.is_file():
        # Look for a sidecar in parent dirs (BIDS inheritance)
        for parent in bold_nii.parents:
            for candidate in parent.glob("*task-rest*_bold.json"):
                json_file = candidate
                break
            if json_file.is_file():
                break
    if json_file.is_file():
        try:
            with open(json_file) as f:
                meta = json.load(f)
            return float(meta.get("RepetitionTime", 0)) or None
        except Exception:
            return None
    return None


def _pearson_seed_map(bold_4d: np.ndarray, seed_ts: np.ndarray,
                       brain_mask: np.ndarray) -> np.ndarray:
    """Voxelwise Pearson correlation between a 1D seed timeseries and a 4D BOLD.
    Returns a 3D map of correlations (NaN outside brain_mask)."""
    if bold_4d.ndim != 4:
        raise ValueError(
            f"_pearson_seed_map expected 4D BOLD, got {bold_4d.ndim}D "
            f"shape {bold_4d.shape}")
    if brain_mask.ndim != 3:
        raise ValueError(
            f"_pearson_seed_map expected 3D mask, got {brain_mask.ndim}D "
            f"shape {brain_mask.shape}")
    if brain_mask.shape != bold_4d.shape[:3]:
        raise ValueError(
            f"_pearson_seed_map mask shape {brain_mask.shape} != "
            f"BOLD {bold_4d.shape[:3]}")
    n_tp = bold_4d.shape[3]
    # Z-score seed
    seed_z = (seed_ts - seed_ts.mean()) / (seed_ts.std(ddof=1) + 1e-12)
    # Z-score each voxel timecourse inside the mask
    data = bold_4d[brain_mask > 0].astype(np.float32)  # (n_vox, n_tp)
    data_z = (data - data.mean(axis=1, keepdims=True)) / \
             (data.std(axis=1, ddof=1, keepdims=True) + 1e-12)
    # Pearson = mean(seed_z * data_z)
    r = (data_z @ seed_z) / (n_tp - 1)
    out = np.full(bold_4d.shape[:3], np.nan, dtype=np.float32)
    out[brain_mask > 0] = r
    return out


def run_seed_analysis(subjects, cov_file, nperm):
    print("\n=== Seed-Based FC (XCP-D denoised input) ===")
    try:
        import nibabel as nib
    except ImportError:
        print("  ERROR: nibabel required (pip install nibabel)")
        return

    out_dir = GROUP_DIR / "fc_seed"
    out_dir.mkdir(parents=True, exist_ok=True)

    # Reference MNI template (for seed mask construction)
    mni_ref = FSLDIR / "data" / "standard" / "MNI152_T1_2mm_brain.nii.gz"
    if not mni_ref.is_file():
        print(f"  ERROR: MNI template not found: {mni_ref}")
        return
    ref_img = nib.load(str(mni_ref))

    # Build seed masks once in MNI space (they'll be resampled per subject)
    seed_dir = out_dir / "seeds_mni2mm"
    seed_dir.mkdir(exist_ok=True)
    seed_mask_imgs = {}
    for name, xyz in SEED_ROIS.items():
        mask_file = seed_dir / f"{name}.nii.gz"
        if not mask_file.is_file():
            img = _build_sphere_mask_mni(ref_img, xyz, SEED_RADIUS_MM)
            if img is None:
                print(f"  WARNING: empty mask for {name}")
                continue
            nib.save(img, str(mask_file))
        seed_mask_imgs[name] = mask_file
    print(f"  Built {len(seed_mask_imgs)} seed masks in MNI 2mm space")

    # Per subject: compute seed FC maps for each ROI
    subject_fc_maps = {name: [] for name in seed_mask_imgs}
    subject_by_seed = {name: [] for name in seed_mask_imgs}

    for subj, grp in subjects:
        xcp_dir = _xcpd_subject_dir(subj)
        bolds = list(xcp_dir.rglob("*task-rest*desc-denoised*bold.nii.gz"))
        # Exclude the smoothed variant (*_desc-denoisedSmoothed_bold.nii.gz)
        # when the plain denoised version is present — XCP-D seed FC should
        # run on unsmoothed residuals.
        plain = [b for b in bolds if "desc-denoisedSmoothed" not in b.name]
        if plain:
            bolds = plain
        if not bolds:
            continue
        # Pick the longest run (largest file = most volumes). PPMI has both
        # a ~240-vol real rest scan (dir-RL) and a ~10-vol reverse-PE run
        # (dir-LR) used for TopUp; we want the former.
        bolds.sort(key=lambda p: p.stat().st_size, reverse=True)
        bold_file = bolds[0]
        # Confirm space
        if "MNI152NLin2009cAsym" not in bold_file.name and \
           "MNI152NLin6Asym" not in bold_file.name:
            print(f"  WARNING: sub-{subj} BOLD not in MNI space — skipped")
            continue

        try:
            bold_img = nib.load(str(bold_file))
            bold_data = bold_img.get_fdata(dtype=np.float32)
        except Exception as e:
            print(f"  WARNING: sub-{subj} load failed: {e}")
            continue

        # Drop singleton dims and require 4D (X,Y,Z,T). Some XCP-D builds
        # emit a single-frame 3D NIfTI for empty runs — skip those.
        bold_data = np.squeeze(bold_data)
        if bold_data.ndim == 3:
            print(f"  WARNING: sub-{subj} BOLD is 3D (only 1 volume?) — skipped")
            continue
        if bold_data.ndim != 4:
            print(f"  WARNING: sub-{subj} BOLD has shape {bold_data.shape} — skipped")
            continue

        # Brain mask: load the authoritative brain mask that matches this BOLD.
        # XCP-D denoised residuals are mean-centred and can legitimately have
        # values of zero and negative values, so thresholding a single frame
        # would carve arbitrary holes out of the brain. Look for the BIDS
        # sibling mask first (XCP-D/fMRIPrep write *desc-brain_mask.nii.gz),
        # then fall back to fMRIPrep's output. If neither exists, skip the
        # subject — don't invent a mask.
        mask_candidates = [
            bold_file.with_name(
                bold_file.name.replace("_desc-denoised_bold.nii.gz",
                                        "_desc-brain_mask.nii.gz")),
            bold_file.with_name(
                bold_file.name.replace("_desc-denoised_bold.nii.gz",
                                        "_mask.nii.gz")),
        ]
        # Also look in the matching fMRIPrep directory (may be on ephemeral
        # in the 05b layout — see _fmriprep_subject_dir)
        fprep_dir = _fmriprep_subject_dir(subj)
        mask_candidates.extend(sorted(fprep_dir.rglob(
            "*task-rest*space-MNI152NLin2009cAsym*desc-brain_mask.nii.gz")))
        mask_candidates.extend(sorted(fprep_dir.rglob(
            "*task-rest*space-MNI152NLin6Asym*desc-brain_mask.nii.gz")))

        brain_mask = None
        for mcand in mask_candidates:
            if mcand and Path(mcand).is_file():
                try:
                    mimg = nib.load(str(mcand))
                    # Must match BOLD grid
                    if mimg.shape[:3] == bold_data.shape[:3] and \
                       np.allclose(mimg.affine, bold_img.affine, atol=1e-3):
                        brain_mask = (mimg.get_fdata() > 0.5).astype(np.int16)
                        break
                    # Otherwise try to resample
                    from nibabel.processing import resample_from_to
                    mresamp = resample_from_to(
                        mimg, (bold_data.shape[:3], bold_img.affine), order=0)
                    brain_mask = (mresamp.get_fdata() > 0.5).astype(np.int16)
                    break
                except Exception:
                    continue

        if brain_mask is None or brain_mask.sum() == 0:
            print(f"  WARNING: sub-{subj} no usable brain mask — skipped")
            continue

        # Force 3D. XCP-D/fMRIPrep sometimes emit masks with a singleton
        # time dim (X,Y,Z,1) — without squeeze, boolean-indexing a 4D BOLD
        # with a 4D mask collapses to 1D and crashes _pearson_seed_map.
        brain_mask = np.squeeze(brain_mask)
        if brain_mask.ndim != 3 or brain_mask.shape != bold_data.shape[:3]:
            print(f"  WARNING: sub-{subj} mask shape {brain_mask.shape} "
                  f"incompatible with BOLD {bold_data.shape[:3]} — skipped")
            continue

        for seed_name, mask_file in seed_mask_imgs.items():
            # Resample seed mask to BOLD grid (nearest)
            seed_img = nib.load(str(mask_file))
            if seed_img.shape != bold_data.shape[:3] or \
               not np.allclose(seed_img.affine, bold_img.affine):
                # Resample via nibabel.processing (avoid subprocess to fsl)
                try:
                    from nibabel.processing import resample_from_to
                    seed_resampled = resample_from_to(
                        seed_img, (bold_data.shape[:3], bold_img.affine), order=0)
                    seed_data = seed_resampled.get_fdata()
                except Exception:
                    continue
            else:
                seed_data = seed_img.get_fdata()

            seed_mask = seed_data > 0.5
            if seed_mask.sum() < 3:
                continue

            # Seed timeseries = mean over seed voxels
            seed_ts = bold_data[seed_mask].mean(axis=0)
            if np.std(seed_ts) < 1e-6:
                continue

            # True Pearson correlation map
            r_map = _pearson_seed_map(bold_data, seed_ts, brain_mask)
            # Fisher z
            z_map = np.arctanh(np.clip(r_map, -0.999, 0.999))
            z_map[np.isnan(z_map)] = 0

            out_file = out_dir / seed_name / f"sub-{subj}_fcz.nii.gz"
            out_file.parent.mkdir(parents=True, exist_ok=True)
            nib.save(nib.Nifti1Image(z_map, bold_img.affine, bold_img.header),
                     str(out_file))
            subject_fc_maps[seed_name].append(str(out_file))
            subject_by_seed[seed_name].append((subj, grp))

    # Group-level randomise for each seed
    for seed_name, file_list in subject_fc_maps.items():
        if len(file_list) < 10:
            print(f"  {seed_name}: {len(file_list)} subjects — skipping")
            continue
        print(f"\n  {seed_name}: {len(file_list)} subjects → randomise")

        # Build design for this seed's kept subjects
        try:
            kept, design, ev_names, group_labels = load_covariates_strict(
                Path(cov_file), subject_by_seed[seed_name])
        except Exception as e:
            print(f"    ERROR: {e}")
            continue

        # Reorder file_list to match kept
        subj_to_file = dict(zip([s for s, _ in subject_by_seed[seed_name]], file_list))
        ordered_files = [subj_to_file[s] for s, _ in kept if s in subj_to_file]

        merged = out_dir / f"all_{seed_name}_fcz.nii.gz"
        subprocess.run(
            ["fslmerge", "-t", str(merged)] + ordered_files, check=True)

        # Write FSL design files
        design_prefix = out_dir / f"design_{seed_name}"
        _write_fsl_design(design, ev_names, group_labels,
                          str(design_prefix), kept)

        # Brain mask from the 4D
        mask_file = out_dir / f"mask_{seed_name}.nii.gz"
        subprocess.run(["fslmaths", str(merged), "-abs", "-Tmean",
                        "-thr", "0.001", "-bin", str(mask_file)], check=False)

        # Randomise with TFCE (3D, not TBSS 2D)
        print(f"    Running randomise ({nperm} perms, -T)...")
        subprocess.run([
            "randomise", "-i", str(merged),
            "-o", str(out_dir / f"seed_{seed_name}"),
            "-d", f"{design_prefix}_design.mat",
            "-t", f"{design_prefix}_design.con",
            "-f", f"{design_prefix}_design.fts",
            "-m", str(mask_file),
            "-n", str(nperm),
            "-T", "-V",
        ], check=False)


def _write_fsl_design_continuous(design, ev_names, contrasts, prefix, kept):
    """Write FSL VEST .mat/.con/.fts files for a CONTINUOUS design.

    contrasts is a list of (name, np.ndarray) coming from
    build_continuous_design. An F-test spanning all contrasts is written
    as the omnibus test.
    """
    n, n_ev = design.shape
    with open(f"{prefix}_design.mat", "w") as f:
        f.write(f"/NumWaves\t{n_ev}\n/NumPoints\t{n}\n")
        f.write("/PPheights\t" + "\t".join(
            f"{max(np.max(np.abs(design[:,i])), 1e-6):.6f}" for i in range(n_ev)
        ) + "\n")
        f.write("\n/Matrix\n")
        for row in design:
            f.write("\t".join(f"{v:.6f}" for v in row) + "\n")

    with open(f"{prefix}_design.con", "w") as f:
        f.write(f"/NumWaves\t{n_ev}\n/NumContrasts\t{len(contrasts)}\n")
        for i, (name, _) in enumerate(contrasts, 1):
            f.write(f"/ContrastName{i}\t{name}\n")
        f.write("\n/Matrix\n")
        for _, c in contrasts:
            f.write("\t".join(f"{v:.6f}" for v in c) + "\n")

    # F-test spans all 3 primary contrasts
    with open(f"{prefix}_design.fts", "w") as f:
        f.write(f"/NumWaves\t{len(contrasts)}\n/NumContrasts\t1\n")
        f.write("\n/Matrix\n")
        f.write("\t".join("1" for _ in contrasts) + "\n")

    with open(f"{prefix}_subjects.txt", "w") as f:
        for sid, grp in kept:
            f.write(f"sub-{sid}\t{grp}\n")


def run_seed_continuous(subjects, cov_file, nperm,
                        extended_covs: list[str] | None = None,
                        quip_binary: bool = False):
    """Continuous seed-FC analysis using apathy/quip/interaction regressors.

    REUSES the per-subject FC maps produced by run_seed_analysis if they
    exist at GROUP_DIR/fc_seed/<seed>/sub-<id>_fcz.nii.gz. Otherwise
    computes them fresh (calling the same machinery as run_seed_analysis).
    Outputs go to GROUP_DIR/fc_seed_continuous/ so the categorical results
    are preserved.
    """
    print("\n=== Continuous Seed-Based FC (apathy × ICB) ===")
    try:
        import nibabel as nib  # noqa: F401
    except ImportError:
        print("  ERROR: nibabel required")
        return

    out_dir = GROUP_DIR / "fc_seed_continuous"
    out_dir.mkdir(parents=True, exist_ok=True)
    categorical_seed_dir = GROUP_DIR / "fc_seed"

    # Check whether per-subject FC maps already exist from the categorical
    # run. If so, we reuse them. If not, fall back to re-running the heavy
    # per-subject computation via run_seed_analysis machinery (expensive).
    reuse = categorical_seed_dir.is_dir() and any(
        (categorical_seed_dir / name).is_dir()
        for name in SEED_ROIS)
    if reuse:
        print(f"  Reusing FC maps from {categorical_seed_dir}/")
    else:
        print(f"  No pre-existing FC maps — running seed computation first")
        # Call the heavy function; it writes to fc_seed/
        run_seed_analysis(subjects, cov_file, nperm=1000)
        # After that, maps exist
        reuse = True

    # Build the continuous design against subjects with FC maps
    kept_candidates: list[tuple[str, str]] = []
    for subj, grp in subjects:
        # Pick one seed just to test for the subject's map
        any_seed = next(iter(SEED_ROIS))
        f = categorical_seed_dir / any_seed / f"sub-{subj}_fcz.nii.gz"
        if f.is_file():
            kept_candidates.append((subj, grp))

    if len(kept_candidates) < 10:
        print(f"  Only {len(kept_candidates)} subjects with FC maps — aborting")
        return

    try:
        kept, design, ev_names, contrasts, diag = build_continuous_design(
            Path(cov_file), kept_candidates, include_depression=True,
            extended_covs=extended_covs, quip_binary=quip_binary)
    except Exception as e:
        print(f"  ERROR: {e}")
        return

    print(f"  {diag['n_kept']} subjects in continuous design "
          f"(dropped {diag['n_dropped']})")
    print(f"  VIFs: {diag['vifs_primary']}")

    # Save the design once (same for all seeds)
    design_prefix = out_dir / "design"
    _write_fsl_design_continuous(design, ev_names, contrasts,
                                  str(design_prefix), kept)
    np.savetxt(out_dir / "design_matrix.tsv", design,
               delimiter="\t", header="\t".join(ev_names), comments="")
    with open(out_dir / "design_diagnostics.json", "w") as f:
        json.dump(diag, f, indent=2, default=str)

    kept_set = {s for s, _ in kept}

    for seed_name in SEED_ROIS:
        seed_src = categorical_seed_dir / seed_name
        if not seed_src.is_dir():
            print(f"  {seed_name}: no FC map directory — skipping")
            continue

        # Build ordered file list matching `kept`
        ordered_files = []
        for sid, _ in kept:
            f = seed_src / f"sub-{sid}_fcz.nii.gz"
            if f.is_file():
                ordered_files.append(str(f))
        if len(ordered_files) < 10:
            print(f"  {seed_name}: only {len(ordered_files)} FC maps in kept "
                  f"subject list — skipping")
            continue
        if len(ordered_files) != len(kept):
            print(f"  {seed_name}: {len(ordered_files)}/{len(kept)} FC maps "
                  f"available — skipping (design has wrong n)")
            continue

        print(f"\n  {seed_name}: {len(ordered_files)} subjects → randomise "
              f"(continuous, {nperm} perms)")

        merged = out_dir / f"all_{seed_name}_fcz.nii.gz"
        subprocess.run(
            ["fslmerge", "-t", str(merged)] + ordered_files, check=True)

        mask_file = out_dir / f"mask_{seed_name}.nii.gz"
        subprocess.run(["fslmaths", str(merged), "-abs", "-Tmean",
                        "-thr", "0.001", "-bin", str(mask_file)], check=False)

        subprocess.run([
            "randomise", "-i", str(merged),
            "-o", str(out_dir / f"seed_{seed_name}_cont"),
            "-d", f"{design_prefix}_design.mat",
            "-t", f"{design_prefix}_design.con",
            "-f", f"{design_prefix}_design.fts",
            "-m", str(mask_file),
            "-n", str(nperm),
            "-T", "-V",
        ], check=False)

    # Summary index
    with open(out_dir / "fc_seed_continuous_summary.json", "w") as f:
        json.dump({
            "n_subjects": diag["n_kept"],
            "ev_names": ev_names,
            "contrasts": [name for name, _ in contrasts],
            "vifs": diag["vifs_primary"],
            "score_ranges": diag["score_ranges"],
        }, f, indent=2, default=str)
    print(f"\n  Results: {out_dir}/")


def _write_fsl_design(design, ev_names, group_labels, prefix, kept):
    """Write FSL VEST .mat/.con/.fts files from a numpy design matrix."""
    n, n_ev = design.shape
    with open(f"{prefix}_design.mat", "w") as f:
        f.write(f"/NumWaves\t{n_ev}\n/NumPoints\t{n}\n")
        f.write("/PPheights\t" + "\t".join(
            f"{np.max(np.abs(design[:,i])):.6f}" for i in range(n_ev)) + "\n")
        f.write("\n/Matrix\n")
        for row in design:
            f.write("\t".join(f"{v:.6f}" for v in row) + "\n")

    contrasts = group_contrast_vectors(ev_names, group_labels)
    with open(f"{prefix}_design.con", "w") as f:
        f.write(f"/NumWaves\t{n_ev}\n/NumContrasts\t{len(contrasts)}\n")
        for i, (name, _) in enumerate(contrasts, 1):
            f.write(f"/ContrastName{i}\t{name}\n")
        f.write("\n/Matrix\n")
        for _, c in contrasts:
            f.write("\t".join(f"{v:.6f}" for v in c) + "\n")

    with open(f"{prefix}_design.fts", "w") as f:
        f.write(f"/NumWaves\t{len(contrasts)}\n/NumContrasts\t1\n")
        f.write("\n/Matrix\n")
        frow = np.zeros(len(contrasts))
        used_pairs = set()
        for k, (name, _) in enumerate(contrasts):
            pair = tuple(sorted(name.split("_gt_")))
            if pair not in used_pairs:
                frow[k] = 1
                used_pairs.add(pair)
        f.write("\t".join(f"{v:.0f}" for v in frow) + "\n")

    with open(f"{prefix}_subjects.txt", "w") as f:
        for sid, grp in kept:
            f.write(f"sub-{sid}\t{grp}\n")


# ======================================================================
# 3. Structural connectome analysis
# ======================================================================
def run_structural_connectome(subjects, cov_file, nperm,
                              control_depression=False, t_thresh=3.1):
    dep_label = " (+ depression)" if control_depression else ""
    print(f"\n=== Structural Connectome Group Analysis{dep_label} ===")
    suffix = "_depctrl" if control_depression else ""
    if t_thresh != 3.1:
        suffix += f"_t{t_thresh:.1f}".replace(".", "p")
    out_dir = GROUP_DIR / f"structural_connectome{suffix}"
    out_dir.mkdir(parents=True, exist_ok=True)

    from scipy import stats

    # Step 1: verify every subject's tractography output was built with the
    # SAME parcellation atlas. Subjects with a different atlas tag are
    # excluded with an error — refusing to silently mix atlases.
    data = {}
    atlases_seen = {}
    excluded_atlas = []
    for subj, grp in subjects:
        d = DERIV_DIR / "tractography" / f"sub-{subj}" / SESSION
        atlas_file = d / "atlas.txt"
        if atlas_file.is_file():
            atlas_name = atlas_file.read_text().strip()
        else:
            # Legacy output without atlas tag: assume DK to keep backward
            # compatibility with subjects produced by the old pipeline, but
            # record this so the user can see it in the summary.
            atlas_name = "desikan_killiany_legacy"
        atlases_seen[subj] = atlas_name

    # Determine the dominant atlas and exclude anything else
    from collections import Counter
    atlas_counts = Counter(atlases_seen.values())
    if len(atlas_counts) > 1:
        print(f"  WARNING: multiple atlases detected in tractography outputs:")
        for a, n in atlas_counts.most_common():
            print(f"    {a}: {n}")
        dominant = atlas_counts.most_common(1)[0][0]
        print(f"  Keeping only subjects with atlas = '{dominant}'")
    else:
        dominant = next(iter(atlas_counts)) if atlas_counts else None

    if dominant is None:
        print("  ERROR: no subjects with structural connectomes")
        return

    for subj, grp in subjects:
        if atlases_seen.get(subj) != dominant:
            excluded_atlas.append(subj)
            continue
        d = DERIV_DIR / "tractography" / f"sub-{subj}" / SESSION
        subj_data = {"group": grp, "atlas": dominant}
        for metric in ["sift2", "count", "fa"]:
            f = d / f"connectome_{metric}.csv"
            if f.is_file():
                try:
                    m = np.loadtxt(f, delimiter=",")
                    subj_data[metric] = m
                except Exception:
                    pass
        if len(subj_data) > 2:  # group + atlas + at least one metric
            data[subj] = subj_data

    if excluded_atlas:
        print(f"  Excluded {len(excluded_atlas)} subjects with non-dominant atlas: "
              f"{excluded_atlas[:5]}...")

    if len(data) < 10:
        print(f"  Only {len(data)} subjects with {dominant} structural connectomes")
        return

    print(f"  Atlas: {dominant} — {len(data)} subjects")

    subjects_ok = [(s, d["group"]) for s, d in data.items()]
    extra_covs = ["depression_factor"] if control_depression else None
    try:
        kept, design, ev_names, group_labels = load_covariates_strict(
            Path(cov_file), subjects_ok, cov_list=extra_covs)
    except Exception as e:
        print(f"  ERROR: {e}")
        return

    if control_depression:
        print(f"  Depression factor included as nuisance covariate")
    kept_set = {s for s, _ in kept}
    cov_cols = [i for i, n in enumerate(ev_names) if not n.startswith("group_")]

    for metric in ["sift2", "count", "fa"]:
        metric_data = [(s, data[s]["group"], data[s][metric])
                       for s in kept_set if metric in data[s]]
        if len(metric_data) < 10:
            continue

        # Within a single atlas, connectome shape MUST be identical across
        # subjects. A shape mismatch here indicates something went wrong
        # (partial FreeSurfer recon, corrupt CSV). Refuse to silently drop.
        shapes = {tuple(m.shape) for _, _, m in metric_data}
        if len(shapes) > 1:
            print(f"  ERROR: inconsistent connectome shapes within atlas for "
                  f"metric={metric}: {shapes}", file=sys.stderr)
            print(f"  Same atlas must yield identical matrix dimensions. "
                  f"Investigate per-subject tractography outputs.", file=sys.stderr)
            return
        ref_shape = metric_data[0][2].shape
        order = {s: i for i, (s, _) in enumerate(kept)}
        metric_data.sort(key=lambda x: order[x[0]])

        n_rois = ref_shape[0]
        mats = np.stack([m for _, _, m in metric_data])
        if metric in ("count", "sift2"):
            mats = np.log1p(mats)

        # Align design to metric_data order
        metric_subjs = [s for s, _, _ in metric_data]
        idx_in_kept = [order[s] for s in metric_subjs]
        design_metric = design[idx_in_kept]

        # Save group means
        for g in group_labels:
            idx = [i for i, (_, gg, _) in enumerate(metric_data) if gg == g]
            if idx:
                np.savetxt(out_dir / f"mean_sc_{metric}_{g}.csv",
                           mats[idx].mean(axis=0), delimiter=",", fmt="%.6f")

        # NBS per pairwise contrast
        print(f"\n  Metric: {metric} ({len(metric_data)} subjects, {n_rois} ROIs)")
        contrasts = group_contrast_vectors(ev_names, group_labels)
        for cname, cvec in contrasts:
            a, b = cname.split("_gt_")
            if a > b:
                continue
            res = run_nbs(mats, design_metric, cvec, cov_cols,
                          t_thresh=t_thresh, nperm=nperm)
            if res is None:
                continue
            n_sig_pos = sum(1 for _, _, p, _ in res["obs_components_pos"] if p < 0.05)
            n_sig_neg = sum(1 for _, _, p, _ in res["obs_components_neg"] if p < 0.05)
            print(f"    {cname}: {n_sig_pos} pos / {n_sig_neg} neg significant components (p<0.05)")
            with open(out_dir / f"nbs_{metric}_{cname}.json", "w") as f:
                json.dump({
                    "obs_components_pos": [
                        {"n_edges": s, "nodes": n, "p": p, "edges": e}
                        for s, n, p, e in res["obs_components_pos"]
                    ],
                    "obs_components_neg": [
                        {"n_edges": s, "nodes": n, "p": p, "edges": e}
                        for s, n, p, e in res["obs_components_neg"]
                    ],
                    "nperm": res["nperm"],
                    "t_thresh": res["t_thresh"],
                }, f, indent=2)


# ======================================================================
# 4. Structural connectome — CONTINUOUS (dimensional) analysis
# ======================================================================
def build_structural_continuous_design(cov_file: Path,
                                       subjects: list[tuple[str, str]],
                                       include_depression: bool = True,
                                       extended_covs: list[str] | None = None,
                                       quip_binary: bool = False):
    """Build a dimensional design matrix for structural connectome analysis.

    Same as build_continuous_design but uses dwi_motion instead of mean_FD
    (diffusion subjects don't have functional motion estimates).

    Parameters
    ----------
    extended_covs : list[str] | None
        Subset of EXTENDED_COVS_VALID to add beyond the lean model.
        Subjects missing any requested covariate are dropped.
    quip_binary : bool
        If True, dichotomise QUIP to presence/absence (>0 vs =0).
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

    required = ["apathy_composite", "quip", "age", "sex", "site", "TIV",
                "dwi_motion"]
    if include_depression:
        required.append("depression_factor")
    required.extend(extended_covs)

    missing_cols = [c for c in required if c not in fieldnames]
    if missing_cols:
        raise RuntimeError(
            f"covariates file missing columns required for structural "
            f"continuous analysis: {missing_cols}\n  header: {fieldnames}")

    # ── Filter subjects ──────────────────────────────────────────────
    # Broadened missing-value tokens (ported from build_continuous_design):
    # catches "NaN"/"nan" that pandas writes, preventing DLASCL errors.
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
            f"only {len(kept)} subjects have all regressors "
            f"populated (need >=10). First few dropped:\n  " +
            "\n  ".join(f"sub-{s}: {r}" for s, r in dropped[:15]))

    n = len(kept)

    def _col(name):
        return np.array([float(rows[sid][name]) for sid, _ in kept],
                        dtype=np.float64)

    # ── Continuous regressors (raw, then demean) ─────────────────────
    apathy = _col("apathy_composite")
    quip = _col("quip")
    if quip_binary:
        n_nonzero = int((quip > 0).sum())
        quip = (quip > 0).astype(np.float64)
        print(f"  [quip_binary] QUIP dichotomised to presence/absence: "
              f"{n_nonzero}/{len(quip)} subjects = 1 (QUIP>0), "
              f"{len(quip)-n_nonzero} = 0",
              file=sys.stderr)
    age = _col("age")
    tiv = _col("TIV")
    dwi_mot = _col("dwi_motion")
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
    dwi_mot_c = dwi_mot - dwi_mot.mean()
    dep_c = (dep - dep.mean()) if dep is not None else None
    disease_dur_c = (disease_dur - disease_dur.mean()) if disease_dur is not None else None
    updrs_c = (updrs - updrs.mean()) if updrs is not None else None
    ledd_c = (ledd - ledd.mean()) if ledd is not None else None
    moca_c = (moca - moca.mean()) if moca is not None else None
    apathy_x_quip = apathy_c * quip_c

    # ── Categorical dummies (sex, site) ──────────────────────────────
    def _dummies(name, drop_first=True):
        raw = [(rows[sid].get(name) or "").strip() for sid, _ in kept]
        levels = sorted(set(v for v in raw if v))
        if len(levels) <= 1:
            return [], np.zeros((n, 0))
        keep_levels = levels[1:] if drop_first else levels
        ref = levels[0] if drop_first else None
        cols, labels = [], []
        for lv in keep_levels:
            c = np.array([1.0 if r == lv else 0.0 for r in raw])
            c = c - c.mean()
            cols.append(c)
            labels.append(f"{name}_{lv}" + (f"_vs_{ref}" if ref else ""))
        return labels, np.column_stack(cols) if cols else np.zeros((n, 0))

    sex_labels, sex_mat = _dummies("sex")
    site_labels, site_mat = _dummies("site")

    # ── Assemble design matrix ───────────────────────────────────────
    ev_names: list[str] = []
    ev_cols: list[np.ndarray] = []

    ev_names.append("apathy_composite_c"); ev_cols.append(apathy_c)
    ev_names.append("quip_c");             ev_cols.append(quip_c)
    ev_names.append("apathy_x_quip");      ev_cols.append(apathy_x_quip)
    if dep_c is not None:
        ev_names.append("depression_factor_c"); ev_cols.append(dep_c)
    ev_names.append("age_c");              ev_cols.append(age_c)
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
    ev_names.append("TIV_c");             ev_cols.append(tiv_c)
    ev_names.append("dwi_motion_c");      ev_cols.append(dwi_mot_c)

    design = np.column_stack(ev_cols)
    n_ev = design.shape[1]

    # ── Sanity: no NaN/Inf allowed in the design ─────────────────────
    if not np.all(np.isfinite(design)):
        bad_cols = [ev_names[i] for i in range(n_ev)
                    if not np.all(np.isfinite(design[:, i]))]
        n_bad_rows = int((~np.isfinite(design).all(axis=1)).sum())
        raise RuntimeError(
            f"design matrix contains NaN/Inf in columns: {bad_cols} "
            f"({n_bad_rows} subjects affected). "
            "Upstream filter likely failed to recognise a missing-value "
            "token in the covariates file.")

    # ── Contrasts (match FC naming for cross-modality comparability) ──
    contrasts = []
    pretty = {"apathy_composite_c": "main_apathy",
              "quip_c": "main_quip",
              "apathy_x_quip": "interaction_apathyxquip"}
    for name in ["apathy_composite_c", "quip_c", "apathy_x_quip"]:
        cvec = np.zeros(len(ev_names))
        cvec[ev_names.index(name)] = 1.0
        contrasts.append((pretty[name], cvec))

    # Differential contrasts (specificity / double dissociation)
    # pos tail = apathy-specific, neg tail = quip-specific
    diff_vec = np.zeros(len(ev_names))
    diff_vec[ev_names.index("apathy_composite_c")] = 1.0
    diff_vec[ev_names.index("quip_c")] = -1.0
    contrasts.append(("diff_apathy_gt_quip", diff_vec))

    # Shared effect (common vulnerability)
    # pos tail = both predict increased SC, neg tail = both predict reduced SC
    shared_vec = np.zeros(len(ev_names))
    shared_vec[ev_names.index("apathy_composite_c")] = 1.0
    shared_vec[ev_names.index("quip_c")] = 1.0
    contrasts.append(("shared_pos", shared_vec))

    return kept, design, ev_names, contrasts


def run_structural_continuous(subjects, cov_file, nperm,
                              control_depression=True, t_thresh=3.1,
                              extended_covs=None, quip_binary=False):
    """Dimensional NBS on structural connectomes: apathy × ICB interaction.

    Uses dwi_motion instead of mean_FD (diffusion subjects).
    Depression_factor is included by default as a nuisance covariate.
    """
    dep_label = " (+ depression)" if control_depression else ""
    print(f"\n=== Structural Connectome Dimensional Analysis{dep_label} ===")
    suffix = "_depctrl" if control_depression else ""
    if t_thresh != 3.1:
        suffix += f"_t{t_thresh:.1f}".replace(".", "p")
    out_dir = GROUP_DIR / f"structural_connectome_continuous{suffix}"
    out_dir.mkdir(parents=True, exist_ok=True)

    # Load connectome data (same as categorical)
    data = {}
    atlases_seen = {}
    for subj, grp in subjects:
        d = DERIV_DIR / "tractography" / f"sub-{subj}" / SESSION
        atlas_file = d / "atlas.txt"
        atlases_seen[subj] = (atlas_file.read_text().strip()
                               if atlas_file.is_file()
                               else "desikan_killiany_legacy")

    from collections import Counter
    atlas_counts = Counter(atlases_seen.values())
    dominant = atlas_counts.most_common(1)[0][0] if atlas_counts else None
    if dominant is None:
        print("  ERROR: no subjects with structural connectomes")
        return

    for subj, grp in subjects:
        if atlases_seen.get(subj) != dominant:
            continue
        d = DERIV_DIR / "tractography" / f"sub-{subj}" / SESSION
        subj_data = {"group": grp}
        for metric in ["sift2", "count", "fa"]:
            f = d / f"connectome_{metric}.csv"
            if f.is_file():
                try:
                    subj_data[metric] = np.loadtxt(f, delimiter=",")
                except Exception:
                    pass
        if len(subj_data) > 1:
            data[subj] = subj_data

    print(f"  Atlas: {dominant} — {len(data)} subjects with connectomes")

    subjs_only = [(s, d["group"]) for s, d in data.items()]
    try:
        kept, design, ev_names, contrasts = build_structural_continuous_design(
            Path(cov_file), subjs_only,
            include_depression=control_depression,
            extended_covs=extended_covs,
            quip_binary=quip_binary)
    except Exception as e:
        print(f"  ERROR: {e}")
        return

    print(f"  Subjects with complete regressors: {len(kept)}")
    kept_set = {s for s, _ in kept}
    order = {s: i for i, (s, _) in enumerate(kept)}

    for metric in ["sift2", "count", "fa"]:
        metric_data = [(s, data[s]["group"], data[s][metric])
                       for s in kept_set if metric in data[s]]
        if len(metric_data) < 10:
            continue

        shapes = {tuple(m.shape) for _, _, m in metric_data}
        if len(shapes) > 1:
            print(f"  ERROR: inconsistent shapes for {metric}: {shapes}")
            return

        metric_data.sort(key=lambda x: order[x[0]])
        n_rois = metric_data[0][2].shape[0]
        mats = np.stack([m for _, _, m in metric_data])
        if metric in ("count", "sift2"):
            mats = np.log1p(mats)

        metric_subjs = [s for s, _, _ in metric_data]
        idx_in_kept = [order[s] for s in metric_subjs]
        design_metric = design[idx_in_kept]

        print(f"\n  Metric: {metric} ({len(metric_data)} subjects, {n_rois} ROIs)")
        for cname, cvec in contrasts:
            print(f"    NBS: {cname} ({nperm} perms)")
            res = _nbs_continuous(mats, design_metric, cvec,
                                  t_thresh=t_thresh, nperm=nperm, seed=42)
            for label, key in [("pos", "obs_components_pos"),
                               ("neg", "obs_components_neg")]:
                comps = res[key][:5]
                if comps:
                    for size, nodes, p, edges in comps:
                        star = " *" if p < 0.05 else ""
                        print(f"      {label}: {size} edges, {len(nodes)} nodes, "
                              f"p={p:.4f}{star}")
            with open(out_dir / f"nbs_{metric}_{cname}.json", "w") as f:
                json.dump({
                    "contrast": cname, "metric": metric,
                    "nperm": res["nperm"], "t_thresh": res["t_thresh"],
                    "null_max_pos_p95": res["null_max_pos_p95"],
                    "null_max_neg_p95": res["null_max_neg_p95"],
                    "obs_components_pos": [
                        {"n_edges": s, "nodes": n, "p": p, "edges": e}
                        for s, n, p, e in res["obs_components_pos"]],
                    "obs_components_neg": [
                        {"n_edges": s, "nodes": n, "p": p, "edges": e}
                        for s, n, p, e in res["obs_components_neg"]],
                }, f, indent=2)

    # Save design matrix and subject list for QA/reproducibility
    np.savetxt(out_dir / "design_matrix.tsv", design,
               delimiter="\t", header="\t".join(ev_names), comments="")
    with open(out_dir / "design_subjects.tsv", "w") as f:
        f.write("subject_id\tgroup\n")
        for sid, grp in kept:
            f.write(f"sub-{sid}\t{grp}\n")


# ======================================================================
# Main
# ======================================================================
def main():
    # Late-rebind: allow overriding FMRIPREP_ROOT from the CLI.
    # Must be declared before any read of the global inside this function.
    global FMRIPREP_ROOT

    parser = argparse.ArgumentParser(description="PPMI Group FC & Connectomics")
    parser.add_argument("--analysis", default="all",
                        choices=["all", "connectome", "seed", "structural",
                                 "connectome_continuous", "seed_continuous",
                                 "structural_continuous", "continuous"])
    parser.add_argument("--control-depression", action="store_true",
                        help="Include depression_factor as nuisance covariate")
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
                        help="Root dir that contains fmriprep_<subj>/sub-<subj>/ "
                             "(ephemeral) or fmriprep/sub-<subj>/ (home)")
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
    tract_subjects = find_tractography_subjects()
    print(f"\nSubjects with XCP-D rest:    {len(xcp_subjects)}")
    print(f"Subjects with tractography:  {len(tract_subjects)}")

    if args.analysis in ("all", "connectome"):
        run_connectome_analysis(xcp_subjects, args.covariates,
                                 args.nperm, args.atlas,
                                 t_thresh=args.t_thresh)

    if args.analysis in ("all", "seed"):
        run_seed_analysis(xcp_subjects, args.covariates, args.nperm)

    if args.analysis in ("all", "structural"):
        run_structural_connectome(tract_subjects, args.covariates, args.nperm,
                                  control_depression=args.control_depression,
                                  t_thresh=args.t_thresh)

    if args.analysis in ("all", "structural_continuous"):
        run_structural_continuous(tract_subjects, args.covariates, args.nperm,
                                  control_depression=args.control_depression,
                                  t_thresh=args.t_thresh,
                                  extended_covs=extended_covs,
                                  quip_binary=args.quip_binary)

    # Continuous (dimensional) analyses — apathy × ICB interaction
    if args.analysis in ("continuous", "connectome_continuous"):
        run_connectome_continuous(xcp_subjects, args.covariates,
                                   args.nperm, args.atlas,
                                   extended_covs=extended_covs,
                                   quip_binary=args.quip_binary,
                                   t_thresh=args.t_thresh)
    if args.analysis in ("continuous", "seed_continuous"):
        run_seed_continuous(xcp_subjects, args.covariates, args.nperm,
                            extended_covs=extended_covs,
                            quip_binary=args.quip_binary)

    print("\n============================================")
    print("Analysis complete.")
    print(f"Results in: {GROUP_DIR}/")


if __name__ == "__main__":
    main()

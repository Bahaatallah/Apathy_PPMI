# Role: Covariate and provenance loader imported by nbs_group_analysis.py at start-up; the connectome_continuous
#     analysis builds its design with nbs_group_analysis.build_continuous_design and does not call the
#     functions of this module.
# Manuscript item: none directly (dependency of nbs_group_analysis.py).
"""
group_covariates.py — Shared covariate + provenance loader

Single source of truth for covariate handling across all group analyses:
  - 07_group_config.sh::build_design_matrix (sh path, for TBSS/VBM voxelwise)
  - nbs_group_analysis.py                          (FC connectome + seed)
  - 11_vbm_group.sh::run_surface            (FreeSurfer surface GLM)
  - 12_morphometry_group.py                 (ROI morphometry stats)

Ensures:
  * Required covariate columns are present and non-missing
  * Processing provenance CSV (if present) is automatically merged as
    categorical nuisance regressors
  * Identical REQUIRED_COVS / HARMONISATION_COLS / policy in every path

Keep this file in sync with 07_group_config.sh — specifically the
REQUIRED_COV_COLS and HARMONISATION_COLS_DEFAULT lists.
"""

from __future__ import annotations

import csv
import os
import sys
from pathlib import Path

# ── Canonical lists (must match 07_group_config.sh) ──────────────────
# Modality-specific REQUIRED column policy. Selected by REQUIRED_MODALITY
# env var (default "structural" — most permissive — see 07_group_config.sh).
#   structural: BASE only (age, sex, site, TIV)
#   diffusion:  BASE + dwi_motion
#   functional: BASE + mean_FD
_BASE_REQUIRED = ["age", "sex", "site", "TIV"]
_REQUIRED_BY_MODALITY = {
    "structural": _BASE_REQUIRED,
    "diffusion":  _BASE_REQUIRED + ["dwi_motion"],
    "functional": _BASE_REQUIRED + ["mean_FD"],
}
REQUIRED_MODALITY = os.environ.get("REQUIRED_MODALITY", "structural")
REQUIRED_COVS = _REQUIRED_BY_MODALITY.get(REQUIRED_MODALITY, _BASE_REQUIRED)

HARMONISATION_COLS_DEFAULT = [
    "structural_bias_correction",
    "structural_freesurfer",
    "structural_flair_registered",
    "diffusion_topup",
    "diffusion_denoising",
    "diffusion_gibbs_removal",
    "diffusion_dwi_bias_correction",
    "diffusion_eddy_binary",
    "diffusion_registration_path",
    "diffusion_pe_override",
    "diffusion_scanner_manufacturer",
    "functional_topup",
    "functional_ica_method",
    "functional_bold_space",
    "functional_scanner_manufacturer",
]

CATEGORICAL_COVS = {"sex", "site", "scanner_model", "handedness"}

# Optional-covariate sparsity threshold (percentage of cohort with missing
# values above which an optional covariate is dropped from the design)
MAX_SPARSE_PCT = 10.0


def _strip_id(sid: str) -> str:
    if sid is None:
        return ""
    return str(sid).replace("sub-", "").strip()


def load_covariates_dict(cov_file: str | os.PathLike) -> tuple[dict, list[str]]:
    """Load covariates.csv into a dict keyed by subject_id (no 'sub-' prefix).

    Returns (cov_rows, cov_col_names). Fails loudly if required columns are
    missing — caller is expected to pre-check existence of the file.
    """
    cov_rows: dict[str, dict] = {}
    with open(cov_file) as f:
        reader = csv.DictReader(f)
        cov_cols = [c for c in (reader.fieldnames or [])
                    if c not in ("subject_id", "group", "participant_id")]
        for row in reader:
            sid = _strip_id(row.get("subject_id", row.get("participant_id", "")))
            if sid:
                cov_rows[sid] = row

    missing_required = [c for c in REQUIRED_COVS if c not in cov_cols]
    if missing_required:
        raise ValueError(
            f"covariates.csv missing required columns: {missing_required}\n"
            f"  Available: {cov_cols}")

    return cov_rows, cov_cols


def merge_provenance(cov_rows: dict, cov_cols: list[str],
                      deriv_dir: str | os.PathLike | None = None,
                      prov_csv_override: str | None = None,
                      verbose: bool = True) -> tuple[dict, list[str], list[str]]:
    """Merge processing_provenance.csv into cov_rows (covariates.csv wins).

    Returns (cov_rows, cov_col_names, prov_cols_added).

    Provenance CSV is located via:
      1. prov_csv_override argument
      2. PROV_CSV environment variable
      3. <deriv_dir>/provenance/processing_provenance.csv
    """
    prov_csv = prov_csv_override or os.environ.get("PROV_CSV", "")
    if not prov_csv and deriv_dir:
        prov_csv = str(Path(deriv_dir) / "provenance" / "processing_provenance.csv")

    prov_cols_added: list[str] = []
    if not prov_csv or not os.path.isfile(prov_csv):
        if verbose:
            print(f"  Provenance CSV not found (searched {prov_csv or 'none'})",
                  file=sys.stderr)
        return cov_rows, cov_cols, prov_cols_added

    with open(prov_csv) as f:
        reader = csv.DictReader(f)
        file_cols = [c for c in (reader.fieldnames or [])
                     if c not in ("subject", "session", "bids_subdir",
                                  "subject_id", "participant_id", "group")]
        # Only include provenance columns matching the current modality
        # (e.g. functional analysis should not include structural_* or
        # diffusion_* provenance columns — they are often perfectly
        # collinear with site dummies, causing rank deficiency)
        modality = REQUIRED_MODALITY.lower()
        to_merge = [c for c in HARMONISATION_COLS_DEFAULT
                    if c in file_cols and c.startswith(modality + "_")]
        if verbose and to_merge:
            skipped = [c for c in HARMONISATION_COLS_DEFAULT
                       if c in file_cols and not c.startswith(modality + "_")]
            if skipped:
                print(f"  Provenance: skipping {len(skipped)} non-{modality} "
                      f"columns: {skipped}", file=sys.stderr)

        for row in reader:
            sid = _strip_id(row.get("subject") or row.get("subject_id", ""))
            if not sid:
                continue
            if sid not in cov_rows:
                cov_rows[sid] = {}
            for c in to_merge:
                val = row.get(c, "")
                # covariates.csv values take precedence (do not overwrite)
                existing = cov_rows[sid].get(c, "")
                if not existing or not existing.strip():
                    cov_rows[sid][c] = val

        prov_cols_added = to_merge

    for c in prov_cols_added:
        if c not in cov_cols:
            cov_cols.append(c)

    if verbose and prov_cols_added:
        print(f"  Provenance integration: merged {len(prov_cols_added)} columns "
              f"from {prov_csv}", file=sys.stderr)

    return cov_rows, cov_cols, prov_cols_added


def filter_by_required(subjects: list[tuple[str, str]],
                       cov_rows: dict) -> list[tuple[str, str]]:
    """Drop subjects missing any REQUIRED covariate (identical policy to
    07_group_config.sh::build_design_matrix)."""
    kept = []
    for sid, grp in subjects:
        sid = _strip_id(sid)
        if sid not in cov_rows:
            print(f"  DROP sub-{sid}: no covariates row", file=sys.stderr)
            continue
        missing = [c for c in REQUIRED_COVS
                   if not cov_rows[sid].get(c, "").strip()
                   or cov_rows[sid][c].strip().lower() == "na"]
        if missing:
            print(f"  DROP sub-{sid}: missing required {missing}", file=sys.stderr)
            continue
        kept.append((sid, grp))
    return kept


def build_nuisance_matrix(kept_subjects: list[tuple[str, str]],
                           cov_rows: dict,
                           cov_list: list[str],
                           prov_cols_added: list[str]):
    """Build a numpy (n, n_ev) nuisance matrix.

    Treats REQUIRED_COVS and prov_cols_added + user-specified cov_list
    uniformly, with:
      - Continuous: mean-impute then demean
      - Categorical: dummy-code with a dedicated "missing" level rather
        than silently mapping to the reference. Prevents systematic bias
        from treating missing data as the reference category.
      - Zero-variance or >MAX_SPARSE_PCT sparsity → skip
    """
    import numpy as np

    n = len(kept_subjects)
    all_categoricals = set(CATEGORICAL_COVS) | set(prov_cols_added)
    required_set = set(REQUIRED_COVS)

    evs: list[tuple[str, "np.ndarray"]] = []

    for cname in cov_list:
        raw = [cov_rows[sid].get(cname, "").strip() for sid, _ in kept_subjects]
        missing_mask = [v == "" or v.lower() == "na" for v in raw]
        pct_missing = 100.0 * sum(missing_mask) / n

        # Sparsity gate (REQUIRED never skipped here — filter_by_required
        # already guaranteed they are present)
        if cname not in required_set and pct_missing > MAX_SPARSE_PCT:
            print(f"  SKIP optional covariate '{cname}': "
                  f"{pct_missing:.1f}% missing (threshold {MAX_SPARSE_PCT}%)",
                  file=sys.stderr)
            continue

        if cname in all_categoricals:
            # Explicit "missing" level for empty/na values so they do NOT
            # get silently coded as the reference level.
            levels_raw = [v if (v and v.lower() != "na") else "__missing__"
                          for v in raw]
            levels = sorted(set(levels_raw))
            if len(levels) < 2:
                print(f"  SKIP {cname}: zero variance (1 level)",
                      file=sys.stderr)
                continue
            ref = levels[0]
            for lv in levels[1:]:
                col = np.array([1.0 if v == lv else 0.0 for v in levels_raw])
                if col.std() < 1e-10:
                    continue
                col = col - col.mean()
                evs.append((f"{cname}_{lv}_vs_{ref}", col))
        else:
            try:
                non_missing = [float(v) for v in raw
                               if v and v.lower() != "na"]
                if not non_missing:
                    print(f"  SKIP {cname}: no numeric values",
                          file=sys.stderr)
                    continue
                col_mean = sum(non_missing) / len(non_missing)
                vals = []
                for v in raw:
                    if v and v.lower() != "na":
                        vals.append(float(v))
                    else:
                        vals.append(col_mean)
                col = np.array(vals)
            except ValueError as e:
                raise ValueError(
                    f"non-numeric value in continuous covariate {cname}: {e}")
            if col.std() < 1e-10:
                print(f"  SKIP {cname}: zero variance", file=sys.stderr)
                continue
            col = col - col.mean()
            evs.append((cname, col))

    if not evs:
        return np.zeros((n, 0)), []

    ev_names = [e[0] for e in evs]
    mat = np.column_stack([e[1] for e in evs])
    return mat, ev_names


def load_covariates_with_provenance(
        cov_file: str | os.PathLike,
        subjects: list[tuple[str, str]],
        deriv_dir: str | os.PathLike | None = None,
        extra_cov_list: list[str] | None = None,
        prov_csv_override: str | None = None,
        verbose: bool = True):
    """High-level one-shot loader used by FC/surface/morphometry scripts.

    Returns:
        kept_subjects:   [(sid, grp), ...] after REQUIRED-covariate filter
        design:          numpy array (n, n_ev) — group one-hot + nuisance
        ev_names:        list[str]
        group_labels:    list[str]
        nuisance_cols:   list[int]  — indices of non-group EVs
    """
    import numpy as np

    cov_file = Path(cov_file)
    if not cov_file.is_file():
        raise FileNotFoundError(
            f"Covariates file not found: {cov_file}\n"
            f"  Required columns: {REQUIRED_COVS}")

    cov_rows, cov_cols = load_covariates_dict(cov_file)
    cov_rows, cov_cols, prov_added = merge_provenance(
        cov_rows, cov_cols,
        deriv_dir=deriv_dir,
        prov_csv_override=prov_csv_override,
        verbose=verbose)

    # Default covariate list: REQUIRED + scanning factors + provenance.
    # Clinical variables (updrs, ledd, moca, bdi, disease duration, ...)
    # are INTENTIONALLY NOT auto-added for phenotype-group analyses
    # because they are likely correlated with group membership itself
    # and would absorb the phenotype effect. Pass them via
    # extra_cov_list if you explicitly want them adjusted for.
    SCANNING_FACTORS = ["mean_FD", "dwi_motion", "scanner_model", "handedness"]
    cov_list = list(REQUIRED_COVS)
    for c in SCANNING_FACTORS:
        if c in cov_cols and c not in cov_list:
            cov_list.append(c)
    for c in prov_added:
        if c not in cov_list:
            cov_list.append(c)
    if extra_cov_list:
        for c in extra_cov_list:
            if c and c not in cov_list and c in cov_cols:
                cov_list.append(c)

    kept = filter_by_required(subjects, cov_rows)
    if not kept:
        raise RuntimeError("no subjects have complete required covariates")

    # Recompute group labels after filter (match 07_group_config.sh fix)
    group_labels = sorted({g for _, g in kept})
    if len(group_labels) < 2:
        raise RuntimeError(
            f"fewer than 2 groups survived covariate filtering: {group_labels}")

    n = len(kept)
    if verbose:
        print(f"  {n} subjects retained ({len(group_labels)} groups: {group_labels})",
              file=sys.stderr)

    # Build group one-hot + nuisance EVs
    group_mat = np.column_stack([
        np.array([1.0 if kept[i][1] == g else 0.0 for i in range(n)])
        for g in group_labels
    ])

    nuisance_mat, nuisance_names = build_nuisance_matrix(
        kept, cov_rows, cov_list, prov_added)

    ev_names = [f"group_{g}" for g in group_labels] + nuisance_names
    design = np.column_stack([group_mat, nuisance_mat]) \
             if nuisance_mat.shape[1] > 0 else group_mat

    n_groups = len(group_labels)

    # ── Rank check: drop linearly dependent nuisance columns ────────
    # Provenance dummies can be perfectly collinear with site dummies,
    # causing rank deficiency.  Identify and remove redundant columns
    # from the nuisance block (never from the group columns).
    rank = np.linalg.matrix_rank(design)
    n_ev = design.shape[1]
    if rank < n_ev:
        n_redundant = n_ev - rank
        if verbose:
            print(f"  RANK CHECK: design {n} × {n_ev}, rank {rank} "
                  f"— {n_redundant} redundant column(s)", file=sys.stderr)
        # QR with column pivoting: last (n_ev - rank) pivoted columns
        # are the linearly dependent ones
        from scipy.linalg import qr
        _, R, piv = qr(design, pivoting=True)
        diag_R = np.abs(np.diag(R))
        # Columns whose diagonal R value is near zero are redundant
        tol = max(design.shape) * np.max(diag_R) * np.finfo(float).eps
        redundant_pivoted = set(piv[i] for i in range(n_ev)
                                if i >= rank or diag_R[min(i, len(diag_R)-1)] < tol)
        # Only drop nuisance columns (indices >= n_groups), never group
        to_drop = sorted(c for c in redundant_pivoted if c >= n_groups)
        if to_drop:
            dropped_names = [ev_names[i] for i in to_drop]
            if verbose:
                print(f"  Dropping {len(to_drop)} redundant nuisance EVs: "
                      f"{dropped_names}", file=sys.stderr)
            keep = [i for i in range(n_ev) if i not in to_drop]
            design = design[:, keep]
            ev_names = [ev_names[i] for i in keep]
            # Verify
            new_rank = np.linalg.matrix_rank(design)
            if verbose:
                print(f"  After drop: {design.shape[1]} EVs, rank {new_rank}",
                      file=sys.stderr)

    nuisance_cols = list(range(n_groups, design.shape[1]))

    return kept, design, ev_names, group_labels, nuisance_cols

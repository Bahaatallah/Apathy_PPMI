#!/bin/bash
# -----------------------------------------------------------------------
# structural_preprocessing.sh — Structural MRI preprocessing
# Usage: bash structural_preprocessing.sh <subject_id> [bids_subdir]
#   e.g. bash structural_preprocessing.sh XXXXXX Apathy_BL
#
# Steps: reorient → crop → N4 → BET → FAST → FIRST → FLIRT → FNIRT → SIENAX
#
# Scientific rationale:
#   - reorient2std ensures consistent L/R orientation across scanners
#   - robustfov removes neck tissue that confuses BET
#   - N4BiasFieldCorrection (ANTs): superior to FAST -B for PD with enlarged ventricles
#   - BET -R: iterative brain extraction on bias-corrected image
#   - FAST: probabilistic tissue segmentation on brain-extracted (for BBR in func)
#   - FIRST: subcortical volumes on WHOLE-HEAD image (does internal registration)
#   - FLIRT 12-DOF on whole-head → FNIRT on whole-head: consistent coordinate space
#   - FNIRT uses non-brain-extracted T1 (needs skull intensity gradients)
#   - SIENAX: normalised brain volume (standard call, does its own BET)
#   - FreeSurfer: uses best bias-corrected whole-head image
# -----------------------------------------------------------------------
set -euo pipefail

# ── Environment setup ─────────────────────────────────────────────────
# Works on both local Mac and HPC (module system)
if command -v module &>/dev/null; then
    module load tools/prod 2>/dev/null || true
    module load ANTs/2.3.5-foss-2021a 2>/dev/null || true
    module load FreeSurfer/7.4.1-centos8_x86_64 2>/dev/null || true
    module load FSL/6.0.7.17 2>/dev/null || true
    # Fix 1: ANTs foss-2021a pollutes PYTHONPATH with old numpy → breaks FSL Python
    unset PYTHONPATH 2>/dev/null || true
    # Fix 2: FreeSurfer MCR ships old libstdc++ that breaks ANTs shared libs
    export LD_LIBRARY_PATH=$(echo "$LD_LIBRARY_PATH" | tr ':' '\n' | grep -v MCR | tr '\n' ':' | sed 's/:$//')
fi

export FSLDIR="${FSLDIR:-$HOME/fsl}"
export PATH="$FSLDIR/bin:$PATH"
source "$FSLDIR/etc/fslconf/fsl.sh" 2>/dev/null || true

# Conda fallback for local Mac if ANTs not in PATH
if ! command -v N4BiasFieldCorrection &>/dev/null; then
    CONDA_DIR="${CONDA_DIR:-$HOME/miniconda3}"
    if [ -f "$CONDA_DIR/etc/profile.d/conda.sh" ]; then
        source "$CONDA_DIR/etc/profile.d/conda.sh"
        conda activate ppmi_mri 2>/dev/null || true
    fi
fi

if [ $# -lt 1 ]; then echo "Usage: $0 <subject_id> [bids_subdir]"; exit 1; fi
SUBJECT="$1"
BIDS_SUBDIR="${2:-}"  # e.g. "Apathy_BL" — optional, auto-detected if not provided

# Session timepoint (default BL). For 24M follow-up set SESSION=ses-24M in env.
SESSION="${SESSION:-ses-BL}"

# ── Path setup (parameterised for HPC portability) ────────────────────
# Use symlink to avoid spaces in path (FSL cannot handle spaces)
BIDS_ROOT="${BIDS_ROOT:-$HOME/PPMI_link/MRI_Data/BIDS}"
DERIV_DIR="${DERIV_DIR:-$HOME/PPMI_link/MRI_Data/derivatives}"

# Auto-detect BIDS subdirectory if not provided
if [ -z "$BIDS_SUBDIR" ]; then
    for candidate in Apathy_BL Both_BL ICB_BL Neither_BL ALL_24M; do
        if [ -d "${BIDS_ROOT}/${candidate}/sub-${SUBJECT}" ]; then
            BIDS_SUBDIR="$candidate"
            break
        fi
    done
    if [ -z "$BIDS_SUBDIR" ]; then
        echo "ERROR: sub-${SUBJECT} not found in any BIDS subdirectory"
        exit 1
    fi
fi

SUB_DIR="${BIDS_ROOT}/${BIDS_SUBDIR}/sub-${SUBJECT}/${SESSION}"
OUT_DIR="${DERIV_DIR}/structural/sub-${SUBJECT}/${SESSION}"
QC_DIR="${OUT_DIR}/qc"
LOG_FILE="${OUT_DIR}/structural_preproc.log"

# ── Idempotency check ────────────────────────────────────────────────
if [ -f "${OUT_DIR}/T1w_brain.nii.gz" ] && [ -f "${OUT_DIR}/T1w2MNI_warp.nii.gz" ]; then
    echo "[$(date)] Structural preprocessing already complete for sub-${SUBJECT}, skipping."
    echo "  Delete ${OUT_DIR}/ to force re-run."
    exit 0
fi
mkdir -p "${OUT_DIR}" "${QC_DIR}"

# Start logging
exec > >(tee -a "$LOG_FILE") 2>&1
echo "[$(date)] Starting structural preprocessing for sub-${SUBJECT} (${BIDS_SUBDIR})"

# ── Step 1: Locate T1w ──────────────────────────────────────────────────
T1W=$(find "${SUB_DIR}/anat" -name "*_T1w.nii.gz" 2>/dev/null | head -1)
if [ -z "$T1W" ]; then echo "ERROR: No T1w found for sub-${SUBJECT}"; exit 1; fi
echo "  T1w: $T1W"

# ── Step 2: Reorient to standard (MNI) orientation ─────────────────────
T1_REORI="${OUT_DIR}/T1w_reori.nii.gz"
fslreorient2std "$T1W" "$T1_REORI"

# ── Step 3: Crop FOV to remove neck/lower head ─────────────────────────
T1_CROP="${OUT_DIR}/T1w_crop.nii.gz"
robustfov -i "$T1_REORI" -r "$T1_CROP"

# ── Step 3b: N4 bias field correction (ANTs — superior to FAST -B) ────
# N4 iterative correction handles the spatially varying bias field better
# than BET -B, especially for 3T PD scans with large ventricles.
# Falls back to BET -B if ANTs is not installed.
T1_FOR_BET="$T1_CROP"
if command -v N4BiasFieldCorrection &>/dev/null; then
    echo "[$(date)] Running N4 bias field correction (ANTs CLI) ..."
    T1_N4="${OUT_DIR}/T1w_N4.nii.gz"
    N4BiasFieldCorrection -d 3 -i "$T1_CROP" -o "$T1_N4" \
        -s 4 -b "[200]" -c "[50x50x50x50,0.0]" 2>/dev/null
    T1_FOR_BET="$T1_N4"
    echo "  N4 correction applied"
elif python3 -c "import ants" 2>/dev/null; then
    echo "[$(date)] Running N4 bias field correction (ANTsPy) ..."
    T1_N4="${OUT_DIR}/T1w_N4.nii.gz"
    python3 -c "
import ants
img = ants.image_read('$T1_CROP')
corrected = ants.n4_bias_field_correction(img, shrink_factor=4, convergence={'iters': [50,50,50,50], 'tol': 0.0})
ants.image_write(corrected, '$T1_N4')
print('  ANTsPy N4 correction applied')
" 2>/dev/null && T1_FOR_BET="$T1_N4" || echo "  ANTsPy N4 failed — BET -B will handle bias"
else
    echo "  N4BiasFieldCorrection not found — BET -B will handle bias"
fi

# ── Step 4: Brain extraction ───────────────────────────────────────────
# -f 0.3: slightly conservative for PD brains (may have enlarged ventricles)
# -R: robust brain centre estimation (iterative)
# -B: bias field correction + neck cleanup (fallback if no N4)
T1_BRAIN="${OUT_DIR}/T1w_brain"
if [ "$T1_FOR_BET" = "$T1_CROP" ]; then
    # No N4 — use BET -B for bias correction
    bet "$T1_CROP" "$T1_BRAIN" -f 0.3 -R -B -m
else
    # N4 already corrected — skip BET's inferior bias correction
    bet "$T1_FOR_BET" "$T1_BRAIN" -f 0.3 -R -m
fi

# QC: Generate brain extraction overlay for visual inspection
slicer "${T1_CROP}" "${T1_BRAIN}.nii.gz" -s 2 \
    -x 0.35 "${QC_DIR}/bet_x1.png" -x 0.45 "${QC_DIR}/bet_x2.png" -x 0.55 "${QC_DIR}/bet_x3.png" \
    -y 0.35 "${QC_DIR}/bet_y1.png" -y 0.45 "${QC_DIR}/bet_y2.png" -y 0.55 "${QC_DIR}/bet_y3.png" \
    -z 0.35 "${QC_DIR}/bet_z1.png" -z 0.45 "${QC_DIR}/bet_z2.png" -z 0.55 "${QC_DIR}/bet_z3.png" 2>/dev/null || true
pngappend "${QC_DIR}/bet_x1.png" + "${QC_DIR}/bet_x2.png" + "${QC_DIR}/bet_x3.png" \
          + "${QC_DIR}/bet_y1.png" + "${QC_DIR}/bet_y2.png" + "${QC_DIR}/bet_y3.png" \
          + "${QC_DIR}/bet_z1.png" + "${QC_DIR}/bet_z2.png" + "${QC_DIR}/bet_z3.png" \
          "${QC_DIR}/bet_overlay.png" 2>/dev/null || true
echo "  QC: Check ${QC_DIR}/bet_overlay.png for brain extraction quality"

# ── Step 5: Tissue segmentation (FAST) ─────────────────────────────────
# -t 1: T1-weighted input; -n 3: 3 tissue classes (CSF, GM, WM)
# Outputs: T1w_pve_0 (CSF), T1w_pve_1 (GM), T1w_pve_2 (WM)
FAST_DIR="${OUT_DIR}/fast"; mkdir -p "$FAST_DIR"
fast -t 1 -n 3 -o "${FAST_DIR}/T1w" "${T1_BRAIN}.nii.gz"

# ── Step 6: Subcortical segmentation (FIRST) ──────────────────────────
# Key for PD: caudate, putamen, pallidum, thalamus, hippocampus, amygdala
# CRITICAL: FIRST must use the WHOLE-HEAD image, not brain-extracted.
# FIRST does its own internal MNI registration and needs skull context.
# Using brain-extracted input causes boundary artifacts in subcortical structures.
#
# NOTE: Do NOT pass -b here. The -b flag tells FIRST the input is
# already brain-extracted, causing it to skip its internal BET step.
# When paired with whole-head input this confuses the shape models,
# producing shifted segmentations with extreme L-R asymmetry (QC showed
# L_Amyg=367 vs R_Amyg=2493 for the same subject). Without -b, FIRST
# runs its own BET internally and produces symmetric, well-fitted
# segmentations.
FIRST_DIR="${OUT_DIR}/first"; mkdir -p "$FIRST_DIR"
run_first_all -i "$T1_FOR_BET" -o "${FIRST_DIR}/T1w" 2>/dev/null || {
    echo "  WARNING: FIRST failed (may need more memory) — continuing without subcortical segmentation"
}

# ── Step 7: Linear registration to MNI (FLIRT) ────────────────────────
# CRITICAL: FLIRT affine must be estimated from the SAME image family as FNIRT.
# Since FNIRT uses the whole-head T1 (needs skull for warping), FLIRT must also
# use whole-head images to avoid coordinate-space mismatch in the affine.
# Use -cost normmi (normalised mutual information) for whole-head registration.
FLIRT_MAT="${OUT_DIR}/T1w2MNI_affine.mat"
flirt -in "$T1_FOR_BET" -ref "$FSLDIR/data/standard/MNI152_T1_2mm.nii.gz" \
      -out "${OUT_DIR}/T1w2MNI_linear.nii.gz" -omat "$FLIRT_MAT" \
      -cost normmi -dof 12 -searchrx -90 90 -searchry -90 90 -searchrz -90 90

# Save the affine mat with a standard name for other pipelines
cp "$FLIRT_MAT" "${OUT_DIR}/T1w2MNI.mat"

# ── Step 8: Non-linear registration to MNI (FNIRT) ────────────────────
# Uses the same whole-head image as FLIRT (consistent coordinate space).
# MUST use 2mm MNI template to match the T1_2_MNI152_2mm config file.
FNIRT_WARP="${OUT_DIR}/T1w2MNI_warp.nii.gz"
fnirt --in="$T1_FOR_BET" --ref="$FSLDIR/data/standard/MNI152_T1_2mm.nii.gz" \
      --aff="${FLIRT_MAT}" --cout="${FNIRT_WARP}" \
      --iout="${OUT_DIR}/T1w2MNI_nonlinear.nii.gz" --config=T1_2_MNI152_2mm

# Also compute the inverse warp (MNI → native) for bringing atlases to native space
invwarp --ref="$T1_FOR_BET" --warp="${FNIRT_WARP}" --out="${OUT_DIR}/MNI2T1w_warp.nii.gz"

# QC: Registration overlay
slicer "${OUT_DIR}/T1w2MNI_nonlinear.nii.gz" "$FSLDIR/data/standard/MNI152_T1_2mm.nii.gz" -s 2 \
    -x 0.35 "${QC_DIR}/reg_x1.png" -x 0.55 "${QC_DIR}/reg_x2.png" \
    -y 0.35 "${QC_DIR}/reg_y1.png" -y 0.55 "${QC_DIR}/reg_y2.png" \
    -z 0.35 "${QC_DIR}/reg_z1.png" -z 0.55 "${QC_DIR}/reg_z2.png" 2>/dev/null || true
pngappend "${QC_DIR}/reg_x1.png" + "${QC_DIR}/reg_x2.png" \
          + "${QC_DIR}/reg_y1.png" + "${QC_DIR}/reg_y2.png" \
          + "${QC_DIR}/reg_z1.png" + "${QC_DIR}/reg_z2.png" \
          "${QC_DIR}/reg_overlay.png" 2>/dev/null || true
echo "  QC: Check ${QC_DIR}/reg_overlay.png for registration quality"

# ── Step 9: SIENAX (brain volume estimation) ──────────────────────────
# Standard call: sienax does its own internal BET and registration.
# -r: regional volumes; -B "... -f 0.3 -R": pass BET options for consistency
# Input is the best whole-head image (N4-corrected if available).
sienax "$T1_FOR_BET" -o "${OUT_DIR}/sienax" -r -B "-f 0.3 -R" 2>/dev/null || {
    echo "  WARNING: SIENAX failed — continuing"
}

# ── Step 10: FreeSurfer recon-all ─────────────────────────────────────
# Use the best whole-head bias-corrected image (N4 > crop > raw).
# FreeSurfer prefers whole-head input for its own skull-stripping.
if command -v recon-all &>/dev/null; then
    echo "[$(date)] Running FreeSurfer recon-all ..."
    export SUBJECTS_DIR="${DERIV_DIR}/freesurfer"; mkdir -p "$SUBJECTS_DIR"
    if [ ! -f "${SUBJECTS_DIR}/sub-${SUBJECT}_${SESSION}/mri/aparc+aseg.mgz" ]; then
        # Use N4-corrected whole-head if available, otherwise cropped
        FS_INPUT="$T1_FOR_BET"
        echo "  FreeSurfer input: $FS_INPUT"
        recon-all -subjid "sub-${SUBJECT}_${SESSION}" -i "$FS_INPUT" -all -threads 4
    else
        echo "  FreeSurfer already complete for sub-${SUBJECT}"
    fi
else
    echo "  FreeSurfer not available — skipping recon-all"
fi

# ── Step 11: Register FLAIR to T1w if available ──────────────────────
# FLAIR→T1 registration must use WHOLE-HEAD images on both sides:
#   - FLAIR is inherently whole-head (no skull stripping before registration)
#   - T1 whole-head provides consistent tissue boundaries for cross-modal alignment
#   - normmi cost function is appropriate for cross-modal (FLAIR ≠ T1 contrast)
FLAIR=$(find "${SUB_DIR}/anat" -name "*_FLAIR.nii.gz" 2>/dev/null | head -1)
if [ -n "$FLAIR" ]; then
    echo "  Co-registering FLAIR → T1w (whole-head, normmi)"
    FLAIR_REORI="${OUT_DIR}/FLAIR_reori.nii.gz"
    fslreorient2std "$FLAIR" "$FLAIR_REORI"
    flirt -in "$FLAIR_REORI" -ref "$T1_FOR_BET" \
          -out "${OUT_DIR}/FLAIR_to_T1w.nii.gz" -omat "${OUT_DIR}/FLAIR2T1w.mat" \
          -cost normmi -dof 6
fi

# ── Write provenance JSON (for harmonisation tracking) ───────────────
BIAS_METHOD="BET_B"
[ -f "${OUT_DIR}/T1w_N4.nii.gz" ] && BIAS_METHOD="N4_CLI"
# Check if ANTsPy was used (N4 file exists but N4BiasFieldCorrection CLI wasn't available)
command -v N4BiasFieldCorrection &>/dev/null || { [ -f "${OUT_DIR}/T1w_N4.nii.gz" ] && BIAS_METHOD="N4_ANTsPy"; }
FS_DONE=false
[ -f "${DERIV_DIR}/freesurfer/sub-${SUBJECT}_${SESSION}/mri/aparc+aseg.mgz" ] && FS_DONE=true
FLAIR_DONE=false
[ -f "${OUT_DIR}/FLAIR_to_T1w.nii.gz" ] && FLAIR_DONE=true
# Detect scanner manufacturer from JSON sidecar
SCANNER="unknown"
T1_JSON="${T1W%.nii.gz}.json"
if [ -f "$T1_JSON" ]; then
    SCANNER=$(python3 -c "
import json
with open('$T1_JSON') as f: j = json.load(f)
print(j.get('Manufacturer', 'unknown'))
" 2>/dev/null || echo "unknown")
fi

python3 -c "
import json, datetime
prov = {
    'subject': '${SUBJECT}',
    'session': 'BL',
    'bids_subdir': '${BIDS_SUBDIR}',
    'pipeline': 'structural',
    'timestamp': datetime.datetime.now().isoformat(),
    'processing': {
        'bias_correction': '${BIAS_METHOD}',
        'freesurfer': ${FS_DONE},
        'flair_registered': ${FLAIR_DONE},
        'scanner_manufacturer': '${SCANNER}'
    }
}
with open('${OUT_DIR}/provenance.json', 'w') as f:
    json.dump(prov, f, indent=2)
" 2>/dev/null || true

# ── Cleanup: remove intermediates to save storage ────────────────────
# Keep: T1w_brain, T1w_brain_mask, T1w2MNI_warp, MNI2T1w_warp, T1w2MNI.mat,
#       fast/, first/, sienax/, FLAIR_to_T1w, provenance.json, qc/
# Delete: reoriented, cropped, N4 whole-head (reconstructable from raw BIDS)
if [ "${CLEANUP_INTERMEDIATES:-true}" = "true" ]; then
    echo "[$(date)] Cleaning up intermediate files ..."
    BYTES_BEFORE=$(du -sb "${OUT_DIR}" 2>/dev/null | cut -f1)
    rm -f "${OUT_DIR}/T1w_reori.nii.gz"
    rm -f "${OUT_DIR}/T1w_crop.nii.gz"
    rm -f "${OUT_DIR}/T1w_N4.nii.gz"
    rm -f "${OUT_DIR}/T1w2MNI_linear.nii.gz"       # FLIRT output (have the warp)
    rm -f "${OUT_DIR}/T1w2MNI_nonlinear.nii.gz"     # FNIRT output (have the warp)
    rm -f "${OUT_DIR}/T1w2MNI_affine.mat"            # duplicate of T1w2MNI.mat
    rm -f "${OUT_DIR}/FLAIR_reori.nii.gz"
    rm -f "${OUT_DIR}/bold2T1_init.nii.gz" "${OUT_DIR}/bold2T1_init.mat"
    BYTES_AFTER=$(du -sb "${OUT_DIR}" 2>/dev/null | cut -f1)
    SAVED=$(( (BYTES_BEFORE - BYTES_AFTER) / 1048576 ))
    echo "  Freed ~${SAVED} MB"
fi

echo "[$(date)] Structural preprocessing COMPLETE for sub-${SUBJECT}"

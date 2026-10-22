#!/bin/bash
# -----------------------------------------------------------------------
# structural_preprocessing.sh — Structural MRI preprocessing
# Usage: bash structural_preprocessing.sh <subject_id> [bids_subdir]
#   e.g. bash structural_preprocessing.sh XXXXXX Apathy_BL
#
# Steps: reorient → crop → N4 → FreeSurfer recon-all
#
# Scientific rationale:
#   - reorient2std ensures consistent L/R orientation across scanners
#   - robustfov removes neck tissue that confuses skull stripping
#   - N4BiasFieldCorrection (ANTs): superior to FAST -B for PD with enlarged ventricles
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
if [ -f "${DERIV_DIR}/freesurfer/sub-${SUBJECT}_${SESSION}/mri/aparc+aseg.mgz" ]; then
    echo "[$(date)] Structural preprocessing already complete for sub-${SUBJECT}, skipping."
    echo "  Delete ${DERIV_DIR}/freesurfer/sub-${SUBJECT}_${SESSION}/ to force re-run."
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
# Without ANTs, FreeSurfer's own bias correction is relied on.
T1_FOR_FS="$T1_CROP"
if command -v N4BiasFieldCorrection &>/dev/null; then
    echo "[$(date)] Running N4 bias field correction (ANTs CLI) ..."
    T1_N4="${OUT_DIR}/T1w_N4.nii.gz"
    N4BiasFieldCorrection -d 3 -i "$T1_CROP" -o "$T1_N4" \
        -s 4 -b "[200]" -c "[50x50x50x50,0.0]" 2>/dev/null
    T1_FOR_FS="$T1_N4"
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
" 2>/dev/null && T1_FOR_FS="$T1_N4" || echo "  ANTsPy N4 failed; FreeSurfer will handle bias"
else
    echo "  N4BiasFieldCorrection not found; FreeSurfer will handle bias"
fi

# ── Step 10: FreeSurfer recon-all ─────────────────────────────────────
# Use the best whole-head bias-corrected image (N4 > crop > raw).
# FreeSurfer prefers whole-head input for its own skull-stripping.
if command -v recon-all &>/dev/null; then
    echo "[$(date)] Running FreeSurfer recon-all ..."
    export SUBJECTS_DIR="${DERIV_DIR}/freesurfer"; mkdir -p "$SUBJECTS_DIR"
    if [ ! -f "${SUBJECTS_DIR}/sub-${SUBJECT}_${SESSION}/mri/aparc+aseg.mgz" ]; then
        # Use N4-corrected whole-head if available, otherwise cropped
        FS_INPUT="$T1_FOR_FS"
        echo "  FreeSurfer input: $FS_INPUT"
        recon-all -subjid "sub-${SUBJECT}_${SESSION}" -i "$FS_INPUT" -all -threads 4
    else
        echo "  FreeSurfer already complete for sub-${SUBJECT}"
    fi
else
    echo "  FreeSurfer not available — skipping recon-all"
fi

# ── Write provenance JSON (for harmonisation tracking) ───────────────
BIAS_METHOD="none"
[ -f "${OUT_DIR}/T1w_N4.nii.gz" ] && BIAS_METHOD="N4_CLI"
# Check if ANTsPy was used (N4 file exists but N4BiasFieldCorrection CLI wasn't available)
command -v N4BiasFieldCorrection &>/dev/null || { [ -f "${OUT_DIR}/T1w_N4.nii.gz" ] && BIAS_METHOD="N4_ANTsPy"; }
FS_DONE=false
[ -f "${DERIV_DIR}/freesurfer/sub-${SUBJECT}_${SESSION}/mri/aparc+aseg.mgz" ] && FS_DONE=true
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
        'freesurfer': '${FS_DONE}' == 'true',
        'scanner_manufacturer': '${SCANNER}'
    }
}
with open('${OUT_DIR}/provenance.json', 'w') as f:
    json.dump(prov, f, indent=2)
" 2>/dev/null || true

# ── Cleanup: remove intermediates to save storage ────────────────────
# Keep: FreeSurfer input (T1w_N4, or T1w_crop without N4), provenance.json
# Delete: T1w_reori, and T1w_crop when N4 was applied (reconstructable from raw BIDS)
if [ "${CLEANUP_INTERMEDIATES:-true}" = "true" ]; then
    echo "[$(date)] Cleaning up intermediate files ..."
    BYTES_BEFORE=$(du -sb "${OUT_DIR}" 2>/dev/null | cut -f1)
    rm -f "${OUT_DIR}/T1w_reori.nii.gz"
    [ "$T1_FOR_FS" = "$T1_CROP" ] || rm -f "${OUT_DIR}/T1w_crop.nii.gz"
    BYTES_AFTER=$(du -sb "${OUT_DIR}" 2>/dev/null | cut -f1)
    SAVED=$(( (BYTES_BEFORE - BYTES_AFTER) / 1048576 ))
    echo "  Freed ~${SAVED} MB"
fi

echo "[$(date)] Structural preprocessing COMPLETE for sub-${SUBJECT}"

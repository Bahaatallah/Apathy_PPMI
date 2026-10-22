#!/bin/bash
# -----------------------------------------------------------------------
# setup_hpc.sh — One-time HPC setup: pull containers, set up dirs
# Run this ONCE on the HPC login node before submitting any jobs.
#
# Usage: bash setup_hpc.sh
# -----------------------------------------------------------------------
set -euo pipefail

echo "============================================"
echo "PPMI MRI Pipeline — HPC Setup (Imperial CX3)"
echo "============================================"
echo "Date: $(date)"
echo ""

# ── Paths ────────────────────────────────────────────────────────────
# Storage strategy: Home (~325 GB free) for finals only
# Ephemeral (~11 TB) for containers, processing, fMRIPrep work
HOME_DIR="$HOME"
RDS_HOME="/rds/general/user/${USER}/home"
RDS_EPHEMERAL="/rds/general/user/${USER}/ephemeral"

# Containers on Ephemeral (~15 GB — can re-pull if needed)
CONTAINER_DIR="${RDS_EPHEMERAL}/PPMI/containers"
TEMPLATEFLOW_DIR="${RDS_HOME}/.cache/templateflow"
FREESURFER_LICENSE="${HOME_DIR}/freesurfer_license.txt"

# RDS project space
RDS_PROJECT="${RDS_HOME}/PPMI"
EPHEMERAL="${RDS_EPHEMERAL}/PPMI"

mkdir -p "$CONTAINER_DIR"
mkdir -p "$TEMPLATEFLOW_DIR"
mkdir -p "$RDS_PROJECT/MRI_Data/BIDS"
mkdir -p "$RDS_PROJECT/MRI_Data/derivatives"
mkdir -p "$RDS_PROJECT/MRI_Data/logs"
mkdir -p "$EPHEMERAL/work"       # fMRIPrep work dir (large, temporary)
mkdir -p "$EPHEMERAL/derivatives" # Processing intermediates

echo "[1/6] Directories created"
echo "  Containers:   $CONTAINER_DIR"
echo "  TemplateFlow: $TEMPLATEFLOW_DIR"
echo "  Project:      $RDS_PROJECT"
echo "  Ephemeral:    $EPHEMERAL"
echo ""

# ── Load modules ─────────────────────────────────────────────────────
echo "[2/6] Loading modules..."
module load Python/3.12.3-GCCcore-13.3.0 2>/dev/null || \
module load Python/3.11.3-GCCcore-12.3.0 2>/dev/null || \
    echo "WARNING: Could not load Python module"

# ── Pull Singularity containers ──────────────────────────────────────
echo "[3/6] Pulling Singularity containers (this takes 15-30 min)..."
echo "  If a container already exists, it will be skipped."
echo ""

# MRIQC — quality control (run before preprocessing)
MRIQC_SIF="$CONTAINER_DIR/mriqc-24.0.2.sif"
if [ ! -f "$MRIQC_SIF" ]; then
    echo "  Pulling MRIQC 24.0.2..."
    singularity pull "$MRIQC_SIF" docker://nipreps/mriqc:24.0.2
else
    echo "  MRIQC already exists: $MRIQC_SIF"
fi

# fMRIPrep — functional preprocessing (LTS release)
# NOTE: --use-aroma was removed in fMRIPrep 23.1. AROMA denoising is now
#       handled by XCP-D (--nuisance-regressors) or standalone fmripost-aroma.
FMRIPREP_SIF="$CONTAINER_DIR/fmriprep-25.2.5.sif"
if [ ! -f "$FMRIPREP_SIF" ]; then
    echo "  Pulling fMRIPrep 25.2.5 LTS..."
    singularity pull "$FMRIPREP_SIF" docker://nipreps/fmriprep:25.2.5
else
    echo "  fMRIPrep already exists: $FMRIPREP_SIF"
fi

# XCP-D — functional denoising and connectivity (post-fMRIPrep)
XCPD_SIF="$CONTAINER_DIR/xcp_d-26.0.1.sif"
if [ ! -f "$XCPD_SIF" ]; then
    echo "  Pulling XCP-D 26.0.1..."
    singularity pull "$XCPD_SIF" docker://pennlinc/xcp_d:26.0.1
else
    echo "  XCP-D already exists: $XCPD_SIF"
fi

echo ""

# ── TemplateFlow setup ───────────────────────────────────────────────
echo "[4/6] Setting up TemplateFlow templates..."
export TEMPLATEFLOW_HOME="$TEMPLATEFLOW_DIR"
pip install --user templateflow 2>/dev/null || true
python3 -c "
import templateflow.api as tflow
# Download templates needed by fMRIPrep and our pipeline
templates = ['MNI152NLin2009cAsym', 'MNI152NLin6Asym', 'OASIS30ANTs', 'fsLR', 'fsaverage']
for t in templates:
    try:
        tflow.get(t, raise_empty=True)
        print(f'  {t}: OK')
    except Exception as e:
        print(f'  {t}: downloading...')
        try:
            tflow.get(t)
            print(f'  {t}: OK')
        except Exception as e2:
            print(f'  {t}: FAILED ({e2})')
" 2>/dev/null || echo "  WARNING: TemplateFlow setup may be incomplete — fMRIPrep will download on first run"
echo ""

# ── FreeSurfer license ───────────────────────────────────────────────
echo "[5/6] FreeSurfer license check..."
if [ -f "$FREESURFER_LICENSE" ]; then
    echo "  License found: $FREESURFER_LICENSE"
elif [ -f "$HOME/.freesurfer_license" ]; then
    cp "$HOME/.freesurfer_license" "$FREESURFER_LICENSE"
    echo "  Copied from ~/.freesurfer_license"
else
    echo "  ⚠ NO FreeSurfer license found!"
    echo "  Get one FREE from: https://surfer.nmr.mgh.harvard.edu/registration.html"
    echo "  Save as: $FREESURFER_LICENSE"
    echo "  fMRIPrep and FreeSurfer will NOT run without it."
fi
echo ""

# ── Verify modules ───────────────────────────────────────────────────
echo "[6/6] Verifying available modules..."
echo ""
echo "  Module availability:"
for mod in "FSL/6.0.7.17" "FreeSurfer/7.4.1-centos8_x86_64" "ANTs/2.3.5-foss-2021a"; do
    if module avail "$mod" 2>&1 | grep -q "$mod"; then
        echo "    ✓ $mod"
    else
        echo "    ✗ $mod — not found (check with: module avail)"
    fi
done
echo ""

# ── Summary ──────────────────────────────────────────────────────────
echo "============================================"
echo "SETUP COMPLETE"
echo "============================================"
echo ""
echo "Container sizes:"
ls -lh "$CONTAINER_DIR"/*.sif 2>/dev/null || echo "  (none yet)"
echo ""
echo "Next steps:"
echo "  1. Ensure your FreeSurfer license is at: $FREESURFER_LICENSE"
echo "  2. Upload your BIDS data to: $RDS_PROJECT/MRI_Data/BIDS/"
echo "  3. Upload pipeline scripts to: $RDS_PROJECT/MRI_Data/scripts/"
echo "  4. Run: bash submit_pipeline.sh"
echo ""
echo "Storage layout:"
echo "  Home (persistent):      finals, BIDS, FreeSurfer, provenance"
echo "  Ephemeral (temporary):  containers, intermediates, fMRIPrep work"
echo ""
echo "Storage usage:"
echo "  Home:       $(du -sh "$HOME_DIR" 2>/dev/null | cut -f1)"
echo "  Containers: $(du -sh "$CONTAINER_DIR" 2>/dev/null | cut -f1) (on ephemeral)"
echo ""
echo "NOTE: Containers are on ephemeral storage (~15 GB). If ephemeral is"
echo "cleaned, re-run this script to re-pull them before submitting jobs."

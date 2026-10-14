#!/usr/bin/env python3
"""Shared configuration for the PPMI MRI pipeline."""

import os
from pathlib import Path

# ── Paths ──────────────────────────────────────────────────────────────────
# Data lives on OneDrive (Google Drive FUSE is unreliable for batch file I/O)
# Use symlink at ~/PPMI_link for FSL compatibility (no spaces in path)
PROJECT_ROOT = Path(os.environ.get("PPMI_IMAGING_ROOT", "."))
MRI_DATA     = PROJECT_ROOT / "MRI_Data"
BIDS_DIR     = MRI_DATA / "BIDS"
DERIV_DIR    = MRI_DATA / "derivatives"
SCRIPTS_DIR  = MRI_DATA / "scripts"
LOGS_DIR     = MRI_DATA / "logs"

# Symlink path (no spaces — required for FSL tools)
SYMLINK_ROOT = Path(os.path.expanduser("~/PPMI_link/MRI_Data"))
BIDS_SYMLINK = SYMLINK_ROOT / "BIDS"
DERIV_SYMLINK = SYMLINK_ROOT / "derivatives"

# Directories containing raw DICOM data (extracted or zipped) + idaSearch CSVs
DICOM_ZIP_DIRS = [
    MRI_DATA / "Apathy_BL",
    MRI_DATA / "Both_BL",
    MRI_DATA / "ICB_BL",
]

# BIDS subdirectories (subject groups)
BIDS_SUBDIRS = ["Apathy_BL", "Both_BL", "ICB_BL"]

# ── Tool paths ────────────────────────────────────────────────────────────
FSLDIR            = os.environ.get("FSLDIR", os.path.expanduser("~/fsl"))
CONDA_DIR         = Path(os.path.expanduser("~/miniconda3"))
CONDA_ENV         = "ppmi_mri"  # conda env with MRtrix3 + ANTs

# ── Processing parameters ─────────────────────────────────────────────────
# Structural
BET_FRAC          = 0.3
# FIX: Use 2mm for both FLIRT and FNIRT (must match T1_2_MNI152_2mm config)
MNI_REF           = "$FSLDIR/data/standard/MNI152_T1_2mm.nii.gz"
MNI_BRAIN         = "$FSLDIR/data/standard/MNI152_T1_2mm_brain.nii.gz"
FNIRT_CONFIG      = "T1_2_MNI152_2mm"
FREESURFER_VER    = "7.4.1"
FSL_VER           = "6.0.7"

# Diffusion
EDDY_REPOL        = True
EDDY_CNR_MAPS     = True
USE_TOPUP         = True
DWI_BET_FRAC      = 0.2   # less aggressive than structural (standard for diffusion b0)
DWI_MIN_DIM4      = 2     # reject DWI with dim4 <= 1 (corrupted 3D stacking)

# Functional
FUNC_TR           = 2.5     # seconds (read from JSON sidecar when available)
SMOOTH_FWHM       = 6.0     # mm
HP_CUTOFF_S       = 150.0   # high-pass cutoff period in seconds (standard with ICA-AROMA)
ICA_AROMA         = True
N_TRIM_VOLS       = 4       # initial volumes to discard
MIN_BOLD_VOLS     = 50      # reject BOLD with fewer volumes (likely calibration scan)

# HPC — Imperial CX3 (PBS)
HPC_SCHEDULER     = "pbs"
HPC_QUEUES = {
    "small":    "v1_small24",       # ≤24h, small memory
    "medium":   "v1_medium72",      # ≤72h, medium memory
    "large":    "v1_large72",       # ≤72h, large memory
    "gpu":      "v1_gpu72",         # ≤72h, GPU nodes
    "largemem": "v1_largemem72",    # ≤72h, high memory
}
HPC_DEFAULT_QUEUE = "v1_medium72"
HPC_CPUS          = 4
HPC_MEM           = "16gb"
HPC_GPU_MEM       = "16gb"
N_THREADS         = 4

# Containers (Singularity/Apptainer)
CONTAINER_DIR     = "$HOME/containers"
MRIQC_VERSION     = "24.0.2"
FMRIPREP_VERSION  = "25.2.5"   # LTS; --use-aroma removed since 23.1
XCPD_VERSION      = "26.0.1"   # calendar versioning since 26.0.0

# HPC Modules (Imperial CX3)
HPC_MODULES = {
    "fsl":        "FSL/6.0.7.17",
    "freesurfer": "FreeSurfer/7.4.1-centos8_x86_64",
    "ants":       "ANTs/2.3.5-foss-2021a",
    "mrtrix":     "MRtrix/3.0.4-foss-2022b",
    "afni":       "AFNI/24.0.02-foss-2023a",
    "python":     "Python/3.12.3-GCCcore-13.3.0",
}

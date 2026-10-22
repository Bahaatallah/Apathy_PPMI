#!/usr/bin/env python3
"""Shared configuration for the PPMI MRI pipeline (paths used by dicom_to_bids.py)."""

import os
from pathlib import Path

# Paths. PPMI_IMAGING_ROOT holds MRI_Data/ with the BIDS tree and the derivatives.
PROJECT_ROOT = Path(os.environ.get("PPMI_IMAGING_ROOT", "."))
MRI_DATA     = PROJECT_ROOT / "MRI_Data"
BIDS_DIR     = MRI_DATA / "BIDS"
DERIV_DIR    = MRI_DATA / "derivatives"

# Directories containing raw DICOM data (extracted or zipped) and the idaSearch CSVs
DICOM_ZIP_DIRS = [
    MRI_DATA / "Apathy_BL",
    MRI_DATA / "Both_BL",
    MRI_DATA / "ICB_BL",
]

# BIDS subdirectories (subject groups)
BIDS_SUBDIRS = ["Apathy_BL", "Both_BL", "ICB_BL"]

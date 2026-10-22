#!/usr/bin/env python3
"""
dicom_to_bids.py
===================
Convert raw PPMI DICOM downloads to BIDS format using dcm2niix.

Handles both zip archives AND pre-extracted DICOM directory trees.
Expected extracted structure:
    PPMI_2/SubjectID/SequenceName/Date/ImageID/*.dcm

Steps:
  1. Discover DICOM files (from extracted dirs or zip archives)
  2. Parse idaSearch CSV for subject/sequence metadata
  3. Run dcm2niix per subject+sequence to produce NIfTI + JSON sidecars
  4. Organise outputs into BIDS directory structure
  5. Generate participants.tsv and dataset_description.json

Requirements:
  dcm2niix must be on PATH (brew install dcm2niix / module load dcm2niix)

Usage:
  python dicom_to_bids.py                    # Process all subjects
  python dicom_to_bids.py --subject XXXXXX   # Single subject
  python dicom_to_bids.py --dry-run           # Preview only
"""

import argparse
import csv
import json
import logging
import os
import re
import shutil
import subprocess
import sys
import tempfile
import zipfile
from pathlib import Path

logging.basicConfig(level=logging.INFO, format="%(asctime)s  %(levelname)-8s  %(message)s")
log = logging.getLogger(__name__)

# Import paths from config
sys.path.insert(0, str(Path(__file__).resolve().parent))
from config import BIDS_DIR, DICOM_ZIP_DIRS, MRI_DATA

# ── Sequence classification ──────────────────────────────────────────────
# Maps PPMI sequence directory names to (BIDS datatype, BIDS suffix, phase-encoding)
# Order matters — first match wins. More specific patterns must come before generic ones.
SEQUENCE_MAP = [
    # Structural T1
    (r"(?i)(3D[_\s-]*T1|MPRAGE|FSPGR|sT1W_3D_TFE|T1.*weighted|T1.*3D.*Volumetric|SAG.*T1|T1.*MPRAGE|^T1$)",
     "anat", "T1w", None),
    # FLAIR (including _ND non-distortion-corrected variants)
    (r"(?i)(FLAIR|3D.*T2.*FLAIR|CUBE.*FLAIR|Brain_VIEW_FLAIR)",
     "anat", "FLAIR", None),
    # PD/T2 (including bare "T2" from some Philips/GE sites)
    (r"(?i)(^T2$|PD.*T2|Double_TSE|T2W_TSE|AX_PD|AX.*T2|AX.*DUAL_TSE|AX.*DE_TSE|AXIAL.*T2|sT2W.*TSE)",
     "anat", "T2w", None),
    # Neuromelanin / GRE-MT
    (r"(?i)(GRE.*M[TN]|NM.*MT|GRE.*NM|NM.?GRE|MT[C_])",
     "anat", "NM", None),
    # DTI with phase encoding — specific PE patterns before generic DTI
    (r"(?i)(DTI.*(L[_\s/-]+R)|DTI_+LR|2D_DTI_EPI_FAT_SHIFT_LEFT)",
     "dwi", "dwi", "LR"),
    (r"(?i)(DTI.*(R[_\s/-]+L)|DTI_+RL|2D_DTI_EPI_FAT_SHIFT_RIGHT)",
     "dwi", "dwi", "RL"),
    (r"(?i)DTI.*AP|DTI.*revB0.*AP",
     "dwi", "dwi", "AP"),
    (r"(?i)DTI.*PA|DTI.*B0.*PA",
     "dwi", "dwi", "PA"),
    # DTI generic — gated/non-gated/triggered variants, no PE direction specified
    (r"(?i)(^DTI$|DT_SSh_iso.*SENSE|Ax_DTI|^DTI_\d+dir$|DTI.*gated|DTI.*NON.?GATED|REPEAT.*DTI|TRIGGERED.*DTI)",
     "dwi", "dwi", None),
    # rsfMRI with phase encoding
    (r"(?i)(rsfMRI.*RL|RSFMRI.*RL|rsfMRI.*FS-L|RESTING.*RIGHT)",
     "func", "bold", "RL"),
    (r"(?i)(rsfMRI.*LR|RSFMRI.*LR|rsfMRI.*FS-R|RESTING.*LEFT|rsfMRI.*L_R)",
     "func", "bold", "LR"),
    (r"(?i)(rsfMRI.*AP)",
     "func", "bold", "AP"),
    (r"(?i)(rsfMRI.*PA)",
     "func", "bold", "PA"),
    # rsfMRI / BOLD generic
    (r"(?i)(rsfMRI|rsFMRI|BOLD.*REST|BOLD.*RS|ep2d.*rest|RESTING.*STATE|ep2d.*RESTING|ep2d_bold_rest|repeat.*BOLD.*RS|REPEAT.*fMRI)",
     "func", "bold", None),
    # Fieldmaps
    (r"(?i)B0.*rf.*Map",
     "fmap", "fieldmap", None),
]

# Skip these derived/scout sequences entirely
SKIP_PATTERNS = [
    r"(?i)^DTI_FA$",
    r"(?i)^DTI_ADC$",
    r"(?i)^DTI_TRACEW$",
    r"(?i)^(SURVEY|SAG|COR|TRA|CORONAL|TRANSVERSE)$",
    r"(?i)split_1",
]


def classify_sequence(seq_name: str):
    """Classify a PPMI sequence directory name.
    Returns (datatype, suffix, phase_encoding) or None to skip."""
    for pattern in SKIP_PATTERNS:
        if re.search(pattern, seq_name):
            return None
    for pattern, datatype, suffix, pe in SEQUENCE_MAP:
        if re.search(pattern, seq_name):
            return (datatype, suffix, pe)
    log.warning(f"  Unrecognised sequence: '{seq_name}' — skipping")
    return None


def discover_extracted_dicoms(zip_dir: Path, subject_filter: str = None) -> dict:
    """Discover DICOM directories from pre-extracted directory trees.
    Uses `find` for fast bulk discovery on Google Drive (avoids slow per-dir Python I/O).
    Expected: PPMI*/SubjectID/SequenceName/Date/ImageID/
    Returns {(subj, seq_name, img_id): dicom_directory_path}"""
    entries = {}
    ppmi_dirs = sorted(zip_dir.glob("PPMI*"))
    ppmi_dirs = [d for d in ppmi_dirs if d.is_dir()]

    if not ppmi_dirs:
        return entries

    for ppmi_dir in ppmi_dirs:
        log.info(f"Scanning extracted directory: {ppmi_dir.name}/ ...")

        # Use `find` at exactly depth 4 (SubjID/SeqName/Date/ImageID) to get
        # leaf DICOM directories in one system call — orders of magnitude faster
        # than Python iterdir() over Google Drive.
        try:
            result = subprocess.run(
                ["find", str(ppmi_dir), "-mindepth", "4", "-maxdepth", "4", "-type", "d"],
                capture_output=True, text=True, timeout=300
            )
            dir_paths = [p for p in result.stdout.strip().split("\n") if p]
        except (subprocess.TimeoutExpired, OSError) as e:
            log.warning(f"  find failed for {ppmi_dir.name}: {e}, falling back to Python walk")
            dir_paths = []
            for root, dirs, files in os.walk(ppmi_dir):
                depth = Path(root).relative_to(ppmi_dir).parts
                if len(depth) == 4:
                    dir_paths.append(root)
                    dirs.clear()  # don't descend further

        subj_count = 0
        seen_subjs = set()
        for dir_path in sorted(dir_paths):
            p = Path(dir_path)
            parts = p.relative_to(ppmi_dir).parts
            if len(parts) != 4:
                continue
            subj, seq_name, _date, img_id_raw = parts

            if subject_filter and subj != subject_filter:
                continue

            if subj not in seen_subjs:
                seen_subjs.add(subj)
                subj_count += 1

            m = re.search(r"I?(\d+)", img_id_raw)
            img_id = m.group(1) if m else img_id_raw

            key = (subj, seq_name, img_id)
            entries[key] = dir_path  # just the directory — dcm2niix takes dirs

        log.info(f"  Finished {ppmi_dir.name}/: {subj_count} subjects")

    log.info(f"Found {len(entries)} unique subject x sequence x run combinations (extracted)")
    return entries


def discover_zip_dicoms(zip_dir: Path, subject_filter: str = None) -> dict:
    """Discover DICOM files from zip archives.
    Returns {(subj, seq_name, img_id): (zip_path, [dcm_paths_in_zip])}"""
    entries = {}
    for zf_path in sorted(zip_dir.glob("*.zip")):
        log.info(f"Scanning zip: {zf_path.name}")
        with zipfile.ZipFile(zf_path, "r") as zf:
            for name in zf.namelist():
                if not name.lower().endswith(".dcm"):
                    continue
                parts = name.split("/")
                # Expected: PPMI/SubjectID/SequenceName/Date/ImageID/file.dcm
                if len(parts) < 5:
                    continue
                subj = parts[1]
                seq_name = parts[2]
                img_id_dir = parts[4] if len(parts) >= 6 else parts[3]
                m = re.search(r"I?(\d+)", img_id_dir)
                img_id = m.group(1) if m else img_id_dir

                if subject_filter and subj != subject_filter:
                    continue

                key = (subj, seq_name, img_id)
                if key not in entries:
                    entries[key] = (zf_path, [])
                entries[key][1].append(name)
    return entries


def convert_to_bids(zip_dir: Path, bids_dir: Path,
                    subject_filter: str = None, dry_run: bool = False):
    """Convert DICOMs (extracted or zipped) to BIDS with dcm2niix."""
    dcm2niix_bin = shutil.which("dcm2niix")
    if not dcm2niix_bin:
        # Check common pip install locations
        import site
        for d in site.getusersitepackages(), *site.getsitepackages():
            candidate = Path(d).parent.parent / "bin" / "dcm2niix"
            if candidate.exists():
                dcm2niix_bin = str(candidate)
                break
        # Also check user-local bin
        if not dcm2niix_bin:
            candidate = Path.home() / "Library" / "Python" / "3.9" / "bin" / "dcm2niix"
            if candidate.exists():
                dcm2niix_bin = str(candidate)
    if not dry_run and not dcm2niix_bin:
        log.error("dcm2niix not found. Install with: pip install dcm2niix")
        sys.exit(1)
    if dcm2niix_bin:
        log.info(f"Using dcm2niix: {dcm2niix_bin}")

    bids_dir.mkdir(parents=True, exist_ok=True)
    meta = {}   # Skip idaSearch CSV (slow on Google Drive); participants.tsv will have IDs only
    participants = {}

    # Try extracted directories first; only scan zips if no extracted data found
    extracted_entries = discover_extracted_dicoms(zip_dir, subject_filter)

    entries = {}
    if extracted_entries:
        for key, dcm_dir_path in extracted_entries.items():
            entries[key] = (None, dcm_dir_path)   # (None, dir_path) = already on disk
    else:
        log.info("No extracted directories found — scanning zip archives...")
        zip_entries = discover_zip_dicoms(zip_dir, subject_filter)
        for key, (zf_path, dcm_names) in zip_entries.items():
            entries[key] = (zf_path, dcm_names)
        log.info(f"Found {len(entries)} entries from zips")

    if not entries:
        log.warning("No DICOM data found!")
        return

    # Group entries by subject
    from collections import defaultdict
    import time

    subj_entries = defaultdict(list)
    for (subj, seq_name, img_id), (src_zip, dcm_data) in sorted(entries.items()):
        classification = classify_sequence(seq_name)
        if classification is None:
            continue
        subj_entries[subj].append((seq_name, img_id, src_zip, dcm_data, classification))

    run_counters = {}
    total_subj = len(subj_entries)

    for subj_idx, (subj, seq_list) in enumerate(sorted(subj_entries.items()), 1):
        sub_label = f"sub-{subj}"
        ses_label = "ses-BL"

        if dry_run:
            for seq_name, img_id, src_zip, dcm_data, (datatype, suffix, pe_dir) in seq_list:
                run_key = (subj, datatype, suffix, pe_dir)
                if run_key not in run_counters:
                    run_counters[run_key] = 1
                run_num = run_counters[run_key]
                run_counters[run_key] += 1
                pe_label = f"_dir-{pe_dir}" if pe_dir else ""
                run_label = f"_run-{run_num:02d}"
                task_label = "_task-rest" if suffix == "bold" else ""
                bids_fname = f"{sub_label}_{ses_label}{task_label}{pe_label}{run_label}_{suffix}"
                log.info(f"  [DRY] {subj}/{seq_name} → {datatype}/{bids_fname}")
            if subj not in participants:
                participants[subj] = {"age": "", "sex": "", "group": ""}
            continue

        log.info(f"\n  === Subject {subj_idx}/{total_subj}: {subj} ({len(seq_list)} sequences) ===")
        t_subj = time.time()

        for seq_name, img_id, src_zip, dcm_data, (datatype, suffix, pe_dir) in seq_list:
            try:
                run_key = (subj, datatype, suffix, pe_dir)
                if run_key not in run_counters:
                    run_counters[run_key] = 1
                run_num = run_counters[run_key]
                run_counters[run_key] += 1

                pe_label = f"_dir-{pe_dir}" if pe_dir else ""
                run_label = f"_run-{run_num:02d}"
                task_label = "_task-rest" if suffix == "bold" else ""
                bids_subdir = bids_dir / sub_label / ses_label / datatype
                bids_fname = f"{sub_label}_{ses_label}{task_label}{pe_label}{run_label}_{suffix}"
                bids_subdir.mkdir(parents=True, exist_ok=True)

                if src_zip is None:
                    dcm_dir = dcm_data
                else:
                    tmp_dcm_dir = tempfile.mkdtemp(prefix="dcm2bids_")
                    for dcm_name in dcm_data:
                        with zipfile.ZipFile(src_zip, "r") as zf:
                            data = zf.read(dcm_name)
                            (Path(tmp_dcm_dir) / Path(dcm_name).name).write_bytes(data)
                    dcm_dir = tmp_dcm_dir

                log.info(f"    Converting {seq_name} ...")
                sys.stdout.flush()
                cmd = [
                    dcm2niix_bin,
                    "-w", "1",
                    "-z", "y", "-b", "y", "-ba", "y",
                    "-f", bids_fname,
                    "-o", str(bids_subdir),
                    dcm_dir,
                ]
                try:
                    result = subprocess.run(cmd, capture_output=True, text=True, timeout=600)
                    returncode = result.returncode
                except subprocess.TimeoutExpired:
                    log.error(f"    TIMEOUT {seq_name} — skipping")
                    returncode = -1

                if src_zip is not None:
                    shutil.rmtree(tmp_dcm_dir, ignore_errors=True)

                # Remove dcm2niix-generated derived maps (ADC, phase, Eq_1)
                # that are not BIDS-valid entities
                for derived in bids_subdir.glob(f"{bids_fname}*"):
                    if any(tag in derived.name for tag in ["_ADC.", "_ph.", "_Eq_1.", "_TRACEW."]):
                        derived.unlink()

                # Inject TaskName into BOLD JSON sidecars (BIDS requirement)
                if suffix == "bold":
                    for json_file in bids_subdir.glob(f"{bids_fname}*.json"):
                        try:
                            with open(json_file, "r") as jf:
                                jdata = json.load(jf)
                            if "TaskName" not in jdata:
                                jdata["TaskName"] = "rest"
                                with open(json_file, "w") as jf:
                                    json.dump(jdata, jf, indent=2)
                        except (json.JSONDecodeError, OSError):
                            pass

                if returncode != 0:
                    log.warning(f"    WARN {seq_name}: dcm2niix returned {returncode}")
                else:
                    log.info(f"    OK {seq_name} → {bids_fname}")
            except Exception as e:
                log.error(f"    ERROR {seq_name}: {e} — skipping")
                continue

        if subj not in participants:
            participants[subj] = {"age": "", "sex": "", "group": ""}
        log.info(f"  === Done {subj}: {time.time() - t_subj:.0f}s ===")

    # Write participants.tsv
    tsv_path = bids_dir / "participants.tsv"
    with open(tsv_path, "w") as f:
        f.write("participant_id\tage\tsex\tresearch_group\n")
        for subj in sorted(participants.keys()):
            info = participants[subj]
            f.write(f"sub-{subj}\t{info.get('age','')}\t"
                    f"{info.get('sex','')}\t{info.get('group','')}\n")
    log.info(f"Wrote {tsv_path}")

    # Write dataset_description.json
    desc_path = bids_dir / "dataset_description.json"
    desc = {
        "Name": "PPMI Apathy MRI",
        "BIDSVersion": "1.9.0",
        "DatasetType": "raw",
        "License": "PPMI Data Use Agreement",
        "GeneratedBy": [{"Name": "dicom_to_bids.py"}],
    }
    with open(desc_path, "w") as f:
        json.dump(desc, f, indent=2)
    log.info(f"Wrote {desc_path}")

    # Summary
    print("\n" + "=" * 60)
    print("DICOM → BIDS CONVERSION SUMMARY")
    print("=" * 60)
    print(f"Subjects processed: {len(participants)}")
    print(f"Output directory:   {bids_dir}")
    print("=" * 60)


def generate_shell_script(zip_dir: Path, bids_dir: Path,
                          subject_filter: str = None) -> str:
    """Generate a standalone bash script that does the DICOM→BIDS conversion.
    Google Drive FUSE deadlocks when Python child processes do file I/O,
    but standalone shell commands work fine. So Python handles discovery
    and classification, and the shell script does all actual file I/O."""

    dcm2niix_bin = shutil.which("dcm2niix")
    if not dcm2niix_bin:
        import site
        for d in site.getusersitepackages(), *site.getsitepackages():
            candidate = Path(d).parent.parent / "bin" / "dcm2niix"
            if candidate.exists():
                dcm2niix_bin = str(candidate)
                break
        if not dcm2niix_bin:
            candidate = Path.home() / "Library" / "Python" / "3.9" / "bin" / "dcm2niix"
            if candidate.exists():
                dcm2niix_bin = str(candidate)
    if not dcm2niix_bin:
        log.error("dcm2niix not found")
        sys.exit(1)

    bids_dir.mkdir(parents=True, exist_ok=True)
    extracted_entries = discover_extracted_dicoms(zip_dir, subject_filter)

    entries = {}
    if extracted_entries:
        for key, dcm_dir_path in extracted_entries.items():
            entries[key] = dcm_dir_path
    else:
        log.error("No extracted directories found (zip mode not supported for script generation)")
        sys.exit(1)

    from collections import defaultdict
    subj_entries = defaultdict(list)
    for (subj, seq_name, img_id), dcm_dir in sorted(entries.items()):
        classification = classify_sequence(seq_name)
        if classification is None:
            continue
        subj_entries[subj].append((seq_name, img_id, dcm_dir, classification))

    run_counters = {}
    total_subj = len(subj_entries)
    participants = set()

    lines = [
        '#!/bin/bash',
        '# Auto-generated DICOM→BIDS conversion script',
        f'# Generated: {__import__("datetime").datetime.now().isoformat()}',
        f'# Group: {zip_dir.name}  Subjects: {total_subj}',
        f'DCM2NIIX="{dcm2niix_bin}"',
        f'BIDS_DIR="{bids_dir}"',
        '',
        'TOTAL=0; OK=0; FAIL=0',
        'convert_one() {',
        '    local DCM_SRC="$1" BIDS_SUB="$2" BIDS_FNAME="$3" SEQ_NAME="$4"',
        '    local TMP_IN=$(mktemp -d /tmp/dcm2bids_in.XXXXXX)',
        '    local TMP_OUT=$(mktemp -d /tmp/dcm2bids_out.XXXXXX)',
        '',
        '    # Copy DICOMs to local temp (Google Drive FUSE needs sequential cp)',
        '    find "$DCM_SRC" -maxdepth 1 \\( -name "*.dcm" -o -name "*.DCM" \\) -print0 | \\',
        '        while IFS= read -r -d \'\' f; do cp "$f" "$TMP_IN/"; done',
        '    find "$DCM_SRC" -maxdepth 1 -not -name "*.*" -not -type d -print0 | \\',
        '        while IFS= read -r -d \'\' f; do cp "$f" "$TMP_IN/"; done',
        '',
        '    local NFILES=$(ls "$TMP_IN" 2>/dev/null | wc -l | tr -d " ")',
        '    if [ "$NFILES" -eq 0 ]; then',
        '        echo "    SKIP $SEQ_NAME (no DICOMs found)"',
        '        rm -rf "$TMP_IN" "$TMP_OUT"',
        '        return 1',
        '    fi',
        '',
        '    mkdir -p "$BIDS_SUB"',
        '    "$DCM2NIIX" -w 1 -z y -b y -ba y -f "$BIDS_FNAME" -o "$TMP_OUT" "$TMP_IN" >/dev/null 2>&1',
        '    local RC=$?',
        '',
        '    # Move outputs to BIDS dir',
        '    for f in "$TMP_OUT"/*; do',
        '        [ -f "$f" ] && mv "$f" "$BIDS_SUB/"',
        '    done',
        '    rm -rf "$TMP_IN" "$TMP_OUT"',
        '',
        '    TOTAL=$((TOTAL + 1))',
        '    if [ $RC -eq 0 ]; then',
        '        OK=$((OK + 1))',
        '        echo "    OK $SEQ_NAME → $BIDS_FNAME ($NFILES DICOMs)"',
        '    else',
        '        FAIL=$((FAIL + 1))',
        '        echo "    FAIL $SEQ_NAME (exit $RC)"',
        '    fi',
        '    return $RC',
        '}',
        '',
    ]

    for subj_idx, (subj, seq_list) in enumerate(sorted(subj_entries.items()), 1):
        sub_label = f"sub-{subj}"
        ses_label = "ses-BL"
        participants.add(subj)

        lines.append(f'echo "=== Subject {subj_idx}/{total_subj}: {subj} ({len(seq_list)} sequences) ==="')

        for seq_name, img_id, dcm_dir, (datatype, suffix, pe_dir) in seq_list:
            run_key = (subj, datatype, suffix, pe_dir)
            if run_key not in run_counters:
                run_counters[run_key] = 1
            run_num = run_counters[run_key]
            run_counters[run_key] += 1

            pe_label = f"_dir-{pe_dir}" if pe_dir else ""
            run_label = f"_run-{run_num:02d}"
            task_label = "_task-rest" if suffix == "bold" else ""
            bids_subdir = f"$BIDS_DIR/{sub_label}/{ses_label}/{datatype}"
            bids_fname = f"{sub_label}_{ses_label}{task_label}{pe_label}{run_label}_{suffix}"

            lines.append(
                f'convert_one "{dcm_dir}" "{bids_subdir}" "{bids_fname}" "{seq_name}"'
            )
        lines.append('')

    # Write participants.tsv and dataset_description.json
    lines.append('# Write BIDS metadata')
    lines.append('cat > "$BIDS_DIR/participants.tsv" << \'PTSV\'')
    lines.append('participant_id\tage\tsex\tresearch_group')
    for subj in sorted(participants):
        lines.append(f'sub-{subj}\t\t\t')
    lines.append('PTSV')
    lines.append('')

    desc_json = json.dumps({
        "Name": "PPMI Apathy MRI",
        "BIDSVersion": "1.9.0",
        "DatasetType": "raw",
        "License": "PPMI Data Use Agreement",
        "GeneratedBy": [{"Name": "dicom_to_bids.py"}],
    }, indent=2)
    lines.append(f"cat > \"$BIDS_DIR/dataset_description.json\" << 'DESC'")
    lines.append(desc_json)
    lines.append('DESC')
    lines.append('')
    lines.append('echo ""')
    lines.append('echo "============================================================"')
    lines.append('echo "CONVERSION COMPLETE: $OK OK, $FAIL failed out of $TOTAL"')
    lines.append(f'echo "Subjects: {total_subj}  Output: $BIDS_DIR"')
    lines.append('echo "============================================================"')

    return "\n".join(lines) + "\n"


def main():
    parser = argparse.ArgumentParser(description="PPMI DICOM → BIDS conversion")
    parser.add_argument("--subject", type=str, help="Process single subject ID")
    parser.add_argument("--dry-run", action="store_true", help="Preview without converting")
    parser.add_argument("--group", type=str, help="Process single group (e.g. Apathy_BL)")
    parser.add_argument("--generate-script", type=str, metavar="PATH",
                        help="Generate standalone bash script (avoids Google Drive FUSE issues)")
    parser.add_argument("--source-root", type=str, metavar="PATH",
                        help="Override source MRI_Data root (e.g. Google Drive path)")
    args = parser.parse_args()

    # Allow overriding source root while keeping BIDS output on OneDrive
    if args.source_root:
        source_root = Path(args.source_root)
        zip_dirs = [source_root / d.name for d in DICOM_ZIP_DIRS]
    else:
        zip_dirs = list(DICOM_ZIP_DIRS)

    for zip_dir in zip_dirs:
        group_name = zip_dir.name  # e.g. "Apathy_BL"
        if args.group and group_name != args.group:
            continue
        if not zip_dir.exists():
            log.warning(f"Directory not found: {zip_dir}")
            continue

        group_bids = BIDS_DIR / group_name

        if args.generate_script:
            log.info(f"Generating conversion script for {group_name} ...")
            script = generate_shell_script(zip_dir, group_bids,
                                           subject_filter=args.subject)
            script_path = Path(args.generate_script)
            script_path.write_text(script)
            script_path.chmod(0o755)
            log.info(f"Script written to: {script_path}")
            log.info(f"Run it with:  bash {script_path}")
        else:
            log.info(f"\n{'='*60}\nProcessing group: {group_name}\nBIDS output: {group_bids}\n{'='*60}")
            convert_to_bids(zip_dir, group_bids,
                            subject_filter=args.subject,
                            dry_run=args.dry_run)


if __name__ == "__main__":
    main()

#!/bin/bash
# Role: Per-participant PBS chain on Imperial CX3: MRIQC 24.0.2, structural stream
#     (structural_preprocessing.sh), diffusion, fMRIPrep 25.2.5, XCP-D 26.0.1 (linc, 36P,
#     4S256Parcels), provenance (06_collect_provenance.py). Needs the containers, the TemplateFlow
#     cache and the FreeSurfer licence set up by setup_hpc.sh.
# Manuscript item: Methods (MRI acquisition and preprocessing). The diffusion job calls 03_diffusion_preproc.sh (no
#     diffusion analysis is reported) and the provenance job calls 06_collect_provenance.py; neither
#     is in this repository, so run with --steps mriqc,struct,func,xcpd. Participants whose default
#     input was a short calibration series were re-run with the filters written by
#     select_runs_by_length.py (this folder).
# -----------------------------------------------------------------------
# submit_pipeline.sh — Submit the full PPMI MRI pipeline on Imperial CX3
#
# Architecture (PBS job dependencies, per-subject):
#   Job 1: MRIQC        (Singularity)               — QC before processing
#   Job 2: Structural   (FreeSurfer + ANTs modules)  — depends on MRIQC
#   Job 3: Diffusion    (MRtrix3 + FSL GPU)          — depends on Structural
#   Job 4: fMRIPrep     (Singularity)                — depends on Structural
#   Job 5: XCP-D        (Singularity)                — depends on fMRIPrep
#   Job 6: Provenance   (Python)                     — depends on all
#
# Each subject gets its own set of jobs (no PBS arrays — CX3 restricts them).
# Dependencies are per-subject: sub-X's diffusion waits for sub-X's structural.
#
# Usage:
#   bash submit_pipeline.sh [--steps all|mriqc|struct|diff|func|xcpd|prov]
#                                [--subjects sub-XXXXXX,sub-YYYYYY]
#                                [--dry-run]
#
# The --steps flag lets you run individual stages (e.g. after fixing failures).
# Default: all steps with dependencies.
# -----------------------------------------------------------------------
set -euo pipefail

# Session timepoint (default BL). For 24M follow-up, invoke as:
#   SESSION=ses-24M bash submit_pipeline.sh --groups ALL_24M
SESSION="${SESSION:-ses-BL}"

# ── Parse arguments ──────────────────────────────────────────────────
STEPS="all"
SUBJECT_FILTER=""
GROUP_FILTER=""
DRY_RUN=false

while [[ $# -gt 0 ]]; do
    case $1 in
        --steps)     STEPS="$2"; shift 2;;
        --subjects)  SUBJECT_FILTER="$2"; shift 2;;
        --groups)    GROUP_FILTER="$2"; shift 2;;
        --dry-run)   DRY_RUN=true; shift;;
        *)           echo "Unknown option: $1"; exit 1;;
    esac
done

# ── Paths & Storage Strategy ─────────────────────────────────────────
# HOME (~325 GB free): BIDS input, final outputs, FreeSurfer, QC, provenance
# EPHEMERAL (~11 TB): ALL intermediate processing, containers, fMRIPrep work
#
# Per-subject storage estimate (with cleanup):
#   Structural finals:  ~200 MB (brain, mask, warps, FAST, FIRST, sienax)
#   FreeSurfer:         ~1 GB
#   Diffusion finals:   ~500 MB (preprocessed DWI 4D, tensor, all scalars, MNI maps)
#   fMRIPrep output:    ~3 GB   (confounds, preprocessed BOLD, boldref)
#   XCP-D output:       ~1 GB   (connectivity matrices, denoised BOLD, parcellations)
#   MRIQC:              ~50 MB  (reports, IQMs)
#   Total per subject:  ~5.75 GB on Home
#   For 55 subjects:    ~316 GB → fits in current Home (325 GB free)
#   For 100 subjects:   ~575 GB → may need cleanup or quota increase
#
# Strategy: process on EPHEMERAL, copy finals to HOME, delete ephemeral working data

RDS_HOME="/rds/general/user/${USER}/home"
RDS_EPHEMERAL="/rds/general/user/${USER}/ephemeral"

BIDS_ROOT="${BIDS_ROOT:-${RDS_HOME}/PPMI/MRI_Data/BIDS}"
SCRIPTS_DIR="${SCRIPTS_DIR:-$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)}"
LOGS_DIR="${LOGS_DIR:-${RDS_HOME}/PPMI/MRI_Data/logs}"

# Final outputs go to Home (persistent)
DERIV_DIR="${DERIV_DIR:-${RDS_HOME}/PPMI/MRI_Data/derivatives}"

# All processing happens on Ephemeral (fast, large, temporary)
EPHEMERAL_DERIV="${RDS_EPHEMERAL}/PPMI/derivatives"
WORK_DIR="${RDS_EPHEMERAL}/PPMI/work"

# Containers on Ephemeral (15 GB — can re-pull if needed)
CONTAINER_DIR="${RDS_EPHEMERAL}/PPMI/containers"
TEMPLATEFLOW_DIR="${RDS_HOME}/.cache/templateflow"
FS_LICENSE="${RDS_HOME}/freesurfer_license.txt"

mkdir -p "$LOGS_DIR" "$WORK_DIR" "$EPHEMERAL_DERIV" "$DERIV_DIR"

# Also ensure containers exist on ephemeral (copy from home if setup was done there)
if [ ! -d "$CONTAINER_DIR" ] && [ -d "${RDS_HOME}/containers" ]; then
    echo "Copying containers from Home to Ephemeral ..."
    mkdir -p "$CONTAINER_DIR"
    cp -n "${RDS_HOME}/containers/"*.sif "$CONTAINER_DIR/" 2>/dev/null || true
fi

echo "============================================"
echo "PPMI MRI Pipeline — Job Submission"
echo "============================================"
echo "BIDS_ROOT:        $BIDS_ROOT"
echo "DERIV_DIR (home):  $DERIV_DIR"
echo "EPHEMERAL_DERIV:   $EPHEMERAL_DERIV"
echo "WORK_DIR:          $WORK_DIR"
echo "CONTAINER_DIR:     $CONTAINER_DIR"
echo "Steps:             $STEPS"
echo "Dry run:           $DRY_RUN"
echo ""
echo "Storage strategy:"
echo "  Processing on: Ephemeral (11 TB)"
echo "  Finals kept on: Home (~5.75 GB/subject)"
echo ""

# ── Build subject list ───────────────────────────────────────────────
# BIDS subdirs by session:
#   ses-BL:   Apathy_BL, Both_BL, ICB_BL, Neither_BL
#   ses-24M:  ALL_24M (single folder holds all 24M-timepoint subjects)
SUBJ_LIST="$LOGS_DIR/subject_list_${SESSION}.txt"
> "$SUBJ_LIST"
if [ "$SESSION" = "ses-BL" ]; then
    CANDIDATE_SUBDIRS=(Apathy_BL Both_BL ICB_BL Neither_BL)
else
    CANDIDATE_SUBDIRS=(ALL_24M)
fi
for subdir in "${CANDIDATE_SUBDIRS[@]}"; do
    if [ -d "${BIDS_ROOT}/${subdir}" ]; then
        for d in "${BIDS_ROOT}/${subdir}"/sub-*/; do
            [ -d "$d" ] || continue
            subj=$(basename "$d" | sed 's/sub-//')
            # Only include subjects whose session directory actually exists
            if [ -d "${d}${SESSION}" ]; then
                echo "${subj}|${subdir}" >> "$SUBJ_LIST"
            fi
        done
    fi
done
sort -t'|' -k1,1 -u "$SUBJ_LIST" -o "$SUBJ_LIST"

# Filter by group if requested (comma-separated, e.g. --groups Neither_BL)
if [ -n "$GROUP_FILTER" ]; then
    FILTERED="$LOGS_DIR/subject_list_filtered.txt"
    > "$FILTERED"
    IFS=',' read -ra GRPS <<< "$GROUP_FILTER"
    for g in "${GRPS[@]}"; do
        grep "|${g}$" "$SUBJ_LIST" >> "$FILTERED" || true
    done
    SUBJ_LIST="$FILTERED"
fi

# Filter by subject if requested (comma-separated, e.g. --subjects sub-XXXXXX,sub-YYYYYY)
# Use a different filename from the group filter so the two can compose
# when both --groups and --subjects are passed (prior bug: both wrote to
# subject_list_filtered.txt, so --subjects truncated the group-filtered list).
if [ -n "$SUBJECT_FILTER" ]; then
    SUBJ_FILTERED="$LOGS_DIR/subject_list_subj_filtered.txt"
    > "$SUBJ_FILTERED"
    IFS=',' read -ra SUBS <<< "$SUBJECT_FILTER"
    for s in "${SUBS[@]}"; do
        s_id=$(echo "$s" | sed 's/sub-//')
        grep "^${s_id}|" "$SUBJ_LIST" >> "$SUBJ_FILTERED" || true
    done
    SUBJ_LIST="$SUBJ_FILTERED"
fi


N_SUBJ=$(wc -l < "$SUBJ_LIST" | tr -d ' ')
echo "Subjects: $N_SUBJ"
[ "$N_SUBJ" -eq 0 ] && { echo "ERROR: No subjects found"; exit 1; }
echo ""

# ── Validate prerequisites ───────────────────────────────────────────
ERRORS=0
if [[ "$STEPS" == "all" || "$STEPS" == *"func"* || "$STEPS" == *"mriqc"* ]]; then
    if [ ! -f "$CONTAINER_DIR/fmriprep-25.2.5.sif" ]; then
        echo "ERROR: fMRIPrep container not found. Run setup_hpc.sh first."
        ERRORS=$((ERRORS+1))
    fi
    if [ ! -f "$FS_LICENSE" ]; then
        echo "ERROR: FreeSurfer license not found at $FS_LICENSE"
        ERRORS=$((ERRORS+1))
    fi
fi
[ "$ERRORS" -gt 0 ] && exit 1

# ── Helper: submit or print ─────────────────────────────────────────
submit_job() {
    local script="$1"
    local dep="${2:-}"
    if $DRY_RUN; then
        echo -n "  [DRY RUN] Would submit: $(basename "$script")"
        [ -n "$dep" ] && echo -n "  (after $dep)"
        echo ""
        echo "DRY_${RANDOM}"
    else
        local JOB_ID=""
        if [ -n "$dep" ]; then
            JOB_ID=$(qsub -W depend=afterok:"$dep" "$script" 2>&1) || {
                echo "  WARNING: dependency submit failed ($JOB_ID), submitting without dependency" >&2
                JOB_ID=$(qsub "$script" 2>&1) || { echo "  ERROR: qsub failed: $JOB_ID" >&2; return 1; }
            }
        else
            JOB_ID=$(qsub "$script" 2>&1) || { echo "  ERROR: qsub failed: $JOB_ID" >&2; return 1; }
        fi
        # Clean up job ID (remove trailing whitespace/newlines)
        JOB_ID=$(echo "$JOB_ID" | tr -d '[:space:]' | tail -1)
        echo "  Submitted: $(basename "$script") → $JOB_ID" >&2
        echo "$JOB_ID"
    fi
}

# ── Track all job IDs per stage for provenance dependency ────────────
ALL_MRIQC_JOBS=""
ALL_STRUCT_JOBS=""
ALL_DIFF_JOBS=""
ALL_FMRIPREP_JOBS=""
ALL_XCPD_JOBS=""

# ======================================================================
# PER-SUBJECT JOB LOOP
# Each subject gets individual jobs with per-subject dependencies.
# ======================================================================
echo "============================================"
echo "Submitting per-subject jobs..."
echo "============================================"
echo ""

IDX=0
while IFS='|' read -r SUBJECT BIDS_SUBDIR; do
    IDX=$((IDX + 1))
    echo "[$IDX/$N_SUBJ] sub-${SUBJECT} ($BIDS_SUBDIR)"

    SUBJ_MRIQC_JOB=""
    SUBJ_STRUCT_JOB=""
    SUBJ_DIFF_JOB=""
    SUBJ_FMRIPREP_JOB=""
    SUBJ_XCPD_JOB=""

    # ==================================================================
    # MRIQC
    # ==================================================================
    if [[ "$STEPS" == "all" || "$STEPS" == *"mriqc"* ]]; then
        MRIQC_SCRIPT="$LOGS_DIR/job_mriqc_${SUBJECT}.pbs"
        cat > "$MRIQC_SCRIPT" << PBSEOF
#!/bin/bash
#PBS -N mriqc_${SUBJECT}
#PBS -l walltime=12:00:00
#PBS -l select=1:ncpus=4:mem=32gb
#PBS -q v1_medium72
#PBS -o ${LOGS_DIR}/mriqc_${SUBJECT}_${SESSION}.out
#PBS -e ${LOGS_DIR}/mriqc_${SUBJECT}_${SESSION}.err

echo "[\$(date)] MRIQC for sub-${SUBJECT} (${BIDS_SUBDIR})"

BIDS_IN="${BIDS_ROOT}/${BIDS_SUBDIR}"
OUT_DIR="${DERIV_DIR}/mriqc"
WORK="${WORK_DIR}/mriqc_${SUBJECT}_${SESSION}"
mkdir -p "\$OUT_DIR" "\$WORK"

singularity run --cleanenv \
    -B "\$BIDS_IN":/data:ro \
    -B "\$OUT_DIR":/out \
    -B "\$WORK":/work \
    "${CONTAINER_DIR}/mriqc-24.0.2.sif" \
    /data /out participant \
    --participant-label "${SUBJECT}" \
    --no-sub \
    -m T1w bold \
    --nprocs 4 \
    --mem 14 \
    -w /work

rm -rf "\$WORK"
echo "[\$(date)] MRIQC COMPLETE for sub-${SUBJECT}"
PBSEOF
        SUBJ_MRIQC_JOB=$(submit_job "$MRIQC_SCRIPT")
        ALL_MRIQC_JOBS="${ALL_MRIQC_JOBS:+${ALL_MRIQC_JOBS}:}${SUBJ_MRIQC_JOB}"
    fi

    # ==================================================================
    # Structural
    # ==================================================================
    if [[ "$STEPS" == "all" || "$STEPS" == *"struct"* ]]; then
        STRUCT_SCRIPT="$LOGS_DIR/job_struct_${SUBJECT}.pbs"
        cat > "$STRUCT_SCRIPT" << PBSEOF
#!/bin/bash
#PBS -N struct_${SUBJECT}
#PBS -l walltime=24:00:00
#PBS -l select=1:ncpus=4:mem=16gb
#PBS -q v1_medium72
#PBS -o ${LOGS_DIR}/struct_${SUBJECT}_${SESSION}.out
#PBS -e ${LOGS_DIR}/struct_${SUBJECT}_${SESSION}.err

module load tools/prod
module load ANTs/2.3.5-foss-2021a
module load FreeSurfer/7.4.1-centos8_x86_64
module load FSL/6.0.7.17
unset PYTHONPATH 2>/dev/null || true
export LD_LIBRARY_PATH=\$(echo "\$LD_LIBRARY_PATH" | tr ':' '\n' | grep -v MCR | tr '\n' ':' | sed 's/:\$//')
source \$FSLDIR/etc/fslconf/fsl.sh 2>/dev/null || true

export BIDS_ROOT="${BIDS_ROOT}"
export DERIV_DIR="${EPHEMERAL_DERIV}"
export CLEANUP_INTERMEDIATES=true
export FS_LICENSE_FILE="${FS_LICENSE}"
export SESSION="${SESSION}"

echo "[\$(date)] Structural for sub-${SUBJECT} (\${SESSION}, processing on ephemeral)"

bash "${SCRIPTS_DIR}/structural_preprocessing.sh" "${SUBJECT}" "${BIDS_SUBDIR}"

# ── Copy finals from EPHEMERAL → HOME ──────────────────────────────
SRC="${EPHEMERAL_DERIV}/structural/sub-${SUBJECT}/${SESSION}"
DST="${DERIV_DIR}/structural/sub-${SUBJECT}/${SESSION}"
mkdir -p "\$DST/fast" "\$DST/first" "\$DST/sienax" "\$DST/qc"

for f in T1w_brain.nii.gz T1w_brain_mask.nii.gz \
         T1w2MNI_warp.nii.gz MNI2T1w_warp.nii.gz T1w2MNI.mat \
         FLAIR_to_T1w.nii.gz FLAIR2T1w.mat provenance.json; do
    [ -f "\$SRC/\$f" ] && cp "\$SRC/\$f" "\$DST/"
done
cp "\$SRC"/fast/T1w_pve_*.nii.gz "\$DST/fast/" 2>/dev/null || true
cp "\$SRC"/fast/T1w_seg.nii.gz "\$DST/fast/" 2>/dev/null || true
cp "\$SRC"/first/*.nii.gz "\$DST/first/" 2>/dev/null || true
cp "\$SRC"/sienax/report.sienax "\$DST/sienax/" 2>/dev/null || true
cp "\$SRC"/qc/*.png "\$DST/qc/" 2>/dev/null || true

# FreeSurfer — copy only essential files to HOME (stats + surfaces for QC)
# Full recon stays on ephemeral (fMRIPrep reads it from there)
FS_SRC="${EPHEMERAL_DERIV}/freesurfer/sub-${SUBJECT}_${SESSION}"
FS_DST="${DERIV_DIR}/freesurfer/sub-${SUBJECT}_${SESSION}"
if [ -d "\$FS_SRC" ]; then
    mkdir -p "\$FS_DST/stats" "\$FS_DST/mri" "\$FS_DST/label"
    # Stats (TIV from aseg, parcellation volumes) — essential for covariates
    cp "\$FS_SRC"/stats/aseg.stats "\$FS_DST/stats/" 2>/dev/null || true
    cp "\$FS_SRC"/stats/lh.aparc.stats "\$FS_DST/stats/" 2>/dev/null || true
    cp "\$FS_SRC"/stats/rh.aparc.stats "\$FS_DST/stats/" 2>/dev/null || true
    # aparc+aseg for registration / parcellation lookup
    cp "\$FS_SRC"/mri/aparc+aseg.mgz "\$FS_DST/mri/" 2>/dev/null || true
    cp "\$FS_SRC"/mri/brain.mgz "\$FS_DST/mri/" 2>/dev/null || true
    # Keep full recon on ephemeral — fMRIPrep needs it
fi

# NOTE: do NOT rm structural or FreeSurfer from ephemeral —
# diffusion needs structural for T1 registration, fMRIPrep needs FreeSurfer

echo "[\$(date)] Structural COMPLETE — finals on Home (~\$(du -sh "\$DST" 2>/dev/null | cut -f1))"
PBSEOF
        SUBJ_STRUCT_JOB=$(submit_job "$STRUCT_SCRIPT" "$SUBJ_MRIQC_JOB")
        ALL_STRUCT_JOBS="${ALL_STRUCT_JOBS:+${ALL_STRUCT_JOBS}:}${SUBJ_STRUCT_JOB}"
    fi

    # ==================================================================
    # Diffusion (GPU)
    # ==================================================================
    if [[ "$STEPS" == "all" || "$STEPS" == *"diff"* ]]; then
        # Check if DWI data exists before creating job
        DWI_DIR="${BIDS_ROOT}/${BIDS_SUBDIR}/sub-${SUBJECT}/${SESSION}/dwi"
        if [ -d "$DWI_DIR" ] && [ -n "$(find "$DWI_DIR" -name "*_dwi.nii.gz" 2>/dev/null | head -1)" ]; then
            DIFF_SCRIPT="$LOGS_DIR/job_diff_${SUBJECT}.pbs"
            cat > "$DIFF_SCRIPT" << PBSEOF
#!/bin/bash
#PBS -N diff_${SUBJECT}
#PBS -l walltime=12:00:00
#PBS -l select=1:ncpus=8:mem=24gb
#PBS -q v1_medium72
#PBS -o ${LOGS_DIR}/diff_${SUBJECT}_${SESSION}.out
#PBS -e ${LOGS_DIR}/diff_${SUBJECT}_${SESSION}.err

module load tools/prod
module load ANTs/2.3.5-foss-2021a
module load MRtrix/3.0.4-foss-2022b
module load FSL/6.0.7.17
unset PYTHONPATH 2>/dev/null || true
export LD_LIBRARY_PATH=\$(echo "\$LD_LIBRARY_PATH" | tr ':' '\n' | grep -v MCR | tr '\n' ':' | sed 's/:\$//')
source \$FSLDIR/etc/fslconf/fsl.sh 2>/dev/null || true

# Fail loudly on any downstream error (including missing preproc script).
set -euo pipefail

export BIDS_ROOT="${BIDS_ROOT}"
export DERIV_DIR="${EPHEMERAL_DERIV}"
export CLEANUP_INTERMEDIATES=true
export SESSION="${SESSION}"

# Structural output should already be on ephemeral from the preceding job.
# If not (e.g. re-running diffusion alone), symlink from HOME.
STRUCT_EPH="${EPHEMERAL_DERIV}/structural/sub-${SUBJECT}"
STRUCT_HOME="${DERIV_DIR}/structural/sub-${SUBJECT}"
if [ ! -d "\$STRUCT_EPH" ] && [ -d "\$STRUCT_HOME" ]; then
    mkdir -p "${EPHEMERAL_DERIV}/structural"
    ln -sf "\$STRUCT_HOME" "\$STRUCT_EPH"
    echo "  Symlinked structural from HOME → ephemeral"
fi

echo "[\$(date)] Diffusion for sub-${SUBJECT} (GPU, ephemeral)"
# Hard-fail if preproc script is missing rather than silently printing "COMPLETE".
PREPROC_SCRIPT="${SCRIPTS_DIR}/03_diffusion_preproc.sh"
if [ ! -f "\$PREPROC_SCRIPT" ]; then
    echo "ERROR: diffusion preprocessing script not found at \$PREPROC_SCRIPT" >&2
    exit 1
fi
bash "\$PREPROC_SCRIPT" "${SUBJECT}" "${BIDS_SUBDIR}"

# ── Copy lightweight finals from EPHEMERAL → HOME ────────────────
# Skip dwi_preprocessed.nii.gz (4D, ~400MB) — stays on ephemeral
# Only copy MNI scalar maps, dtifit, registrations, provenance
SRC="${EPHEMERAL_DERIV}/diffusion/sub-${SUBJECT}/${SESSION}"
DST="${DERIV_DIR}/diffusion/sub-${SUBJECT}/${SESSION}"
mkdir -p "\$DST/dtifit" "\$DST/qc"

for f in FA_MNI.nii.gz MD_MNI.nii.gz RD_MNI.nii.gz AD_MNI.nii.gz FA2T1.nii.gz \
         b0_brain.nii.gz b0_brain_mask.nii.gz \
         b0toT1_0GenericAffine.mat FA2T1.mat FA2MNI_warp.nii.gz \
         provenance.json; do
    [ -f "\$SRC/\$f" ] && cp "\$SRC/\$f" "\$DST/"
done
# Preprocessed bval/bvec (tiny) but NOT the 4D nifti
cp "\$SRC"/dwi_preprocessed.bval "\$DST/" 2>/dev/null || true
cp "\$SRC"/dwi_preprocessed.bvec "\$DST/" 2>/dev/null || true
# dtifit scalar maps only (skip tensor)
cp "\$SRC"/dtifit/dti_FA.nii.gz "\$DST/dtifit/" 2>/dev/null || true
cp "\$SRC"/dtifit/dti_MD.nii.gz "\$DST/dtifit/" 2>/dev/null || true
cp "\$SRC"/dtifit/dti_RD.nii.gz "\$DST/dtifit/" 2>/dev/null || true
cp "\$SRC"/qc/* "\$DST/qc/" 2>/dev/null || true

# Keep full diffusion on ephemeral (tractography may need it later)
# To reclaim space later: rm -rf ${EPHEMERAL_DERIV}/diffusion/
[ -L "\$STRUCT_EPH" ] && rm -f "\$STRUCT_EPH"

echo "[\$(date)] Diffusion COMPLETE — finals on Home (~\$(du -sh "\$DST" 2>/dev/null | cut -f1))"
PBSEOF
            SUBJ_DIFF_JOB=$(submit_job "$DIFF_SCRIPT" "$SUBJ_STRUCT_JOB")
            ALL_DIFF_JOBS="${ALL_DIFF_JOBS:+${ALL_DIFF_JOBS}:}${SUBJ_DIFF_JOB}"
        else
            echo "  (no DWI data — skipping diffusion)"
        fi
    fi

    # ==================================================================
    # fMRIPrep
    # ==================================================================
    if [[ "$STEPS" == "all" || "$STEPS" == *"func"* ]]; then
        FUNC_DIR="${BIDS_ROOT}/${BIDS_SUBDIR}/sub-${SUBJECT}/${SESSION}/func"
        if [ -d "$FUNC_DIR" ] && [ -n "$(find "$FUNC_DIR" -name "*_bold.nii.gz" 2>/dev/null | head -1)" ]; then
            FMRIPREP_SCRIPT="$LOGS_DIR/job_fmriprep_${SUBJECT}.pbs"
            cat > "$FMRIPREP_SCRIPT" << PBSEOF
#!/bin/bash
#PBS -N fprep_${SUBJECT}
#PBS -l walltime=48:00:00
#PBS -l select=1:ncpus=8:mem=32gb
#PBS -q v1_medium72
#PBS -o ${LOGS_DIR}/fmriprep_${SUBJECT}_${SESSION}.out
#PBS -e ${LOGS_DIR}/fmriprep_${SUBJECT}_${SESSION}.err

echo "[\$(date)] fMRIPrep for sub-${SUBJECT}"

BIDS_IN="${BIDS_ROOT}/${BIDS_SUBDIR}"
# Each subject gets its own output dir to avoid fsaverage race condition
# (multiple jobs writing to the same /out/sourcedata/freesurfer/fsaverage)
EPH_OUT="${EPHEMERAL_DERIV}/fmriprep_${SUBJECT}"
WORK="${WORK_DIR}/fmriprep_${SUBJECT}_${SESSION}"
mkdir -p "\$EPH_OUT" "\$WORK"

# Reuse FreeSurfer from structural pipeline (saves 6-12h per subject)
# Full recon lives on EPHEMERAL (only stats copied to HOME)
# Create per-subject staging dir so fMRIPrep can write fsaverage freely
FS_EPH="${EPHEMERAL_DERIV}/freesurfer"
FS_HOME="${DERIV_DIR}/freesurfer"
FS_STAGE="\${WORK}/fs_stage"
# Always create the bind target so singularity -B doesn't fail on missing source.
# Without this, a missing FreeSurfer recon caused a silent <20s container crash
# (log: "mount source ... fs_stage doesn't exist"), and the pipeline falsely
# reported subject as "processed" because stdout only logs PBS timing.
mkdir -p "\$FS_STAGE"
FS_REUSE=""
# Check ephemeral first (full recon), then HOME (for old subjects)
if [ -d "\${FS_EPH}/sub-${SUBJECT}_${SESSION}/mri" ]; then
    ln -sfn "\${FS_EPH}/sub-${SUBJECT}_${SESSION}" "\${FS_STAGE}/sub-${SUBJECT}"
    FS_REUSE="--fs-subjects-dir /fsdir"
    echo "  Reusing FreeSurfer from ephemeral (staged)"
elif [ -d "\${FS_HOME}/sub-${SUBJECT}_${SESSION}/mri" ]; then
    ln -sfn "\${FS_HOME}/sub-${SUBJECT}_${SESSION}" "\${FS_STAGE}/sub-${SUBJECT}"
    FS_REUSE="--fs-subjects-dir /fsdir"
    echo "  Reusing FreeSurfer from HOME (staged)"
else
    echo "  WARNING: no FreeSurfer recon found for sub-${SUBJECT} — fmriprep will run recon-all in-line (+6-12 h)"
fi

export APPTAINERENV_TEMPLATEFLOW_HOME="/templateflow"

singularity run --cleanenv \
    -B "\$BIDS_IN":/data:ro \
    -B "\$EPH_OUT":/out \
    -B "\$WORK":/work \
    -B "${FS_LICENSE}":/opt/freesurfer/license.txt:ro \
    -B "${TEMPLATEFLOW_DIR}":/templateflow:ro \
    -B "\${FS_STAGE}":/fsdir \
    "${CONTAINER_DIR}/fmriprep-25.2.5.sif" \
    /data /out participant \
    --participant-label "${SUBJECT}" \
    --output-spaces MNI152NLin2009cAsym:res-2 MNI152NLin6Asym:res-2 \
    --fd-spike-threshold 0.5 \
    --skull-strip-t1w force \
    --nprocs 8 \
    --mem-mb 30000 \
    --skip-bids-validation \
    --stop-on-first-crash \
    --notrack \
    \$FS_REUSE \
    -w /work
# NOTE: --skip-bids-validation needed because NM files and empty subjects
# in the BIDS group dirs cause the built-in validator to reject the dataset.
# --use-aroma removed in fMRIPrep 23.1. Denoising handled by XCP-D.

rm -rf "\$WORK"

# ── Keep fMRIPrep output on EPHEMERAL for XCP-D to read ──────────
# Only copy lightweight files to Home now; XCP-D job will do final cleanup
FPREP_DST="${DERIV_DIR}/fmriprep"
mkdir -p "\$FPREP_DST"
cp "\$EPH_OUT/sub-${SUBJECT}.html" "\$FPREP_DST/" 2>/dev/null || true
cp "\$EPH_OUT/dataset_description.json" "\$FPREP_DST/" 2>/dev/null || true

echo "[\$(date)] fMRIPrep COMPLETE — output on Ephemeral (~\$(du -sh "\$EPH_OUT" 2>/dev/null | cut -f1)), waiting for XCP-D"
PBSEOF
            SUBJ_FMRIPREP_JOB=$(submit_job "$FMRIPREP_SCRIPT" "$SUBJ_STRUCT_JOB")
            ALL_FMRIPREP_JOBS="${ALL_FMRIPREP_JOBS:+${ALL_FMRIPREP_JOBS}:}${SUBJ_FMRIPREP_JOB}"
        else
            echo "  (no BOLD data — skipping fMRIPrep)"
        fi
    fi

    # ==================================================================
    # XCP-D
    # ==================================================================
    if [[ "$STEPS" == "all" || "$STEPS" == *"xcpd"* ]]; then
        XCPD_SCRIPT="$LOGS_DIR/job_xcpd_${SUBJECT}.pbs"
        cat > "$XCPD_SCRIPT" << PBSEOF
#!/bin/bash
#PBS -N xcpd_${SUBJECT}
#PBS -l walltime=24:00:00
#PBS -l select=1:ncpus=4:mem=48gb
#PBS -q v1_medium72
#PBS -o ${LOGS_DIR}/xcpd_${SUBJECT}_${SESSION}.out
#PBS -e ${LOGS_DIR}/xcpd_${SUBJECT}_${SESSION}.err

echo "[\$(date)] XCP-D for sub-${SUBJECT}"

# Read fMRIPrep from EPHEMERAL (not Home)
FMRIPREP_EPH="${EPHEMERAL_DERIV}/fmriprep_${SUBJECT}"
EPH_OUT="${EPHEMERAL_DERIV}/xcpd_${SUBJECT}"
WORK="${WORK_DIR}/xcpd_${SUBJECT}_${SESSION}"
mkdir -p "\$EPH_OUT" "\$WORK"

if [ ! -d "\${FMRIPREP_EPH}/sub-${SUBJECT}" ]; then
    echo "No fMRIPrep output for sub-${SUBJECT} on ephemeral — skipping"
    exit 0
fi

singularity run --cleanenv \
    -B "\$FMRIPREP_EPH":/data:ro \
    -B "\$EPH_OUT":/out \
    -B "\$WORK":/work \
    -B "${FS_LICENSE}":/opt/freesurfer/license.txt:ro \
    "${CONTAINER_DIR}/xcp_d-26.0.1.sif" \
    /data /out participant \
    --mode linc \
    --file-format nifti \
    --participant-label "${SUBJECT}" \
    --despike y \
    --lower-bpf 0.01 --upper-bpf 0.08 \
    --fd-thresh 0.5 \
    -p 36P \
    --smoothing 6 \
    --nprocs 4 \
    --mem-mb 14000 \
    -w /work

rm -rf "\$WORK"

# ── Copy only lightweight XCP-D outputs to HOME ─────────────────
# Connectivity matrices + timeseries (needed for group analysis) → HOME
# Denoised BOLD 4D stays on ephemeral (bulk of the ~1GB)
XCPD_SRC="\$EPH_OUT/sub-${SUBJECT}"
XCPD_DST="${DERIV_DIR}/xcp_d/sub-${SUBJECT}"
if [ -d "\$XCPD_SRC" ]; then
    # Copy timeseries, connectivity matrices, QC, and atlased outputs only
    mkdir -p "\$XCPD_DST/${SESSION}/func"
    # Timeseries + correlation matrices (what group analysis reads) — tiny
    cp "\$XCPD_SRC"/${SESSION}/func/*_timeseries.tsv "\$XCPD_DST/${SESSION}/func/" 2>/dev/null || true
    cp "\$XCPD_SRC"/${SESSION}/func/*_relmat.tsv "\$XCPD_DST/${SESSION}/func/" 2>/dev/null || true
    cp "\$XCPD_SRC"/${SESSION}/func/*_stat-*.tsv "\$XCPD_DST/${SESSION}/func/" 2>/dev/null || true
    # Confounds / motion summary
    cp "\$XCPD_SRC"/${SESSION}/func/*_motion.tsv "\$XCPD_DST/${SESSION}/func/" 2>/dev/null || true
    cp "\$XCPD_SRC"/${SESSION}/func/*_outliers.tsv "\$XCPD_DST/${SESSION}/func/" 2>/dev/null || true
    cp "\$XCPD_SRC"/${SESSION}/func/*_design.tsv "\$XCPD_DST/${SESSION}/func/" 2>/dev/null || true
    # Parcellated images (small)
    cp "\$XCPD_SRC"/${SESSION}/func/*_boldref.nii.gz "\$XCPD_DST/${SESSION}/func/" 2>/dev/null || true
    # Skip *_bold.nii.gz (denoised 4D — ~500MB-1GB)
fi
cp "\$EPH_OUT/sub-${SUBJECT}.html" "${DERIV_DIR}/xcp_d/" 2>/dev/null || true
cp "\$EPH_OUT/dataset_description.json" "${DERIV_DIR}/xcp_d/" 2>/dev/null || true

# Full XCP-D + fMRIPrep stay on ephemeral
# To clean up later: rm -rf ${EPHEMERAL_DERIV}/fmriprep_* ${EPHEMERAL_DERIV}/xcpd_*

echo "[\$(date)] XCP-D COMPLETE — lightweight outputs on Home, bulk on Ephemeral"
PBSEOF
        SUBJ_XCPD_JOB=$(submit_job "$XCPD_SCRIPT" "$SUBJ_FMRIPREP_JOB")
        ALL_XCPD_JOBS="${ALL_XCPD_JOBS:+${ALL_XCPD_JOBS}:}${SUBJ_XCPD_JOB}"
    fi

    echo ""
done < "$SUBJ_LIST"

# ======================================================================
# Provenance collection (single job, after all subjects complete)
# ======================================================================
if [[ "$STEPS" == "all" || "$STEPS" == *"prov"* ]]; then
    echo "--- Generating Provenance collection job ---"
    PROV_SCRIPT="$LOGS_DIR/job_provenance.pbs"
    cat > "$PROV_SCRIPT" << PBSEOF
#!/bin/bash
#PBS -N ppmi_prov
#PBS -l walltime=00:30:00
#PBS -l select=1:ncpus=1:mem=4gb
#PBS -q v1_small24
#PBS -o ${LOGS_DIR}/provenance.out
#PBS -e ${LOGS_DIR}/provenance.err

module load tools/prod
module load Python/3.12.3-GCCcore-13.3.0

export DERIV_DIR="${DERIV_DIR}"
cd "${SCRIPTS_DIR}"
python3 06_collect_provenance.py

echo "[\$(date)] Provenance collection COMPLETE"

echo ""
echo "=== STORAGE REPORT ==="
echo "Home derivatives:"
du -sh "${DERIV_DIR}"/* 2>/dev/null || true
echo ""
echo "Ephemeral (should be mostly empty):"
du -sh "${EPHEMERAL_DERIV}"/* 2>/dev/null || echo "  (clean)"
echo ""
quota -s 2>/dev/null | head -10 || true
PBSEOF
    # Provenance depends on the last jobs from diff + xcpd
    PROV_DEP=""
    # Collect all final-stage job IDs
    if [ -n "$ALL_DIFF_JOBS" ]; then
        PROV_DEP="$ALL_DIFF_JOBS"
    fi
    if [ -n "$ALL_XCPD_JOBS" ]; then
        PROV_DEP="${PROV_DEP:+${PROV_DEP}:}${ALL_XCPD_JOBS}"
    fi
    PROV_JOB=$(submit_job "$PROV_SCRIPT" "$PROV_DEP")
    echo ""
fi

# ── Summary ──────────────────────────────────────────────────────────
echo "============================================"
echo "JOB SUBMISSION SUMMARY"
echo "============================================"
echo "Subjects:     $N_SUBJ"
echo ""
echo "Jobs per subject: MRIQC → Structural → Diffusion (GPU)"
echo "                                     → fMRIPrep → XCP-D"
echo "Final job:        Provenance (after all complete)"
echo ""
echo "Storage plan:"
echo "  Home (finals):     ~5.75 GB/subject × $N_SUBJ = ~$((N_SUBJ * 6)) GB"
echo "  Ephemeral (temp):  cleaned after each subject"
echo ""
N_JOBS=$((N_SUBJ * 5 + 1))
echo "Total jobs submitted: ~${N_JOBS} (5 per subject + 1 provenance)"
echo ""
echo "Monitor:  module purge && qstat -u $USER"
echo "Logs:     $LOGS_DIR/"
echo "Storage:  quota -s"

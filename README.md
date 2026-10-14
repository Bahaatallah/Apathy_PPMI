# Apathy, impulse-control behaviours and cognitive decline in Parkinson's disease: analysis code

Code for the main analyses and figures of *Apathy marks brain network and genetic vulnerability to cognitive decline in Parkinson's disease* (Parkinson's Precision Medicine Initiative cohort; primary connectivity sample N = 305): cohort construction and longitudinal mixed-effects models, the bifactor measurement model of global cognition, imaging preprocessing, the Network-Based Statistic (NBS), the connectivity-to-outcome analyses and the scripts that draw Figures 1 to 5. Scripts for the supplementary analyses are available from the corresponding author on request. No participant-level data are included.

Repository: `<GITHUB_URL_PLACEHOLDER>`. Archived release: `<ZENODO_DOI_PLACEHOLDER>`. Licence: MIT (`LICENSE`). Citation metadata: `CITATION.cff`.

## Layout

| File | Role |
|---|---|
| `01_clinical_cognition/Publicaiton_analysis.m` | Master clinical script: cohort, apathy and ICB phenotyping, four-group classification, longitudinal tables with the bifactor g scores; its workspace is the input of every MATLAB script below |
| `02_cfa_bifactor/ppm_cfa_Gemeni_bothModels.R` | Bifactor confirmatory factor analysis; writes `cfa_factor_scores.csv` (`bifactor_g`) |
| `03_imaging_preproc/config.py`, `01_dicom_to_bids.py` | DICOM to BIDS |
| `03_imaging_preproc/05a_setup_hpc.sh`, `05b_submit_pipeline.sh`, `02_structural_preproc.sh` | Containers; MRIQC, structural, fMRIPrep and XCP-D jobs |
| `03_imaging_preproc/make_bids_filters.py` | Run selection by run length (fMRIPrep `--bids-filter-file`) |
| `04_group_nbs/10_fc_group.py`, `group_covariates.py` | NBS (edge-wise GLM, residualised-regressor permutation, component FWE); helper module |
| `05_mediation_moderation/extract_nbs_suprathreshold_fc.py` | Per-participant component features |
| `05_mediation_moderation/Figure4.m` | Cognitive slopes, component-to-slope correlations, mediation |
| `05_mediation_moderation/Figure4_Structural.m` | 42 regional volumes, 24-month hippocampal atrophy |
| `05_mediation_moderation/Figure4_render.m` | Figure 4 |
| `06_figures/Figure1.m`, `Figure2.m`, `Figure3.m` | Figure 1 and Table 1; Figure 2; Figure 5 |
| `06_figures/plot_nbs_combined_revised.py`, `plot_nbs_circos_revised.py`, `atlas_4S256Parcels_labels_xcpd.tsv` | Figure 3; parcel-order label table |

## Software

MATLAB R2024b with the Statistics and Machine Learning Toolbox; R 4.1.1 (lavaan, semTools, semPlot, dplyr, readr, tidyr, readxl, ggplot2, patchwork, corrplot); Python 3 (numpy, scipy, pandas; for Figure 3 also matplotlib, nibabel, nilearn with the Schaefer-2018 and fsaverage5 data, Pillow); dcm2niix; MRIQC 24.0.2, fMRIPrep 25.2.5 and XCP-D 26.0.1 containers; FSL 6.0.7.17, FreeSurfer 7.4.1, ANTs 2.3.5; a PBS Pro cluster.

## Configuration

`PPMI_CLINICAL_ROOT` holds `ppmi_data/` (the curated data cut `PPMI_Curated_Data_Cut_Public_20260223.xlsx` and the instrument CSVs under `download/`, outputs under `MATLAB_Results/`), `MRI_Data/all_structural_measures.csv` and the output tree `revision/`; `PPMI_IMAGING_ROOT` holds `MRI_Data/` (BIDS and derivatives). Unset, both default to the working directory. The cluster scripts default to `$HOME/PPMI/` and `/rds/general/user/$USER/` and accept `BIDS_ROOT`, `DERIV_DIR`, `GROUP_DIR`, `SESSION` and `FSLDIR`.

## Run order

1. `Rscript 02_cfa_bifactor/ppm_cfa_Gemeni_bothModels.R`, then copy `ppmi_data/MATLAB_Results/cfa_factor_scores.csv` to `revision/figures/regenerated/` (the output folder of step 2).
2. MATLAB, with `01_clinical_cognition`, `05_mediation_moderation` and `06_figures` on the path: `Publicaiton_analysis`, then `delete(findall(groot,'Type','figure')); save('master_workspace.mat','-v7.3')`.
3. MATLAB, one fresh session per script on the loaded workspace: `FIG_SUPPRESS_ID_EXPORTS = true; Figure1`, then `Figure2`, then `Figure3`.
4. Cluster: `python3 01_dicom_to_bids.py`; `bash 05a_setup_hpc.sh`; `bash 05b_submit_pipeline.sh --steps mriqc,struct,func,xcpd`.
5. Extended sample only: `python3 make_bids_filters.py --roster <roster>`, then fMRIPrep with `--bids-filter-file filters/<sub>_filter.json` and XCP-D as in step 4.
6. `DERIV_DIR=<derivatives> GROUP_DIR=<run> python3 04_group_nbs/10_fc_group.py --analysis connectome_continuous --nperm 10000 --t-thresh 3.1 --subjects-file <design roster> --covariates <covariate table> --fmriprep-root /nonexistent_neutralized`.
7. `python3 05_mediation_moderation/extract_nbs_suprathreshold_fc.py --nbs-dir <run>/fc_connectome_continuous --covariates <covariate table>`; place the output with `nbs_edge_strengths.csv` in one feature folder.
8. `python3 06_figures/plot_nbs_combined_revised.py --nbs-dir <run>/fc_connectome_continuous --atlas-labels 06_figures/atlas_4S256Parcels_labels_xcpd.tsv --fig-dir <out> --contrasts main_apathy diff_apathy_gt_quip interaction_apathyxquip --out-name Figure3`.
9. MATLAB on the loaded workspace: `outputDir = '<fig4>'; EXT_FC_DIR = '<features>'; EXT_FEATURE_SET = 'P305'; EXT_MEDIATOR = 'diff'; Figure4; Figure4_Structural`, then `Figure4_render('<features>', '<fig4>', '<slopes csv>', '<out stem>')`.

## Inputs prepared outside these scripts

- **Covariate table**: one row per participant with `subject_id, group, age, sex, site, TIV, mean_FD, education_years, disease_duration_months, updrs_iii_off, ledd, moca, apathy_composite, depression_factor, quip, hoehn_yahr, scanner_model, handedness`. Demographics, site and baseline scores from the curated cut; TIV is FreeSurfer eTIV; `mean_FD` is the mean framewise displacement of the XCP-D motion file of the analysed run.
- **Design roster**: lines `sub-<id>|<group>` (group = BIDS folder) for the participants with a usable baseline task-free run on the conventional echo-planar protocol (200 to 260 volumes, TR of at least 2.0 s), complete covariates and adequate parcel coverage (parcels missing in more than 10% of participants blanked for everyone, 241 of 256 nodes kept; participants with more than 50 missing parcels excluded).
- **All-pairs features** `nbs_edge_strengths.csv`: per participant, the mean Fisher-z over all node pairs of each FWE-significant component (`<contrast>_<sign>_comp<k>_mean_fc`) plus the covariate columns.
- **Participant slopes** (`subject_id, slope_g`): the change in g per year of each participant, the fixed time effect plus the participant's random slope from the reduced LME of `Figure4.m` section 2.
- **Structural tables** (baseline and 24 months): one row per participant with FreeSurfer aseg volumes (`fs_vol_*`), aparc volumes and thickness (`fs_lh_*`, `fs_rh_*`) and eTIV (`fs_EstimatedTotalIntraCranialVol`) from the outputs of `02_structural_preproc.sh`.
- **Atlas label table** (included): `index, label, hemisphere, atlas_block` in the parcel order of the XCP-D 4S256Parcels outputs; any 4S256Parcels relmat or timeseries TSV can be passed to `--atlas-labels` instead.

Outputs of these scripts carry participant identifiers (rosters, feature and slope tables, `Figure1.m` exports without `FIG_SUPPRESS_ID_EXPORTS`) and must not be redistributed.

## Data access

PPMI data (www.ppmi-info.org/access-data-specimens/download-data; RRID:SCR_006431) cannot be redistributed; access requires application to PPMI. The scripts expect the curated cut and instrument CSVs under `$PPMI_CLINICAL_ROOT/ppmi_data/` and the imaging data (DICOM converted to BIDS) under `$PPMI_IMAGING_ROOT/MRI_Data/`.

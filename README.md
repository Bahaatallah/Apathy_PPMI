# Apathy, impulse-control behaviours and cognitive decline in Parkinson's disease: analysis code

Code for the analyses and figures of *Apathy marks brain network and genetic vulnerability to cognitive decline in Parkinson's disease* (Parkinson's Precision Medicine Initiative cohort; primary connectivity sample N = 305). The repository holds only the code behind the reported results; scripts for the supplementary analyses are available from the corresponding author on request. No participant-level data are included.

Repository: https://github.com/Bahaatallah/Apathy_PPMI. Licence: MIT (`LICENSE`). Citation metadata: `CITATION.cff`.

## Layout

```
behavioural/                         clinical, cognitive and genetic analyses (MATLAB, R)
  cognition_bifactor_cfa.R           bifactor measurement model of global cognition (g), reliability and longitudinal invariance
  cohort_phenotyping.m               cohort, apathy and ICB phenotyping, four groups, long-format modelling tables
  figure1_cohort_characterisation.m  Figure 1, Table 1, Supplementary Figure S2 and Table S1
  figure2_cognitive_trajectories.m   Figure 2, linear mixed-effects models of cognitive decline and their sensitivity refits
  figure5_genetic_moderation.m       Figure 5, genetic moderation models
imaging/
  preprocessing/                     DICOM to BIDS, MRIQC, FreeSurfer, fMRIPrep and XCP-D (cluster)
    config.py                        paths
    dicom_to_bids.py                 DICOM to BIDS conversion
    setup_hpc.sh                     containers and folders (run once)
    submit_pipeline.sh               per-participant job chain
    structural_preprocessing.sh      T1w preparation and FreeSurfer recon-all (structural job)
    select_runs_by_length.py         fMRIPrep bids-filter files selecting the longest task-free run
  connectivity/                      Network-Based Statistic and Figure 3 (Python)
    nbs_group_analysis.py            NBS: edge-wise GLM, residualised-regressor permutation, component FWE
    extract_component_features.py    per-participant component features (mean FC over supra-threshold edges)
    figure3_connectivity_components.py  Figure 3 (surface maps and chord diagrams)
    circos_helpers.py                chord-diagram helpers imported by the Figure 3 script
    atlas_4S256Parcels_labels.tsv    parcel-order label table of the 4S256Parcels atlas
  connectivity_outcomes/             connectivity to cognition and structure, Figure 4 (MATLAB)
    fc_cognitive_slope_mediation.m   cognitive slopes, component-to-slope correlations, mediation
    fc_structural_correlates.m       42 regional volumes, 24-month hippocampal atrophy
    figure4_connectivity_outcomes.m  Figure 4
```

## Software

MATLAB R2024b with the Statistics and Machine Learning Toolbox; R 4.1.1 (lavaan, dplyr, readr, readxl); Python 3 (numpy, scipy, pandas; for Figure 3 also matplotlib, nibabel, nilearn with the Schaefer-2018 and fsaverage5 data, Pillow); dcm2niix; MRIQC 24.0.2, fMRIPrep 25.2.5 and XCP-D 26.0.1 containers; FSL 6.0.7.17, FreeSurfer 7.4.1, ANTs 2.3.5; a PBS Pro cluster.

## Configuration

`PPMI_CLINICAL_ROOT` holds `ppmi_data/` (the curated data cut `PPMI_Curated_Data_Cut_Public_20260223.xlsx` and the instrument CSVs under `download/`, outputs under `MATLAB_Results/`), `MRI_Data/all_structural_measures.csv` and the output tree `revision/`; `PPMI_IMAGING_ROOT` holds `MRI_Data/` (BIDS and derivatives). Unset, both default to the working directory. The cluster scripts default to `$HOME/PPMI/` and `/rds/general/user/$USER/` and accept `BIDS_ROOT`, `DERIV_DIR`, `GROUP_DIR`, `SESSION` and `FSLDIR`.

## Run order

Behavioural analyses:

1. `Rscript behavioural/cognition_bifactor_cfa.R`, then copy `ppmi_data/MATLAB_Results/cfa_factor_scores.csv` to `revision/figures/regenerated/` (the output folder of step 2).
2. MATLAB, with `behavioural/` and `imaging/connectivity_outcomes/` on the path: `cohort_phenotyping`, then `delete(findall(groot,'Type','figure')); save('master_workspace.mat','-v7.3')`.
3. MATLAB, one fresh session per script on the loaded workspace: `figure1_cohort_characterisation`, then `figure2_cognitive_trajectories`, then `figure5_genetic_moderation`.

Imaging preprocessing (cluster):

4. `python3 imaging/preprocessing/dicom_to_bids.py`; `bash imaging/preprocessing/setup_hpc.sh`; `bash imaging/preprocessing/submit_pipeline.sh --steps mriqc,struct,func,xcpd`.
5. Extended sample only: `python3 imaging/preprocessing/select_runs_by_length.py --roster <roster>`, then fMRIPrep with `--bids-filter-file filters/<sub>_filter.json` and XCP-D as in step 4.

Connectivity:

6. `DERIV_DIR=<derivatives> GROUP_DIR=<run> python3 imaging/connectivity/nbs_group_analysis.py --analysis connectome_continuous --nperm 10000 --t-thresh 3.1 --subjects-file <design roster> --covariates <covariate table> --fmriprep-root /nonexistent_neutralized`. The nested covariate models of the Supplement add `--include-extended --extended-covs <list>`, the binary-QUIP model `--quip-binary`.
7. `python3 imaging/connectivity/extract_component_features.py --nbs-dir <run>/fc_connectome_continuous --covariates <covariate table>`; place the output with `nbs_edge_strengths.csv` in one feature folder.
8. `python3 imaging/connectivity/figure3_connectivity_components.py --nbs-dir <run>/fc_connectome_continuous --atlas-labels imaging/connectivity/atlas_4S256Parcels_labels.tsv --fig-dir <out> --contrasts main_apathy diff_apathy_gt_quip interaction_apathyxquip --out-name Figure3`.

Connectivity to outcomes:

9. MATLAB on the loaded workspace: `outputDir = '<fig4>'; EXT_FC_DIR = '<features>'; EXT_FEATURE_SET = 'P305'; EXT_MEDIATOR = 'diff'; fc_cognitive_slope_mediation; fc_structural_correlates`, then `figure4_connectivity_outcomes('<features>', '<fig4>', '<fig4>/participant_slopes_g.csv', '<out stem>')`. `EXT_FEATURE_SET = 'X343'` repeats the chain on the extended sample.

Output file names inside the scripts (for example `Fig3_Genetics_Motivation_Cognition` for Figure 5) are those of the manuscript's figure files.

## Inputs prepared outside these scripts

- **Covariate table**: one row per participant with `subject_id, group, age, sex, site, TIV, mean_FD, education_years, disease_duration_months, updrs_iii_off, ledd, moca, apathy_composite, depression_factor, quip, hoehn_yahr, scanner_model, handedness`. Demographics, site and baseline scores from the curated cut; TIV is FreeSurfer eTIV; `mean_FD` is the mean framewise displacement of the XCP-D motion file of the analysed run.
- **Design roster**: lines `sub-<id>|<group>` (group = BIDS folder) for the participants with a usable baseline task-free run on the conventional echo-planar protocol (200 to 260 volumes, TR of at least 2.0 s), complete covariates and adequate parcel coverage (parcels missing in more than 10% of participants blanked for everyone, 241 of 256 nodes kept; participants with more than 50 missing parcels excluded).
- **All-pairs features** `nbs_edge_strengths.csv`: per participant, the mean Fisher-z over all node pairs of each FWE-significant component (`<contrast>_<sign>_comp<k>_mean_fc`) plus the covariate columns.
- **Structural tables** (baseline and 24 months): one row per participant with FreeSurfer aseg volumes (`fs_vol_*`), aparc volumes and thickness (`fs_lh_*`, `fs_rh_*`) and eTIV (`fs_EstimatedTotalIntraCranialVol`) from the outputs of `structural_preprocessing.sh`.
- **Atlas label table** (included): `index, label, hemisphere, atlas_block` in the parcel order of the XCP-D 4S256Parcels outputs; any 4S256Parcels relmat or timeseries TSV can be passed to `--atlas-labels` instead.

Outputs of these scripts carry participant identifiers (rosters, feature and slope tables) and must not be redistributed.

## Data access

PPMI data (www.ppmi-info.org/access-data-specimens/download-data; RRID:SCR_006431) cannot be redistributed; access requires application to PPMI. The scripts expect the curated cut and instrument CSVs under `$PPMI_CLINICAL_ROOT/ppmi_data/` and the imaging data (DICOM converted to BIDS) under `$PPMI_IMAGING_ROOT/MRI_Data/`.

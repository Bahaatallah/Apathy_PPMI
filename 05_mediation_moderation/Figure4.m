% Role: Connectivity-to-outcome chain, step 1 of 3 for Figure 4: participant-specific cognitive
%     slopes (reduced LME, section 2), path a and path b partial correlations, the mediation
%     bootstrap and its sensitivity models; writes Table_Fig4_Sensitivity.csv,
%     Table_Fig4_Mediation.csv, Table_Fig4_Mediation_Sensitivity.csv and the other Table_Fig4_*
%     tables to outputDir. Step 2 is Figure4_Structural.m (same workspace), step 3 Figure4_render.m.
%     The figure this script draws itself (Fig4_FC_Apathy_CogDecline) is not the manuscript figure,
%     which Figure4_render.m draws.
% Run: on the workspace left by 01_clinical_cognition/Publicaiton_analysis.m, after setting
%     outputDir (output folder), EXT_FC_DIR (folder holding nbs_suprathreshold_fc.csv and
%     nbs_edge_strengths.csv), EXT_FEATURE_SET = 'P305' and EXT_MEDIATOR = 'diff' (the primary
%     configuration).
% Manuscript item: Figure 4a,b values; Results (connectivity and cognitive decline, mediation);
%     Supplementary Table S10.
%% EXT_FEATURE_SET selects the feature configuration (P305, the primary configuration) and EXT_MEDIATOR
%  the component carried forward (diff, the apathy-specific differential component). Okabe-Ito palette
%  shared with the other figure scripts.

fprintf('\n%s\n', repmat('=',1,72));
fprintf('  FIGURE 4: FC MEDIATES APATHY-RELATED COGNITIVE DECLINE\n');
fprintf('%s\n\n', repmat('=',1,72));

%% ========================================================================
%  0.  SETTINGS
% =========================================================================
if ~exist('outputDir','var') || isempty(outputDir), outputDir = pwd; end

% --- FD threshold: 0.3 or 0.4 (selects NBS component set & data path) ---
fdThreshold = 0.3;   % <<< CHANGE HERE: 0.3 or 0.4

% --- FC source for path b (FC -> cognition) ---
%     true  = suprathreshold FC (mean of NBS-flagged edges only — cleaner signal)
%     false = mean FC (mean of ALL edges in connected component — broader)
%     Path a (apathy -> FC) always uses suprathreshold.
%     When true, only apathy components (which have supra columns) are in
%     the primary FDR family; diff/shared fall back to mean FC for the
%     exploratory heatmap (and remain outside the primary family).
useSupraForPathB = true;

assert(ismember(fdThreshold, [0.3, 0.4]), 'fdThreshold must be 0.3 or 0.4');
fprintf('  FD threshold: %.1f\n', fdThreshold);
fprintf('  Path b FC source: %s\n\n', ternary(useSupraForPathB,'suprathreshold','mean'));

% feature CSVs from EXT_FC_DIR (set in the workspace); see header.
assert(exist('EXT_FC_DIR','var')==1 && isfolder(EXT_FC_DIR), 'EXT_FC_DIR must be set to the feature CSV folder');
if ~exist('EXT_FEATURE_SET','var') || isempty(EXT_FEATURE_SET), EXT_FEATURE_SET = 'C1'; end
if ~exist('EXT_MEDIATOR','var')    || isempty(EXT_MEDIATOR),    EXT_MEDIATOR    = 'apneg'; end
assert(ismember(EXT_FEATURE_SET, {'C1','C2','P307','X345','P305','R305','P305c','X343'}) && ismember(EXT_MEDIATOR, {'apneg','diff'}));
fprintf('  Feature set %s, mediator %s, dir %s\n', EXT_FEATURE_SET, EXT_MEDIATOR, EXT_FC_DIR);
basePath = EXT_FC_DIR;
fcPath_mean  = fullfile(basePath, 'nbs_edge_strengths.csv');
fcPath_supra = fullfile(basePath, 'nbs_suprathreshold_fc.csv');

% REVISED (R1.10): deterministic colour-blind-safe palette. Overrides any
% workspace COL so all revised figure scripts share one palette regardless
% of the order they are run in (previously a guarded fallback with a
% different, red/teal palette).
COL.bg=[1 1 1]; COL.txt=[0.15 0.15 0.18]; COL.accent=[0.20 0.40 0.65];
COL.grid=[0.85 0.85 0.85];
COL.apathy=[0.90 0.60 0.00];    % Okabe-Ito orange       (was [0.80 0.30 0.30])
COL.icb=[0.00 0.45 0.70];       % Okabe-Ito blue         (was [0.20 0.55 0.70])
COL.overlap=[0.45 0.20 0.50];   % deep purple            (was [0.55 0.25 0.60])
COL.neither=[0.00 0.60 0.50];   % Okabe-Ito bluish green (cross-script consistency)

scCol = [0.72 0.38 0.20; 0.20 0.40 0.65; 0.17 0.50 0.55; 0.55 0.25 0.60];
nBoot = 5000; rng(42);

%% ========================================================================
%  1.  LOAD DATA
% =========================================================================
fprintf('[1] Loading data...\n');

% Both FC files
assert(exist(fcPath_mean,'file')==2 && exist(fcPath_supra,'file')==2);
fcMean  = readtable(fcPath_mean);
fcSupra = readtable(fcPath_supra);
nFC = height(fcMean);
fprintf('    FC: %d subjects (mean + suprathreshold loaded)\n', nFC);

% Longitudinal table
if exist('lmeT','var'), Tlong=lmeT;
elseif exist('Tpd_lme','var'), Tlong=Tpd_lme;
else, error('No longitudinal table.'); end
fprintf('    Long: %d rows, %d subjects\n', height(Tlong), numel(unique(Tlong.PATNO)));

% Resolve variable names — strict, no fallback. Global cognition uses
% bifactor_g; subdomains use the raw obs composites (obs_att_exec/memory/
% visuo) and bifactor_g is added as a covariate when modelling subdomain
% effects, rather than being absorbed via pre-residualised obs_*_resid.
cogVarG = 'bifactor_g';
assert(ismember(cogVarG, Tlong.Properties.VariableNames), ...
    'Tlong is missing column %s — bifactor factor scores must be merged before running Figure 4.', cogVarG);
timeVar=''; for c={'time','YEAR','year_num'}, if ismember(c{1},Tlong.Properties.VariableNames), timeVar=c{1}; break; end, end
assert(~isempty(timeVar), 'Tlong is missing a time column (time/YEAR/year_num).');

% Domain variables (strict)
cogDomains  = {'global','exec','memory','visuo'};
cogVarNames = {'bifactor_g','obs_att_exec','obs_memory','obs_visuo'};
cogLabels   = {'Global Cognition','Executive','Memory','Visuospatial'};
for di = 1:numel(cogVarNames)
    assert(ismember(cogVarNames{di}, Tlong.Properties.VariableNames), ...
        'Tlong is missing column %s.', cogVarNames{di});
end
nDom = numel(cogVarNames);

% Sex -> numeric
sexFC=nan(nFC,1);
if iscell(fcMean.sex)||isstring(fcMean.sex)||iscategorical(fcMean.sex)
    sStr=string(fcMean.sex);
    sexFC(strcmpi(sStr,'M')|strcmpi(sStr,'Male'))=1;
    sexFC(strcmpi(sStr,'F')|strcmpi(sStr,'Female'))=0;
else, sexFC=double(fcMean.sex); end
%% --- FC component definitions (FD-threshold dependent) ---
% Updated 2026-04-16: Primary analysis N=310 Neither-expanded Freedman-Lane
% (10k perms, t=3.1). Post-Simpson-paradox fix.
% Four significant components total: apathy pos, apathy neg, diff_apathy_gt_quip
% neg, shared_pos pos. QUIP main and interaction are null.
%
% Suprathreshold CSV now generated with the corrected projection-based reduced
% model (matches _nbs_continuous in 10_fc_group.py) and includes all 4 sig
% contrasts (was previously apathy-only due to hardcoded extractor).
% Ergo: path a (predictor -> FC) and DAT -> FC analyses now extend to all 4
% components via supraIdx (rather than apathyIdx-only).
if fdThreshold == 0.3
    % MEAN FC (for path b: FC -> cognition) — from nbs_edge_strengths.csv
    meanFCvars = {'main_apathy_pos_comp1_mean_fc', ...           % Broad hyperconn: SomMot(13)+Default(11)+Ctrl(7)+Sal(6)+Thal(2)+Vis+DA+NAcc+Cereb  (51 edges, p=0.0014)
                  'main_apathy_neg_comp1_mean_fc', ...           % RH Somatomotor + DorsAttn decoupling                                         (14 edges, p=0.048)
                  'diff_apathy_gt_quip_neg_comp1_mean_fc', ...   % Apathy-dominant: SomMot(13)+DA(4) apathy decouples > QUIP                   (17 edges, p=0.037)
                  'shared_pos_pos_comp1_mean_fc'};               % Shared hyperconn: Visual(7)+Ctrl/Default pCun(4)+DA                          (15 edges, p=0.046)
    % SUPRATHRESHOLD FC (for path a: predictor -> FC) — all 4 contrasts now present
    supraFCvars = {'main_apathy_pos_comp1_suprathresh_fc', ...           % Broad hyperconnectivity
                   'main_apathy_neg_comp1_suprathresh_fc', ...           % RH Somatomotor-DA decoupling
                   'diff_apathy_gt_quip_neg_comp1_suprathresh_fc', ...   % Apathy-dominant SomMot-DA
                   'shared_pos_pos_comp1_suprathresh_fc'};               % Visual-precuneus hyperconn
    fcLabels = {'Apathy(+) Comp1: Broad hyperconnectivity (SomMot-DMN-Ctrl-Sal-Thal-Cereb)', ...
                'Apathy(-) Comp1: RH Somatomotor-DorsAttn decoupling', ...
                'Diff(-) Comp1: Apathy-dominant Somatomotor-DorsAttn', ...
                'Shared(+) Comp1: Visual-precuneus hyperconnectivity'};
    fcShort  = {'Ap(+)1: Broad hyper', 'Ap(-)1: RH SM-DA', ...
                'Diff(-)1: Ap-dom SM-DA', 'Sh(+)1: Vis-pCun'};
    apathyIdx  = 1:2;
    quipIdx    = [];       % NULL in N=310 (Simpson's paradox resolved — no independent QUIP FC signature)
    interIdx   = [];       % NULL in N=310
    diffIdx    = 3;
    sharedIdx  = 4;
    fcContrast = [1 1 4 5];   % 1=apathy, 2=quip(n/a), 3=interaction(n/a), 4=diff, 5=shared
    % Panel mapping (all indices are fi = index into meanFCvars)
    idxA = 1;  % Broad apathy hyperconnectivity (Panel A scatter — primary / mediation)
    idxB = 2;  % RH SomMot-DA decoupling        (Panel B scatter -> Supp)
    idxMediation = idxB;  % FC component used as mediator in Section 6
                          % (hypo network: apathy → FC → decline)
    idxStructural = diffIdx;  % Diff(-)1: Ap-dom SM-DA — L hipp structural correlate
                              % (Cook's D for Ap(-)1 dominated by a single participant →
                              %  Ap(-)1 Pearson is outlier-driven; Diff(-)1 is the
                              %  FDR-sig Spearman primary L hipp finding under M1).
    labelA = 'Broad apathy hyperconn';
    labelA_axis = 'Broad apathy hyperconn FC (Fisher z)';
    labelB = 'RH SomMot-DA decoupling';
    labelB_axis = 'RH Somatomotor-DorsAttn FC (Fisher z)';
    % Panel C: domain-specific scatter
    idxC_fc  = 1;  % Broad apathy component in meanFCvars for cognition link
    idxC_dom = 1;  % Executive domain index in cogLabels (update as needed)
    labelC = 'Broad apathy hyperconn';
    labelC_axis = 'Broad apathy hyperconn FC (Fisher z)';
    labelC_dom  = 'Executive';
    % ---- Feature-set adjustments ----------------------------------------
    if strcmp(EXT_FEATURE_SET, 'C2')
        fcLabels{2} = 'Apathy(-) Comp1: SomMot-DorsAttn decoupling (extended-defined, 11e, p=.090, not FWE-sig)';
        fcShort{2}  = 'Ap(-)1 ext (n.s.)';
        meanFCvars(end+1:end+2)  = {'interaction_apathyxquip_neg_comp1_mean_fc', ...
                                    'interaction_apathyxquip_pos_comp1_mean_fc'};
        supraFCvars(end+1:end+2) = {'interaction_apathyxquip_neg_comp1_suprathresh_fc', ...
                                    'interaction_apathyxquip_pos_comp1_suprathresh_fc'};
        fcLabels(end+1:end+2) = {'Int(-) Comp1: DorsAttn-Default-Control (extended-defined, 36e, p=.011)', ...
                                 'Int(+) Comp1: Default-mode (extended-defined, 18e, p=.048)'};
        fcShort(end+1:end+2)  = {'Int(-)1: DA-DMN', 'Int(+)1: DMN'};
        interIdx   = 5:6;
        fcContrast = [fcContrast 3 3];
    end
    % ---- N = 307 PRIMARY (protocol-stratified, corrected FD; 2026-09-07) ----
    %   components of group_fixed_neither_fdcorr307: hyper 48e p=.0025,
    %   hypo 18e p=.029, differential 18e p=.030; shared NOT significant
    %   (absent); interaction neg 20e p=.037 carried as a flagged extra only.
    if strcmp(EXT_FEATURE_SET, 'P307')
        meanFCvars  = {'main_apathy_pos_comp1_mean_fc', 'main_apathy_neg_comp1_mean_fc', ...
                       'diff_apathy_gt_quip_neg_comp1_mean_fc', 'interaction_apathyxquip_neg_comp1_mean_fc'};
        supraFCvars = {'main_apathy_pos_comp1_suprathresh_fc', 'main_apathy_neg_comp1_suprathresh_fc', ...
                       'diff_apathy_gt_quip_neg_comp1_suprathresh_fc', 'interaction_apathyxquip_neg_comp1_suprathresh_fc'};
        fcLabels = {'Apathy(+) Comp1: Broad hyperconnectivity (N=307, 48e, p=.0025)', ...
                    'Apathy(-) Comp1: Somatomotor-Salience decoupling (N=307, 18e, p=.029)', ...
                    'Diff(-) Comp1: Apathy-dominant Somatomotor-DorsAttn (N=307, 18e, p=.030)', ...
                    'Int(-) Comp1: Visual-DorsAttn-Control (N=307, 20e, p=.037, near threshold, not a finding)'};
        fcShort  = {'Ap(+)1: Broad hyper', 'Ap(-)1: SM-Sal', 'Diff(-)1: Ap-dom SM-DA', 'Int(-)1: Vis-DA (n.t.)'};
        apathyIdx = 1:2; quipIdx = []; interIdx = 4; diffIdx = 3; sharedIdx = [];
        fcContrast = [1 1 4 3];
        if ~ismember('interaction_apathyxquip_neg_comp1_suprathresh_fc', fcSupra.Properties.VariableNames)
            % e.g. the covariate-only fdcorr310 run: no interaction component
            meanFCvars(4) = []; supraFCvars(4) = []; fcLabels(4) = []; fcShort(4) = [];
            interIdx = []; fcContrast = [1 1 4];
        end
        idxA = 1; idxB = 2; idxMediation = idxB; idxStructural = diffIdx;
        labelB = 'SomMot-Salience decoupling'; labelB_axis = 'Apathy SomMot-Salience FC (Fisher z)';
    end
    % ---- N = 305 FINAL PRIMARY (protocol-stratified, corrected FD, coverage
    %      rule: 15 dropout parcels blanked, 2 participants dropped;
    %      2026-09-08). Components of group_fixed_neither_fdcorr305_cov:
    %      hyper 32e p=.0057, hypo comp1 16e p=.032 (somatomotor-salience,
    %      16 of the 18 edges of the 307 hypo component), hypo comp2 16e
    %      p=.032 (tied second cluster, dorsal-attention / somatomotor /
    %      salience), differential 18e p=.028, interaction neg 20e p=.026;
    %      shared not significant (absent). ----
    if strcmp(EXT_FEATURE_SET, 'P305')
        meanFCvars  = {'main_apathy_pos_comp1_mean_fc', 'main_apathy_neg_comp1_mean_fc', ...
                       'diff_apathy_gt_quip_neg_comp1_mean_fc', 'interaction_apathyxquip_neg_comp1_mean_fc', ...
                       'main_apathy_neg_comp2_mean_fc'};
        supraFCvars = {'main_apathy_pos_comp1_suprathresh_fc', 'main_apathy_neg_comp1_suprathresh_fc', ...
                       'diff_apathy_gt_quip_neg_comp1_suprathresh_fc', 'interaction_apathyxquip_neg_comp1_suprathresh_fc', ...
                       'main_apathy_neg_comp2_suprathresh_fc'};
        fcLabels = {'Apathy(+) Comp1: Broad hyperconnectivity (N=305, 32e, p=.0057)', ...
                    'Apathy(-) Comp1: Somatomotor-Salience decoupling (N=305, 16e, p=.032)', ...
                    'Diff(-) Comp1: Apathy-dominant Somatomotor-DorsAttn (N=305, 18e, p=.028)', ...
                    'Int(-) Comp1: Visual-DorsAttn-Control (N=305, 20e, p=.026)', ...
                    'Apathy(-) Comp2: DorsAttn-Somatomotor-Salience (N=305, 16e, p=.032, tied second cluster)'};
        fcShort  = {'Ap(+)1: Broad hyper', 'Ap(-)1: SM-Sal', 'Diff(-)1: Ap-dom SM-DA', 'Int(-)1: Vis-DA', 'Ap(-)2: DA-SM-Sal'};
        apathyIdx = [1 2 5]; quipIdx = []; interIdx = 4; diffIdx = 3; sharedIdx = [];
        fcContrast = [1 1 4 3 1];
        idxA = 1; idxB = 2; idxMediation = idxB; idxStructural = diffIdx;
        labelB = 'SomMot-Salience decoupling'; labelB_axis = 'Apathy SomMot-Salience FC (Fisher z)';
    end
    % ---- P305c: as P305 but with the two tied negative-apathy clusters COMBINED
    %      into one feature (union of 32 edges / 30 nodes, features/P305c; the
    %      comp2 column is left unused). 2026-09-08, user request. ----
    if strcmp(EXT_FEATURE_SET, 'P305c')
        meanFCvars  = {'main_apathy_pos_comp1_mean_fc', 'main_apathy_neg_comp1_mean_fc', ...
                       'diff_apathy_gt_quip_neg_comp1_mean_fc', 'interaction_apathyxquip_neg_comp1_mean_fc'};
        supraFCvars = {'main_apathy_pos_comp1_suprathresh_fc', 'main_apathy_neg_comp1_suprathresh_fc', ...
                       'diff_apathy_gt_quip_neg_comp1_suprathresh_fc', 'interaction_apathyxquip_neg_comp1_suprathresh_fc'};
        fcLabels = {'Apathy(+) Comp1: Broad hyperconnectivity (N=305, 32e, p=.0057)', ...
                    'Apathy(-) Combined: clusters 1 + 2 united (N=305, 32e / 30n, each p=.032)', ...
                    'Diff(-) Comp1: Apathy-dominant Somatomotor-DorsAttn (N=305, 18e, p=.028)', ...
                    'Int(-) Comp1: Visual-DorsAttn-Control (N=305, 20e, p=.026)'};
        fcShort  = {'Ap(+)1: Broad hyper', 'Ap(-)c: SM-Sal-DA combined', 'Diff(-)1: Ap-dom SM-DA', 'Int(-)1: Vis-DA'};
        apathyIdx = 1:2; quipIdx = []; interIdx = 4; diffIdx = 3; sharedIdx = [];
        fcContrast = [1 1 4 3];
        idxA = 1; idxB = 2; idxMediation = idxB; idxStructural = diffIdx;
        labelB = 'SomMot-Salience-DorsAttn decoupling (combined)'; labelB_axis = 'Apathy combined hypo FC (Fisher z)';
    end
    % ---- R305: as-reviewed N = 310 component edge sets applied to the 305
    %      (continuity diagnostic; features/R305, corrected covariates). ----
    if strcmp(EXT_FEATURE_SET, 'R305')
        meanFCvars  = {'main_apathy_pos_comp1_mean_fc', 'main_apathy_neg_comp1_mean_fc', ...
                       'diff_apathy_gt_quip_neg_comp1_mean_fc', 'shared_pos_pos_comp1_mean_fc'};
        supraFCvars = {'main_apathy_pos_comp1_suprathresh_fc', 'main_apathy_neg_comp1_suprathresh_fc', ...
                       'diff_apathy_gt_quip_neg_comp1_suprathresh_fc', 'shared_pos_pos_comp1_suprathresh_fc'};
        fcLabels = {'Apathy(+) Comp1: reviewed 310 edge set (51e) on the 305', ...
                    'Apathy(-) Comp1: reviewed 310 edge set (14e) on the 305', ...
                    'Diff(-) Comp1: reviewed 310 edge set (17e) on the 305', ...
                    'Shared(+) Comp1: reviewed 310 edge set (15e) on the 305, not significant at 305'};
        fcShort  = {'Ap(+)1: reviewed edges', 'Ap(-)1: reviewed edges', 'Diff(-)1: reviewed edges', 'Shared(+)1: reviewed edges'};
        apathyIdx = 1:2; quipIdx = []; interIdx = []; diffIdx = 3; sharedIdx = 4;
        fcContrast = [1 1 4 5];
        if ~ismember('shared_pos_pos_comp1_suprathresh_fc', fcSupra.Properties.VariableNames)
            meanFCvars(4) = []; supraFCvars(4) = []; fcLabels(4) = []; fcShort(4) = []; sharedIdx = []; fcContrast = [1 1 4];
        end
        idxA = 1; idxB = 2; idxMediation = idxB; idxStructural = diffIdx;
        labelB = 'RH SomMot-DA decoupling (reviewed edges)'; labelB_axis = 'RH Somatomotor-DorsAttn FC, reviewed edges (Fisher z)';
    end
    % ---- N = 345 EXTENSION (protocol-heterogeneous, corrected FD) ----
    %   components of group_fixed_neither_fdcorr345: hyper 27e p=.013,
    %   differential 22e p=.025, interaction neg 35e p=.013 / pos 18e p=.045;
    %   hypo and shared not significant. If the bridge column (307 Ap(-)
    %   edges applied to the 345, extract_fixed_feature_fc.py) is present it
    %   is appended as an extra apathy component.
    if strcmp(EXT_FEATURE_SET, 'X345')
        hasBridge = ismember('main_apathy_neg_comp1_suprathresh_fc', fcSupra.Properties.VariableNames) && ...
                    ismember('main_apathy_neg_comp1_mean_fc', fcMean.Properties.VariableNames);
        meanFCvars  = {'main_apathy_pos_comp1_mean_fc', 'diff_apathy_gt_quip_neg_comp1_mean_fc', ...
                       'interaction_apathyxquip_neg_comp1_mean_fc', 'interaction_apathyxquip_pos_comp1_mean_fc'};
        supraFCvars = {'main_apathy_pos_comp1_suprathresh_fc', 'diff_apathy_gt_quip_neg_comp1_suprathresh_fc', ...
                       'interaction_apathyxquip_neg_comp1_suprathresh_fc', 'interaction_apathyxquip_pos_comp1_suprathresh_fc'};
        fcLabels = {'Apathy(+) Comp1: Broad hyperconnectivity (N=345, 27e, p=.013)', ...
                    'Diff(-) Comp1: Apathy-dominant Somatomotor-DorsAttn (N=345, 22e, p=.025)', ...
                    'Int(-) Comp1: DorsAttn-Default-SomMot (N=345, 35e, p=.013)', ...
                    'Int(+) Comp1: Default-mode (N=345, 18e, p=.045)'};
        fcShort  = {'Ap(+)1: Broad hyper', 'Diff(-)1: Ap-dom SM-DA', 'Int(-)1: DA-DMN', 'Int(+)1: DMN'};
        apathyIdx = 1; quipIdx = []; interIdx = 3:4; diffIdx = 2; sharedIdx = [];
        fcContrast = [1 4 3 3];
        if hasBridge
            meanFCvars{end+1}  = 'main_apathy_neg_comp1_mean_fc';
            supraFCvars{end+1} = 'main_apathy_neg_comp1_suprathresh_fc';
            fcLabels{end+1} = 'Apathy(-) Comp1: N=307 edges applied to the 345 (bridge, not NBS-significant at 345)';
            fcShort{end+1}  = 'Ap(-)1 307-edges';
            apathyIdx = [1 5]; fcContrast(end+1) = 1;
        end
        idxA = 1; idxB = diffIdx; idxMediation = diffIdx; idxStructural = diffIdx;
        labelB = 'Ap-dom SomMot-DA decoupling'; labelB_axis = 'Apathy-dominant SomMot-DorsAttn FC (Fisher z)';
        if strcmp(EXT_MEDIATOR, 'apneg')
            assert(hasBridge, 'X345 with mediator apneg needs the bridge column (307 edges on 345)');
            idxMediation = 5; idxB = 5;
            labelB = 'SomMot decoupling (307 edges)'; labelB_axis = 'Apathy SomMot FC, 307 edges (Fisher z)';
        end
    end
    % ---- N = 343 EXTENSION under the final rules (2026-09-13; RESULTS 4k):
    %      305 primary + 38 multiband, same 15 blanked parcels, corrected FD.
    %      Only the differential is FWE-significant at 343 (20e p=.026);
    %      hyper 14e p=.053 and interaction 14e p=.061 are above threshold,
    %      hypo lost. The 305 primary component edge sets are carried as
    %      fixed BRIDGE features (features/X343, columns *_p305edges_*), so
    %      the chain keeps the same feature order as P305: 1 hyper (bridge),
    %      2 hypo cluster 1 (bridge), 3 differential (343-defined), 4 the
    %      305 differential edges (bridge), 5 hypo cluster 2 (bridge). ----
    if strcmp(EXT_FEATURE_SET, 'X343')
        need = {'main_apathy_pos_comp1_p305edges_suprathresh_fc', 'main_apathy_neg_comp1_p305edges_suprathresh_fc', ...
                'diff_apathy_gt_quip_neg_comp1_suprathresh_fc', 'diff_apathy_gt_quip_neg_comp1_p305edges_suprathresh_fc', ...
                'main_apathy_neg_comp2_p305edges_suprathresh_fc'};
        assert(all(ismember(need, fcSupra.Properties.VariableNames)), 'X343 needs the 343 differential column and the P305-edge bridge columns (build_features_343.py)');
        meanFCvars  = {'main_apathy_pos_comp1_p305edges_mean_fc', 'main_apathy_neg_comp1_p305edges_mean_fc', ...
                       'diff_apathy_gt_quip_neg_comp1_mean_fc', 'diff_apathy_gt_quip_neg_comp1_p305edges_mean_fc', ...
                       'main_apathy_neg_comp2_p305edges_mean_fc'};
        supraFCvars = need;
        fcLabels = {'Apathy(+) Comp1: 305 hyper edges (32e) on the 343 (bridge; 343 NBS 14e p=.053)', ...
                    'Apathy(-) Comp1: 305 hypo cluster-1 edges (16e) on the 343 (bridge; not significant at 343)', ...
                    'Diff(-) Comp1: Apathy-dominant Somatomotor-DorsAttn (N=343, 20e, p=.026)', ...
                    'Diff(-) 305-edges: the 305 differential edges (18e) on the 343 (bridge)', ...
                    'Apathy(-) Comp2: 305 hypo cluster-2 edges (16e) on the 343 (bridge)'};
        fcShort  = {'Ap(+)1: 305 edges', 'Ap(-)1: 305 edges', 'Diff(-)1: Ap-dom SM-DA', 'Diff(-)1: 305 edges', 'Ap(-)2: 305 edges'};
        apathyIdx = [1 2 5]; quipIdx = []; interIdx = []; diffIdx = 3; sharedIdx = [];
        fcContrast = [1 1 4 4 1];
        idxA = 1; idxB = 2; idxMediation = idxB; idxStructural = diffIdx;
        labelB = 'SomMot-Salience decoupling (305 edges)'; labelB_axis = 'Apathy SomMot-Salience FC, 305 edges (Fisher z)';
    end
    if strcmp(EXT_MEDIATOR, 'diff')
        idxMediation = diffIdx;   % Diff(-)1 as the cognitive mediator for panels a/b
        labelB = 'Ap-dom SomMot-DA decoupling';
        labelB_axis = 'Apathy-dominant SomMot-DorsAttn FC (Fisher z)';
    end
else  % fdThreshold == 0.4
    % MEAN FC (for path b: FC -> cognition)
    meanFCvars = {'main_apathy_pos_comp1_mean_fc', ...           % Visual-Putamen
                  'main_apathy_neg_comp1_mean_fc', ...           % Control-Cerebellum
                  'main_apathy_neg_comp2_mean_fc', ...           % Somatomotor
                  'main_quip_pos_comp1_mean_fc', ...             % Somatomotor
                  'main_quip_neg_comp1_mean_fc', ...             % Hipp-DorsAttn-Visual-Default
                  'main_quip_neg_comp2_mean_fc', ...             % (smaller component)
                  'interaction_apathyxquip_neg_comp1_mean_fc',... % (single neg component)
                  'diff_apathy_gt_quip_pos_comp1_mean_fc', ...   % Apathy-dominant
                  'diff_apathy_gt_quip_neg_comp1_mean_fc', ...   % ICB-dominant
                  'diff_apathy_gt_quip_neg_comp2_mean_fc', ...   % ICB-dominant comp2
                  'shared_pos_pos_comp1_mean_fc', ...            % Shared hyper comp1
                  'shared_pos_pos_comp2_mean_fc', ...            % Shared hyper comp2
                  'shared_pos_pos_comp3_mean_fc', ...            % Shared hyper comp3
                  'shared_pos_neg_comp1_mean_fc'};               % Shared hypo comp1
    % SUPRATHRESHOLD FC (for path a: apathy -> FC)
    supraFCvars = {'main_apathy_pos_comp1_suprathresh_fc', ...   % Visual-Putamen
                   'main_apathy_neg_comp1_suprathresh_fc', ...   % Control-Cerebellum
                   'main_apathy_neg_comp2_suprathresh_fc'};      % Somatomotor
    fcLabels = {'Apathy(+) Comp1: Visual-Putamen', ...
                'Apathy(-) Comp1: Control-Cerebellum', ...
                'Apathy(-) Comp2: Somatomotor', ...
                'ICB(+) Comp1: Somatomotor', ...
                'ICB(-) Comp1: Hipp-DorsAttn-Visual-Default', ...
                'ICB(-) Comp2', ...
                'Interact(-) Comp1', ...
                'Diff(+) Comp1: Apathy-dominant', ...
                'Diff(-) Comp1: ICB-dominant', ...
                'Diff(-) Comp2: ICB-dominant', ...
                'Shared(+) Comp1', ...
                'Shared(+) Comp2', ...
                'Shared(+) Comp3', ...
                'Shared(-) Comp1'};
    fcShort = {'Ap(+)1: Vis-Pu','Ap(-)1: Ctrl-Cereb','Ap(-)2: SM', ...
               'ICB(+)1: SM', ...
               'ICB(-)1: Hipp-DA','ICB(-)2', ...
               'Int(-)1', ...
               'Diff(+)1: Ap-dom','Diff(-)1: ICB-dom','Diff(-)2', ...
               'Sh(+)1','Sh(+)2','Sh(+)3','Sh(-)1'};
    apathyIdx  = 1:3;
    quipIdx    = 4:6;
    interIdx   = 7;
    diffIdx    = 8:10;
    sharedIdx  = 11:14;
    fcContrast = [1 1 1 2 2 2 3 4 4 4 5 5 5 5];
    % Panel mapping (all indices are fi = index into meanFCvars)
    idxA = 2;  % Control-Cerebellum (Panel A scatter — primary / mediation)
    idxB = 3;  % Somatomotor        (Panel B scatter -> Supp)
    idxMediation = idxA;  % FC component used as mediator in Section 6
    idxStructural = diffIdx(1);  % first diff component — placeholder; verify if 0.4 used
    labelA = 'Ctrl-Cereb';
    labelA_axis = 'Control-Cerebellar FC (Fisher z)';
    labelB = 'Somatomotor';
    labelB_axis = 'Somatomotor FC (Fisher z)';
    % Panel C: domain-specific scatter
    idxC_fc  = 2;  % Control-Cereb component index in meanFCvars
    idxC_dom = 1;  % Executive domain index in cogLabels (update as needed)
    labelC = 'Ctrl-Cereb';
    labelC_axis = 'Control-Cerebellar FC (Fisher z)';
    labelC_dom  = 'Executive';
end
nAllFC     = numel(meanFCvars);
nApFC      = numel(apathyIdx);
% --- supraIdx: which fi indices have suprathreshold columns available ---
%     Used for path-a / DAT path-a / their FDR families. Distinct from
%     apathyIdx (which marks contrasts that are *semantically* apathy).
%     Post-2026-04-16: at FD=0.3 all 4 sig components have supra columns;
%     FD=0.4 still only has supra for the first 3 (main_apathy comps).
if fdThreshold == 0.3
    supraIdx = 1:numel(supraFCvars);   % 4 (C1) or 6 (C2)
else
    supraIdx = 1:3;
end
nSupraFC   = numel(supraIdx);
assert(numel(supraFCvars) == nSupraFC, ...
    'supraFCvars length must equal numel(supraIdx).');
idxDAT_fi  = idxA;               % Primary component index for DAT scatter (fi)

% --- Build fcData matrix (nFC x nAllFC) — one column per FC component. ---
%     useSupraForPathB flag selects supra (any component in supraIdx) vs
%     mean FC. All downstream path-b analyses use fcData(:, fi).
fcData    = nan(height(fcMean), nAllFC);
fcSource  = repmat({'mean'}, 1, nAllFC);
for fi = 1:nAllFC
    ci_supra = find(supraIdx == fi, 1);
    if useSupraForPathB && ~isempty(ci_supra) && ci_supra <= numel(supraFCvars)
        fcData(:, fi) = fcSupra.(supraFCvars{ci_supra});
        fcSource{fi}  = 'supra';
    else
        fcData(:, fi) = fcMean.(meanFCvars{fi});
        fcSource{fi}  = 'mean';
    end
end
fprintf('  FC source per component: %s\n', strjoin(fcSource, ', '));

%% ========================================================================
%  2.  EXTRACT SLOPES — 4 DOMAINS × REDUCED LME (+ OLS fallback)
% =========================================================================
fprintf('\n[2] Extracting cognitive slopes (reduced LME)...\n');

longPATNO = double(Tlong.PATNO);
fcPATNO   = double(fcMean.subject_id);
slopeAll  = nan(nFC, nDom);  % reduced LME slopes

% Covariates
covCands = {'age','sex','educ','nmiss','duration','ledd','updrs3_score','depfactor'};
covPresent = {};
for ci=1:numel(covCands)
    if ismember(covCands{ci},Tlong.Properties.VariableNames), covPresent{end+1}=covCands{ci}; end %#ok
end
covStr=''; if ~isempty(covPresent), covStr=[' + ' strjoin(covPresent,' + ')]; end

okBase = ~isnan(double(Tlong.(timeVar)));
for ci=1:numel(covPresent)
    vn=covPresent{ci};
    if iscategorical(Tlong.(vn)), okBase=okBase&~isundefined(Tlong.(vn));
    else, okBase=okBase&~isnan(double(safeCol(Tlong,vn))); end
end

for di = 1:nDom
    if isempty(cogVarNames{di}), continue; end
    okD = okBase & ~isnan(double(safeCol(Tlong,cogVarNames{di})));
    TlmeD = Tlong(okD,:);
    fmla = sprintf('%s ~ %s%s + (1 + %s | PATNO)', cogVarNames{di}, timeVar, covStr, timeVar);
    try
        mdlD = fitlme(TlmeD, fmla, 'FitMethod','REML');
        [slopes,nM] = extractBLUPsFromModel(mdlD, fcPATNO, nFC);
        slopeAll(:,di) = slopes;
        fprintf('    %-15s  AIC=%7.1f  matched=%d\n', cogLabels{di}, mdlD.ModelCriterion.AIC, nM);
    catch ME
        fprintf('    %-15s  FAILED: %s\n', cogLabels{di}, ME.message);
    end
end

% OLS fallback
for si=1:nFC
    sIdx=longPATNO==fcPATNO(si); if sum(sIdx)<2, continue; end
    tSub=Tlong(sIdx,:); tVals=double(tSub.(timeVar));
    for di=1:nDom
        if ~isnan(slopeAll(si,di))||isempty(cogVarNames{di}), continue; end
        yV=double(safeCol(tSub,cogVarNames{di})); ok=~isnan(tVals)&~isnan(yV);
        if sum(ok)>=2, b=[ones(sum(ok),1),tVals(ok)]\yV(ok); slopeAll(si,di)=b(2); end
    end
end
fprintf('    Slopes (raw): %d/%d subjects\n', sum(~isnan(slopeAll(:,1))), nFC);

% Primary = raw slopes (NoWinsor). slopeAllRaw kept as alias for compatibility.
slopeAllRaw = slopeAll;

% --- Sensitivity-only: 1st/99th percentile winsorised copy for §14 -------
fprintf('    Computing winsorised slope copy (sensitivity only) at 1st/99th percentile...\n');
slopeAllWinsor = slopeAll;
for di = 1:nDom
    v = slopeAllWinsor(:,di);
    ok = ~isnan(v);
    if sum(ok) < 10, continue; end
    p1  = prctile(v(ok), 1);
    p99 = prctile(v(ok), 99);
    nClip = sum(v(ok)<p1 | v(ok)>p99);
    v(ok & v<p1)  = p1;
    v(ok & v>p99) = p99;
    slopeAllWinsor(:,di) = v;
    if nClip > 0
        fprintf('      %-15s  clipped %d values (1st=%.4f, 99th=%.4f)\n', ...
            cogLabels{di}, nClip, p1, p99);
    end
end

fprintf('    Slopes (final, primary=raw): %d/%d subjects\n', sum(~isnan(slopeAll(:,1))), nFC);

% --- Per-subject baseline cognition values (for partial-correlation control) ---
%     Empirical check showed BLUP slopes still correlate strongly with
%     baseline values even though the LME random-effects covariance is
%     supposed to handle this:
%        Corr(BLUP cog slope, baseline bifactor_g):  +0.47
%        Corr(BLUP DAT slope, baseline DAT_mean):    -0.78
%        Corr(BLUP apathy slope, baseline apathy):   +0.23
%     With sparse longitudinal data the random-effects covariance is
%     under-determined and the BLUP shrinkage doesn't decouple level
%     from rate cleanly. We therefore add the relevant baseline as an
%     explicit covariate whenever the outcome is a BLUP slope.
fprintf('    Extracting per-subject baseline cognition values...\n');
baselineCog = nan(nFC, nDom);
longPATNO = double(Tlong.PATNO);
fcPATNO   = double(fcMean.subject_id);
tValAll   = double(Tlong.(timeVar));
for di = 1:nDom
    if isempty(cogVarNames{di}), continue; end
    yV = double(safeCol(Tlong, cogVarNames{di}));
    for si = 1:nFC
        rows = (longPATNO == fcPATNO(si)) & ~isnan(yV) & ~isnan(tValAll);
        if ~any(rows), continue; end
        yRows = yV(rows);  tRows = tValAll(rows);
        [~, ord] = sort(tRows);
        yRows = yRows(ord);
        baselineCog(si,di) = yRows(1);   % first available visit
    end
end
baselineG     = baselineCog(:,1);
baselineExec  = baselineCog(:,2);
baselineMem   = baselineCog(:,3);
baselineVisuo = baselineCog(:,4);
fprintf('      bifactor_g baseline:   %d/%d non-NaN\n', sum(~isnan(baselineG)), nFC);
fprintf('      obs_att_exec baseline: %d/%d non-NaN\n', sum(~isnan(baselineExec)), nFC);
fprintf('      obs_memory baseline:   %d/%d non-NaN\n', sum(~isnan(baselineMem)), nFC);
fprintf('      obs_visuo baseline:    %d/%d non-NaN\n', sum(~isnan(baselineVisuo)), nFC);

% --- Per-subject baseline clinical covariates (for sensitivity analysis) -
%     Pull baseline (min-time row per subject). Column names vary across
%     tables, so try common variants and report which was found.

% Define: display name -> candidate column names in Tlong / TpdSub
covCandsMap = { ...
    'duration',  {'duration','disease_duration','PD_duration','pd_duration'}; ...
    'LEDD',      {'ledd','LEDD','ledd_baseline'}; ...
    'UPDRS3',    {'updrs3_score','UPDRS3','updrs3','mds_updrs3'}; ...
    'QUIP',      {'quip','QUIP','quip_score','QUIP_score','quip_total'}; ...
    'depfactor', {'depfactor','depression_factor','depfctr','GDS_depfactor','dep_score'}; ...
    'educ',      {'educ','EDUCYRS','education','education_years','EDUC'}};

tableCandidates = {'Tlong'};
if exist('TpdSub','var'), tableCandidates{end+1} = 'TpdSub'; end
fprintf('    Looking up baseline covariates in: %s\n', strjoin(tableCandidates,', '));

baselineDur     = nan(nFC,1);
baselineLEDD    = nan(nFC,1);
baselineUPDRS3  = nan(nFC,1);
baselineQUIP    = nan(nFC,1);
baselineDEP     = nan(nFC,1);
baselineEduc    = nan(nFC,1);
storeMap = {'baselineDur','baselineLEDD','baselineUPDRS3','baselineQUIP','baselineDEP','baselineEduc'};

for ri = 1:size(covCandsMap,1)
    displayName = covCandsMap{ri,1};
    candNames   = covCandsMap{ri,2};
    found = false;
    for tableName = tableCandidates
        T = eval(tableName{1});
        if ~istable(T), continue; end
        if ~ismember('PATNO', T.Properties.VariableNames), continue; end
        for cn = candNames
            if ismember(cn{1}, T.Properties.VariableNames)
                % Extract baseline value per subject from T
                tmp = nan(nFC,1);
                tPATNO = double(T.PATNO);
                if ismember(timeVar, T.Properties.VariableNames)
                    tYr = double(T.(timeVar));
                elseif ismember('YEAR', T.Properties.VariableNames)
                    tYr = double(T.YEAR);
                else
                    tYr = zeros(height(T),1);
                end
                tVals = double(safeCol(T, cn{1}));
                for si = 1:nFC
                    pidx = tPATNO == fcPATNO(si);
                    if ~any(pidx), continue; end
                    yrS = tYr(pidx); valS = tVals(pidx);
                    [~,ord] = sort(yrS);
                    valS = valS(ord);
                    k = find(~isnan(valS),1);
                    if ~isempty(k), tmp(si) = valS(k); end
                end
                if sum(~isnan(tmp)) > sum(~isnan(eval(storeMap{ri})))
                    eval(sprintf('%s = tmp;', storeMap{ri}));
                    fprintf('      %-10s found as ''%s'' in %s  (n=%d non-NaN)\n', ...
                        displayName, cn{1}, tableName{1}, sum(~isnan(tmp)));
                    found = true;
                    break;
                end
            end
        end
        if found, break; end
    end
    if ~found
        fprintf('      %-10s NOT FOUND (tried: %s)\n', displayName, ...
            strjoin(candNames,', '));
    end
end

fprintf('    Baseline covariates: duration n=%d, LEDD n=%d, UPDRS-3 n=%d, QUIP n=%d, dep n=%d, educ n=%d\n', ...
    sum(~isnan(baselineDur)), sum(~isnan(baselineLEDD)), sum(~isnan(baselineUPDRS3)), ...
    sum(~isnan(baselineQUIP)), sum(~isnan(baselineDEP)), sum(~isnan(baselineEduc)));

% ====================================================================== %
%   PRIMARY ANALYSIS SPECIFICATION                                       %
% ---------------------------------------------------------------------- %
%   All primary subject-level FC partial-correlation tests in Sections   %
%   3, 4, 5, 5b (iii), and 10 use:                                       %
%                                                                        %
%     - Spearman partial r        (robust to monotonic non-linearity     %
%                                   and slope outliers)                  %
%     - Minimal covariates        FD + TIV                               %
%       Aligns the primary partial-correlation framework with the        %
%       Section 6 mediation framework (which has always used FD+TIV +    %
%       OLS); reviewers reasonably expected one consistent primary       %
%       covariate set across §3–§10 and §6.                              %
%                                                                        %
%   Sensitivity analyses (Section 14) report the broader NBS-matched     %
%   covariate set (FD + TIV + age + sex + QUIP + dep + apathy×QUIP)      %
%   plus rank-based / ExtCovars variants for robustness.                 %
% ====================================================================== %
sexNum = sexFC;                                           % M=1 / F=0
primaryCovs   = [fcMean.mean_FD, fcMean.TIV];
primaryCovLbl = {'FD','TIV'};
primaryMethod = 'Spearman';   % Capitalised to match MATLAB partialcorr's 'Type' values
nPrimaryComplete = sum(all(~isnan(primaryCovs),2));
fprintf('\n  Primary: %s partial r | %s\n', primaryMethod, strjoin(primaryCovLbl,', '));
fprintf('  Subjects complete on primary covariates: %d / %d\n', nPrimaryComplete, nFC);

% --- Merge DAT values from baseline table --------------------------------
fprintf('\n    Merging DAT values from bl table...\n');
fcMean.DAT_caudate = nan(nFC,1);
fcMean.DAT_putamen = nan(nFC,1);
if exist('bl','var') && ismember('PATNO', bl.Properties.VariableNames)
    blPATNO = double(bl.PATNO);
    for si = 1:nFC
        row = find(blPATNO == fcPATNO(si), 1);
        if ~isempty(row)
            if ismember('MIA_CAUDATE_BILAT', bl.Properties.VariableNames)
                fcMean.DAT_caudate(si) = double(bl.MIA_CAUDATE_BILAT(row));
            end
            if ismember('MIA_PUTAMEN_BILAT', bl.Properties.VariableNames)
                fcMean.DAT_putamen(si) = double(bl.MIA_PUTAMEN_BILAT(row));
            end
        end
    end
    nDAT = sum(~isnan(fcMean.DAT_caudate));
    fprintf('    DAT matched: %d / %d (caudate: %.2f +/- %.2f)\n', ...
        nDAT, nFC, nanmean(fcMean.DAT_caudate), nanstd(fcMean.DAT_caudate));
    hasDAT = nDAT >= 20;
else
    fprintf('    WARNING: bl table not found or missing PATNO. DAT analyses skipped.\n');
    hasDAT = false;
end

%% ========================================================================
%  3.  PATH a — APATHY → FC (suprathreshold, NBS confirmation)
%      pathA indexed by fi (1..nAllFC); supraIdx entries populated.
%      supraFCvars must be in the same order as supraIdx.
%      FDR family: across all components in supraIdx (pP_FDR field).
% =========================================================================
fprintf('\n[3] Path a: Apathy -> FC (suprathreshold)...\n');

emptyPathA = struct('rP',NaN,'pP',NaN,'pP_FDR',NaN,'beta',NaN,'beta_p',NaN, ...
                    'n',NaN,'x',[],'y',[],'mdl',[]);
pathA = repmat(emptyPathA, 1, nAllFC);

fprintf('    Primary: %s partial r | %s\n', primaryMethod, strjoin(primaryCovLbl,', '));
fprintf('    FDR family: %d components with suprathreshold FC (supraIdx)\n', nSupraFC);
for ci = 1:nSupraFC
    fi  = supraIdx(ci);
    fcV = fcSupra.(supraFCvars{ci});
    apV = fcMean.apathy_composite;

    % Primary: Spearman partial r | FD + TIV (matches primary spec at L442)
    covsPathA = [fcMean.mean_FD, fcMean.TIV];
    [rP_prim, pP_prim] = partialcorr(apV, fcV, covsPathA, 'Type', primaryMethod, 'Rows', 'complete');
    nP_prim = sum(~isnan(apV) & ~isnan(fcV) & all(~isnan(covsPathA),2));
    % Sensitivity: Pearson partial r | FD+TIV (matches initial minimal analysis)
    ok = ~isnan(fcV)&~isnan(apV)&~isnan(fcMean.mean_FD)&~isnan(fcMean.TIV);
    rP_basic = partialcorr(fcV(ok), apV(ok), [fcMean.mean_FD(ok), fcMean.TIV(ok)], 'Rows','complete');
    nP_basic = sum(ok); df=nP_basic-4;
    tS = rP_basic*sqrt(df/(1-rP_basic^2));
    pP_basic = 2*(1-tcdf(abs(tS),df));

    mdl = fitlm([apV(ok),fcMean.mean_FD(ok),fcMean.TIV(ok)], fcV(ok));

    % Primary stored in rP, pP, n
    pathA(fi).rP = rP_prim;  pathA(fi).pP = pP_prim;  pathA(fi).n = nP_prim;
    % Sensitivity for record-keeping
    pathA(fi).rP_basic = rP_basic;  pathA(fi).pP_basic = pP_basic;
    pathA(fi).beta = mdl.Coefficients.Estimate(2);
    pathA(fi).beta_p = mdl.Coefficients.pValue(2);
    pathA(fi).x = apV(ok); pathA(fi).y = fcV(ok); pathA(fi).mdl = mdl;
end

% --- Path-a FDR across supraIdx family ---
pathA_pVec = [pathA(supraIdx).pP];
pathA_fdr  = fdr_bh(pathA_pVec);
for k = 1:nSupraFC, pathA(supraIdx(k)).pP_FDR = pathA_fdr(k); end

% Print with FDR
fprintf('    %-45s  %8s  %8s  %8s   [%-18s]\n', '', 'rho', 'p_unc', 'p_FDR', 'basic Pearson r, p');
fprintf('    %s\n', repmat('-',1,100));
for ci = 1:nSupraFC
    fi = supraIdx(ci);
    fprintf('    %-45s  %+7.3f  %.4f  %.4f %s   [%+5.3f, %.4f]\n', ...
        fcLabels{fi}, pathA(fi).rP, pathA(fi).pP, pathA(fi).pP_FDR, ...
        pToStars(pathA(fi).pP_FDR), pathA(fi).rP_basic, pathA(fi).pP_basic);
end

%% ========================================================================
%  4.  PATH b — FC × GLOBAL COGNITION SLOPE
%      Primary = Spearman partial r | FD + TIV.
%      Sensitivity values (Pearson | FD+TIV) retained in compR.rP_basic
%      for direct comparison with conventional scatter residualization.
% =========================================================================
fprintf('\n[4] Path b: FC -> Global cognition slope\n');
fprintf('    Primary: %s partial r | %s\n', primaryMethod, strjoin(primaryCovLbl,', '));

% --- Primary FDR family ---
%     With useSupraForPathB=true, the primary family is the set of
%     components actually using suprathreshold FC (now all 4 at FD=0.3,
%     post-2026-04-16 supra-CSV fix; previously was apathyIdx-only).
if useSupraForPathB
    primaryIdx = supraIdx;
else
    primaryIdx = 1:nAllFC;
end

% Pre-allocate compR with all fields, length = nAllFC (indexed by fi)
emptyEntry = struct('r',NaN,'p',NaN,'rP',NaN,'pP',NaN,'pFDR',NaN, ...
                    'n',NaN,'nP',NaN,'x',[],'y',[],'mdl',[], ...
                    'rP_DAT',NaN,'pP_DAT',NaN, ...
                    'rP_basic',NaN,'pP_basic',NaN);
compR = repmat(emptyEntry, 1, nAllFC);

for fi = primaryIdx
    fcV = fcData(:, fi);
    sV  = slopeAll(:,1);  % global cognition slope

    % Primary: Spearman partial r | FD + TIV
    % Outcome is global cog slope → control for baseline bifactor_g to
    % remove residual baseline-rate coupling the LME RE-cov didn't decouple.
    covs_b = [primaryCovs];
    [rP_prim, pP_prim] = partialcorr(fcV, sV, covs_b, 'Type', primaryMethod, 'Rows', 'complete');
    n_prim = sum(~isnan(fcV) & ~isnan(sV) & all(~isnan(covs_b),2));
    % Sensitivity: Pearson partial r | FD + TIV (kept for scatter residualization)
    [r_raw, p_raw, rP_basic, pP_basic, n_basic, nP_basic, x_raw, y_raw, mdl_raw] = ...
        fcSlopeCorr(fcV, sV, fcMean.mean_FD, fcMean.TIV);

    compR(fi).rP  = rP_prim;    compR(fi).pP  = pP_prim;    compR(fi).nP = n_prim;
    compR(fi).r   = r_raw;      compR(fi).p   = p_raw;      compR(fi).n  = n_basic;
    compR(fi).rP_basic = rP_basic;  compR(fi).pP_basic = pP_basic;
    compR(fi).x   = x_raw;      compR(fi).y   = y_raw;      compR(fi).mdl = mdl_raw;
end

% FDR across primary family (on primary Spearman NBS-matched p-values)
rawP = [compR(primaryIdx).pP];
fdrP = fdr_bh(rawP);
for k = 1:numel(primaryIdx)
    compR(primaryIdx(k)).pFDR = fdrP(k);
end

fprintf('\n  Primary (FDR across %d component(s)):\n', numel(primaryIdx));
fprintf('  %-45s  %9s  %9s  %9s   [%-18s]\n', ...
    '', 'rho', 'p_unc', 'p_FDR', 'basic Pearson r, p');
fprintf('  %s\n', repmat('-',1,100));
for fi = primaryIdx
    fprintf('  %-45s  %+8.3f  %.4f  %.4f %s   [%+5.3f, %.4f]\n', fcLabels{fi}, ...
        compR(fi).rP, compR(fi).pP, compR(fi).pFDR, pToStars(compR(fi).pFDR), ...
        compR(fi).rP_basic, compR(fi).pP_basic);
end

% sigGlobal now contains fi indices directly (not positions into primaryIdx)
sigGlobal = primaryIdx(fdrP < 0.05);
fprintf('  FDR-significant: %d component(s)\n', numel(sigGlobal));

%% ========================================================================
%  5.  DOMAIN SPECIFICITY  — moved to Figure4_Supp_Domain.m
%      The full FC × subdomain heatmap analysis (4 variants of LME × FC
%      partialcorr, family-wise FDR-BH, 2×2 figure) now lives in
%      Figure4_Supp_Domain.m. Run it after this script's §1–§3 finish.
%      The variables rMat, pMat, pFDR_domain, domScatter that this section
%      used to produce were not consumed by any plot in this file (only
%      the SuppFig S1D heatmap, which is also moved).
% =========================================================================

%% ========================================================================
%  5c. NBS STRUCTURAL-CONNECTOME METRICS × Global Cog Slope (STANDALONE)
%      Source: cognition_ready.tsv — per-subject NBS metrics + covariates.
%      Cohort = TSV PATNOs (NOT restricted to FC).
%      Outcome: cog-slope LME BLUP re-extracted on full Tlong.
%      Covariates (all from cognition_ready.tsv):
%          TIV  +  dwi_motion  +  site (dummies)
%        + whole_conn_excl_apathy_t**                   (if column present)
%        + hubvol_<LABEL>                              (Hubs family only)
%      age/sex/educ/duration absorbed by the LME, not re-added (matches §4).
%      FDR-BH within each variable family.
% =========================================================================
fprintf('\n[5c] NBS structural-connectome metrics x cog slope (STANDALONE)\n');

cogFile = fullfile(getenv('PPMI_IMAGING_ROOT'), 'MRI_Analysis/cognition_ready.tsv');
if isfile(cogFile)
    Tcog = readtable(cogFile,'FileType','text','Delimiter','\t','TreatAsMissing','NA');
    acmPATNO = double(Tcog.PATNO); nAcm = numel(acmPATNO);
    fprintf('  Cohort: %d subjects from cognition_ready.tsv.\n', nAcm);

    % --- Threshold selector ('31' or '35') ------------------------------
    acmThreshold = '31';
    if ~exist('acmThreshold','var') || isempty(acmThreshold), acmThreshold = '35'; end
    assert(any(strcmp(acmThreshold,{'31','35'})), 'acmThreshold must be ''31'' or ''35''.');
    sfx = ['_t' acmThreshold];
    fprintf('  Threshold: t = %s.%s   (suffix %s)\n', acmThreshold(1), acmThreshold(2), sfx);

    famNames   = {'Primary','Hubs','NetPairs','Canonical'};
    primStems  = {'circuit_strength'};
    hubStems   = {'hub_L_IsthmusCing','hub_L_CaudAntCing','hub_R_Putamen','hub_R_SupFrontal', ...
                  'hub_L_InfParietal','hub_L_Precuneus','hub_R_Hippocampus','hub_R_LatOrbFront', ...
                  'hub_R_Insula','hub_R_RostMidFront'};
    canonStems = {'circuit_frontostriatal','circuit_dmn','circuit_limbic','circuit_salience'};
    if strcmp(acmThreshold,'31')
        netStems = {'net_Control_Control','net_Control_Default','net_Control_SalVentAttn', ...
                    'net_Control_Somatomotor','net_Control_Subcortical','net_Default_Default', ...
                    'net_Default_DorsAttn','net_Default_Limbic','net_Default_SalVentAttn', ...
                    'net_Default_Subcortical','net_Default_Visual','net_SalVentAttn_Subcortical', ...
                    'net_Subcortical_Subcortical','net_Subcortical_Visual','net_Visual_Visual'};
    else  % '35'
        netStems = {'net_Control_Control','net_Control_Default','net_Control_SalVentAttn', ...
                    'net_Default_Default','net_Default_DorsAttn','net_Default_Limbic', ...
                    'net_Default_SalVentAttn','net_Default_Subcortical','net_Default_Visual', ...
                    'net_Limbic_Subcortical','net_Subcortical_Subcortical'};
    end
    appendSfx = @(c) cellfun(@(s)[s sfx], c, 'UniformOutput', false);
    fams = {appendSfx(primStems), appendSfx(hubStems), appendSfx(netStems), appendSfx(canonStems)};

    % --- Cog slope: LME BLUP re-fit on full Tlong, extracted for cohort --
    okBaseL = ~isnan(double(Tlong.(timeVar))) & ~isnan(double(safeCol(Tlong, cogVarNames{1})));
    for ci=1:numel(covPresent)
        vn=covPresent{ci};
        if iscategorical(Tlong.(vn)), okBaseL=okBaseL&~isundefined(Tlong.(vn));
        else, okBaseL=okBaseL&~isnan(double(safeCol(Tlong,vn))); end
    end
    fmla_acm = sprintf('%s ~ %s%s + (1 + %s | PATNO)', cogVarNames{1}, timeVar, covStr, timeVar);
    mdlAcm = fitlme(Tlong(okBaseL,:), fmla_acm, 'FitMethod','REML');
    [sV, nBlup] = extractBLUPsFromModel(mdlAcm, acmPATNO, nAcm);
    fprintf('  Cog BLUP: %d/%d cohort subjects with longitudinal cog data.\n', nBlup, nAcm);

    % --- Base covariates from cognition_ready.tsv ------------------------
    TIVv = double(Tcog.TIV);
    dwi  = double(Tcog.dwi_motion);
    siteRaw = Tcog.site;
    if isnumeric(siteRaw), siteCat = categorical(siteRaw);
    else,                  siteCat = categorical(string(siteRaw));
    end
    siteDum = dummyvar(siteCat);
    if size(siteDum,2) > 1, siteDum(:,1) = []; end   % drop reference level

    wcEx_col = sprintf('whole_connectome_mean_excl_apathy_t%s', acmThreshold);
    if ismember(wcEx_col, Tcog.Properties.VariableNames)
        Z_base = [TIVv, dwi, siteDum, double(Tcog.(wcEx_col))];
        wcStr  = ' + whole_conn_excl_apathy';
    else
        Z_base = [TIVv, dwi, siteDum];
        wcStr  = '';
        fprintf('  NOTE: %s not in TSV; running without circuit-specificity covariate.\n', wcEx_col);
    end
    fprintf('  Full Z: TIV + dwi_motion + site_dummies(k=%d)%s [+ hubvol for Hubs].\n', ...
        size(siteDum,2), wcStr);

    % --- Minimal Z (sensitivity reference): TIV + dwi_motion only --------
    Z_min = [TIVv, dwi];
    fprintf('  Min  Z: TIV + dwi_motion (sensitivity reference).\n');
    fprintf('  Complete-case ceilings: full %d/%d   min %d/%d.\n', ...
        sum(~isnan(sV) & all(~isnan(Z_base),2)), nAcm, ...
        sum(~isnan(sV) & all(~isnan(Z_min ),2)), nAcm);

    % --- Per-family Spearman partial + FDR (both minimal and full Z) -----
    acmRows = {};
    for gi = 1:numel(fams)
        vars = fams{gi}; nV = numel(vars);
        rhoM  = nan(nV,1); pM  = nan(nV,1); nM = zeros(nV,1);
        rhoF  = nan(nV,1); pF  = nan(nV,1); nF = zeros(nV,1);
        isHubs = strcmp(famNames{gi}, 'Hubs');
        for vi = 1:nV
            varname = vars{vi};
            if ~ismember(varname, Tcog.Properties.VariableNames), continue; end
            xV = double(Tcog.(varname));

            % Minimal Z (same for all variables; hub volume NOT appended)
            okM = ~isnan(xV) & ~isnan(sV) & all(~isnan(Z_min),2);
            nM(vi) = sum(okM);
            if nM(vi) >= 10
                [rhoM(vi), pM(vi)] = partialcorr(xV(okM), sV(okM), Z_min(okM,:), 'Type','Spearman');
            end

            % Full Z (base + per-hub volume for Hubs family)
            Z_full = Z_base;
            if isHubs
                hv = regexprep(varname, '^hub_', 'hubvol_');
                hv = regexprep(hv, '_t\d+$', '');
                if ismember(hv, Tcog.Properties.VariableNames)
                    Z_full = [Z_full, double(Tcog.(hv))]; %#ok
                end
            end
            okF = ~isnan(xV) & ~isnan(sV) & all(~isnan(Z_full),2);
            nF(vi) = sum(okF);
            if nF(vi) >= 10
                [rhoF(vi), pF(vi)] = partialcorr(xV(okF), sV(okF), Z_full(okF,:), 'Type','Spearman');
            end
        end
        % FDR-BH within family for each Z separately (NaN-safe)
        pFDR_M = nan(nV,1);  vM = ~isnan(pM); if any(vM), pFDR_M(vM) = fdr_bh(pM(vM)); end
        pFDR_F = nan(nV,1);  vF = ~isnan(pF); if any(vF), pFDR_F(vF) = fdr_bh(pF(vF)); end

        fprintf('\n  %s family (n_var=%d, FDR within each Z):\n', famNames{gi}, nV);
        fprintf('  %-32s   %7s %7s %8s    %7s %7s %8s    %s\n', ...
            'Variable','rho_min','p_min','pFDR_min','rho_full','p_full','pFDR_full','n(min/full)');
        for vi = 1:nV
            if isnan(rhoM(vi)),  s1='     —'; s1p='     —'; s2='     —'; starM=' ';
            else, s1=sprintf('%+7.3f',rhoM(vi)); s1p=sprintf('%7.4f',pM(vi)); s2=sprintf('%7.4f',pFDR_M(vi)); starM=pToStars(pFDR_M(vi)); end
            if isnan(rhoF(vi)),  s3='     —'; s3p='     —'; s4='     —'; starF=' ';
            else, s3=sprintf('%+7.3f',rhoF(vi)); s3p=sprintf('%7.4f',pF(vi)); s4=sprintf('%7.4f',pFDR_F(vi)); starF=pToStars(pFDR_F(vi)); end
            fprintf('  %-32s   %s %s %s%s   %s %s %s%s   %d/%d\n', ...
                vars{vi}, s1, s1p, s2, starM, s3, s3p, s4, starF, nM(vi), nF(vi));
            acmRows(end+1,:) = {famNames{gi}, vars{vi}, ...
                rhoM(vi), pM(vi), pFDR_M(vi), nM(vi), ...
                rhoF(vi), pF(vi), pFDR_F(vi), nF(vi)}; %#ok
        end
    end
    acmOutCsv = sprintf('Table_Fig4_ApathyCircuit_x_CogSlope_t%s.csv', acmThreshold);
    writetable(cell2table(acmRows,'VariableNames', ...
        {'Family','Variable', ...
         'rho_min','p_min','pFDR_min','n_min', ...
         'rho_full','p_full','pFDR_full','n_full'}), ...
        fullfile(outputDir, acmOutCsv));
    fprintf('\n  -> %s saved (rho_min + rho_full per row).\n', acmOutCsv);
else
    fprintf('  %s not found; skipping.\n', cogFile);
end

%% ========================================================================
%  5d. FC × COG SLOPE — UPWARD-MATCHED COVARIATES (sensitivity vs §5c)
%      Mirrors §5c full Z so FC and NBS analyses are stringency-matched.
%      Min:  mean_FD + TIV                              (= §4 current primary)
%      Full: mean_FD + TIV + site_dummies + brain_atrophy_index
%                          + whole_FC_mean_excl_apathy_t** (if column present)
%      brain_atrophy_index (and whole_FC) come from cognition_ready.tsv —
%      full-Z sample drops to FC ∩ Tcog (~280); min-Z keeps all 310.
%      FDR-BH across all FC components for each Z separately.
% =========================================================================
fprintf('\n[5d] FC x cog slope under upward-matched covariates (sensitivity)\n');
% Site for FC pulled from PPMI curated baseline (covers all FC subjects;
% cognition_ready.tsv only has 280/310, which would drop the sample).
if ~exist('Tcurated','var') || ~istable(Tcurated)
    curFile = fullfile(getenv('PPMI_CLINICAL_ROOT'), 'ppmi_data/PPMI_Curated_Data_Cut_Public_20260223.xlsx');
    Tcurated = readtable(curFile, 'Sheet','20260105');
end
[~, uIdx_d] = unique(Tcurated.PATNO, 'first');
Tbl_d = Tcurated(uIdx_d, :);
if iscell(Tbl_d.PATNO)||isstring(Tbl_d.PATNO), cPATNO_d = str2double(string(Tbl_d.PATNO));
else, cPATNO_d = double(Tbl_d.PATNO); end
[tfFC, locFC] = ismember(fcPATNO, cPATNO_d);
fprintf('  PATNO match (site source = Tcurated): %d/%d FC subjects.\n', sum(tfFC), nFC);

siteFC_raw = nan(nFC,1);
siteFC_raw(tfFC) = double(Tbl_d.SITE(locFC(tfFC)));
% Categorical w/ NaN handled: subjects with missing site get all-zero dummy
% rows (implicitly reference) — still kept in the analysis.
siteCat_d = categorical(siteFC_raw);
siteDumFC = dummyvar(siteCat_d);
if size(siteDumFC,2) > 1, siteDumFC(:,1) = []; end   % drop reference level

% --- brain_atrophy_index + whole_FC_mean_excl_apathy (from Tcog) --------
%   These live in cognition_ready.tsv; FC subjects not in Tcog get NaN
%   (drops the full-Z sample, but min-Z still uses all 310).
baiFC = nan(nFC,1);
if exist('Tcog','var') && istable(Tcog)
    [tfCog, locCog] = ismember(fcPATNO, double(Tcog.PATNO));
    if ismember('brain_atrophy_index', Tcog.Properties.VariableNames)
        baiFC(tfCog) = double(Tcog.brain_atrophy_index(locCog(tfCog)));
    end
end
wFCex = nan(nFC,1); wFCex_src = '';
wFCCands = {sprintf('whole_FC_mean_excl_apathy_t%s', acmThreshold), ...
            'whole_FC_mean_excl_apathy'};
% Prefer fcSupra (row-aligned with fcMean, 310 subjects). Fall back to Tcog.
if exist('fcSupra','var') && istable(fcSupra)
    for ci = 1:numel(wFCCands)
        if ismember(wFCCands{ci}, fcSupra.Properties.VariableNames)
            wFCex = double(fcSupra.(wFCCands{ci}));
            wFCex_src = [wFCCands{ci} ' (fcSupra)'];
            break;
        end
    end
end
if isempty(wFCex_src) && exist('Tcog','var') && istable(Tcog)
    for ci = 1:numel(wFCCands)
        if ismember(wFCCands{ci}, Tcog.Properties.VariableNames)
            wFCex(tfCog) = double(Tcog.(wFCCands{ci})(locCog(tfCog)));
            wFCex_src = [wFCCands{ci} ' (Tcog)'];
            break;
        end
    end
end
if isempty(wFCex_src)
    fprintf('  NOTE: no whole_FC_mean_excl_apathy column found in fcSupra or Tcog; running Full Z without it.\n');
    Z_full_extra = [baiFC];
    wfcStr = '';
else
    Z_full_extra = [baiFC, wFCex];
    wfcStr = sprintf(' + %s', wFCex_src);
end

Z_min_FC  = [fcMean.mean_FD, fcMean.TIV];
Z_full_FC = [fcMean.mean_FD, fcMean.TIV, siteDumFC, Z_full_extra];
fprintf('  Min  Z: mean_FD + TIV.\n');
fprintf('  Full Z: mean_FD + TIV + site_dummies(k=%d) + brain_atrophy_index%s.\n', ...
    size(siteDumFC,2), wfcStr);

sV_d = slopeAll(:,1);
rhoMd = nan(nAllFC,1); pMd = nan(nAllFC,1); nMd = zeros(nAllFC,1);
rhoFd = nan(nAllFC,1); pFd = nan(nAllFC,1); nFd = zeros(nAllFC,1);
for fi = 1:nAllFC
    fcV = fcData(:,fi);
    okM = ~isnan(fcV) & ~isnan(sV_d) & all(~isnan(Z_min_FC ),2); nMd(fi)=sum(okM);
    if nMd(fi)>=10, [rhoMd(fi), pMd(fi)] = partialcorr(fcV(okM), sV_d(okM), Z_min_FC(okM,:),  'Type','Spearman'); end
    okF = ~isnan(fcV) & ~isnan(sV_d) & all(~isnan(Z_full_FC),2); nFd(fi)=sum(okF);
    if nFd(fi)>=10, [rhoFd(fi), pFd(fi)] = partialcorr(fcV(okF), sV_d(okF), Z_full_FC(okF,:),'Type','Spearman'); end
end
vMd = ~isnan(pMd); pFDR_Md = nan(nAllFC,1); if any(vMd), pFDR_Md(vMd) = fdr_bh(pMd(vMd)); end
vFd = ~isnan(pFd); pFDR_Fd = nan(nAllFC,1); if any(vFd), pFDR_Fd(vFd) = fdr_bh(pFd(vFd)); end

fprintf('\n  FC components (FDR within each Z):\n');
fprintf('  %-40s   %7s %7s %8s    %7s %7s %8s    %s\n', ...
    'Component','rho_min','p_min','pFDR_min','rho_full','p_full','pFDR_full','n(min/full)');
fcRows = {};
for fi = 1:nAllFC
    if isnan(rhoMd(fi)), s1='     —'; s1p='     —'; s2='     —'; sm=' ';
    else, s1=sprintf('%+7.3f',rhoMd(fi)); s1p=sprintf('%7.4f',pMd(fi)); s2=sprintf('%7.4f',pFDR_Md(fi)); sm=pToStars(pFDR_Md(fi)); end
    if isnan(rhoFd(fi)), s3='     —'; s3p='     —'; s4='     —'; sf=' ';
    else, s3=sprintf('%+7.3f',rhoFd(fi)); s3p=sprintf('%7.4f',pFd(fi)); s4=sprintf('%7.4f',pFDR_Fd(fi)); sf=pToStars(pFDR_Fd(fi)); end
    fprintf('  %-40s   %s %s %s%s   %s %s %s%s   %d/%d\n', ...
        fcLabels{fi}, s1, s1p, s2, sm, s3, s3p, s4, sf, nMd(fi), nFd(fi));
    fcRows(end+1,:) = {fcLabels{fi}, rhoMd(fi), pMd(fi), pFDR_Md(fi), nMd(fi), ...
                                      rhoFd(fi), pFd(fi), pFDR_Fd(fi), nFd(fi)}; %#ok
end
fcOutCsv = sprintf('Table_Fig4_FC_x_CogSlope_upwardMatched_t%s.csv', acmThreshold);
writetable(cell2table(fcRows,'VariableNames', ...
    {'Component','rho_min','p_min','pFDR_min','n_min','rho_full','p_full','pFDR_full','n_full'}), ...
    fullfile(outputDir, fcOutCsv));
fprintf('\n  -> %s saved.\n', fcOutCsv);

%% ========================================================================
%  5b. DAT ANALYSES — Dopaminergic integrity × FC × Cognition
% =========================================================================
if hasDAT
    fprintf('\n[5b] DAT analyses...\n');

    fcMean.DAT_mean = nanmean([fcMean.DAT_caudate, fcMean.DAT_putamen], 2);
    fprintf('    DAT composite (mean caudate+putamen): n=%d\n', sum(~isnan(fcMean.DAT_mean)));

    % --- (i) Sensitivity: FC → Global Slope | NBS-matched + DAT ----------
    %     Primary path b used NBS-matched covariates; add DAT on top.
    fprintf('\n    (i) Sensitivity: FC -> Global Slope | primary (%s) + DAT\n', strjoin(primaryCovLbl,', '));
    fprintf('    %-35s  %8s  %8s  %8s  %8s\n', ...
        '', 'r_p(prim)', 'p(prim)', 'r_p+DAT', 'p(+DAT)');
    fprintf('    %s\n', repmat('-',1,78));

    % Pre-allocate all DAT structs to length nAllFC (indexed by fi).
    emptySens = struct('rP0',NaN,'pP0',NaN,'rPd',NaN,'pPd',NaN,'pPd_FDR',NaN);
    emptyDatFC = struct('rP',NaN,'pP',NaN,'pFDR',NaN,'n',NaN, ...
                        'dat_vals',[],'fc_vals',[]);
    datSens = repmat(emptySens, 1, nAllFC);
    datFC   = repmat(emptyDatFC, 1, nAllFC);

    for ci = 1:nSupraFC
        fi = supraIdx(ci);
        fcV = fcData(:, fi);
        sV  = slopeAll(:,1);

        rP0 = compR(fi).rP;  pP0 = compR(fi).pP;   % primary Spearman | FD+TIV+baselineG
        % Cog slope outcome: include baselineG (matches §4 primary), then add DAT.
        covsPlusDAT = [primaryCovs, baselineG, fcMean.DAT_mean];
        [rPd, pPd] = partialcorr(fcV, sV, covsPlusDAT, 'Type', primaryMethod, 'Rows', 'complete');
        nD = sum(~isnan(fcV) & ~isnan(sV) & all(~isnan(covsPlusDAT),2));
        compR(fi).rP_DAT = rPd; compR(fi).pP_DAT = pPd;
        datSens(fi).rP0 = rP0; datSens(fi).pP0 = pP0;
        datSens(fi).rPd = rPd; datSens(fi).pPd = pPd;

        fprintf('    %-35s  %+8.3f  %.4f%s  %+8.3f  %.4f%s  n=%d\n', ...
            fcLabels{fi}, rP0, pP0, pToStars(pP0), rPd, pPd, pToStars(pPd), nD);
    end

    % FDR across supraIdx components for DAT-adjusted results
    sensP = [datSens(supraIdx).pPd];
    fdrDAT = fdr_bh(sensP);
    for k = 1:nSupraFC
        datSens(supraIdx(k)).pPd_FDR = fdrDAT(k);
    end

    % --- (ii) DAT × FC interaction on global slope -------------------------
    fprintf('\n    (ii) DAT x FC interaction on global cog slope\n');
    for ci = 1:nSupraFC
        fi = supraIdx(ci);
        fcV = fcData(:, fi); sV=slopeAll(:,1);
        okD = ~isnan(fcV)&~isnan(sV)&~isnan(fcMean.DAT_mean)&~isnan(fcMean.mean_FD)&~isnan(fcMean.TIV);
        if sum(okD)<20, continue; end
        fc_d=fcV(okD); dat_d=fcMean.DAT_mean(okD); sl_d=sV(okD);
        cov_d=[fcMean.mean_FD(okD), fcMean.TIV(okD)];
        mdl_int = fitlm([fc_d,dat_d,fc_d.*dat_d,cov_d], sl_d, ...
            'VarNames',{'FC','DAT','FC_x_DAT','FD','TIV','Slope'});
        intRow = strcmp(mdl_int.CoefficientNames,'FC_x_DAT');
        fprintf('    %-35s  FC_x_DAT: b=%+.4f  p=%.4f %s\n', ...
            fcLabels{fi}, mdl_int.Coefficients.Estimate(intRow), ...
            mdl_int.Coefficients.pValue(intRow), pToStars(mdl_int.Coefficients.pValue(intRow)));
    end

    % --- (iii) DAT → FC: does dopamine drive FC disruption? ----------------
    %     Primary: Spearman partial r | age + sex + FD + eTIV.
    %     Covariate set matches the cross-sectional FreeSurfer ROI convention
    %     (Figure4_Structural.m M1 Primary), so DAT-FC and DAT-volume effects
    %     are on the same covariate footing for structural/functional cross-
    %     reference. fcMean.TIV ≡ eTIV (same FreeSurfer-derived value).
    %     FDR family: across all components in supraIdx.
    covsDATfc    = [fcMean.mean_FD, fcMean.TIV, fcMean.age, sexNum];
    covsDATfcLbl = {'FD','eTIV','age','sex'};

    fprintf('\n    (iii) DAT -> FC components (supraIdx) | primary (%s)\n', ...
        strjoin(covsDATfcLbl,', '));
    fprintf('    %-35s  %8s  %8s  %8s\n', '', 'rho', 'p_unc', 'p_FDR');
    fprintf('    %s\n', repmat('-',1,65));

    for ci = 1:nSupraFC
        fi = supraIdx(ci);
        fcV = fcData(:, fi);
        % Primary: Spearman partial r | age + sex + FD + eTIV (FS-volume parity)
        [rPdat, pPdat] = partialcorr(fcMean.DAT_mean, fcV, covsDATfc, ...
            'Type', primaryMethod, 'Rows', 'complete');
        nD = sum(~isnan(fcMean.DAT_mean) & ~isnan(fcV) & all(~isnan(covsDATfc),2));
        datFC(fi).rP = rPdat; datFC(fi).pP = pPdat; datFC(fi).n = nD;
        % Sensitivity: Pearson partial on the same covariate set
        okB = ~isnan(fcV)&~isnan(fcMean.DAT_mean) & all(~isnan(covsDATfc),2);
        datFC(fi).rP_basic = partialcorr(fcMean.DAT_mean(okB), fcV(okB), ...
            covsDATfc(okB,:), 'Rows','complete');
        dfB = sum(okB) - 2 - size(covsDATfc,2);
        tB  = datFC(fi).rP_basic*sqrt(dfB/(1-datFC(fi).rP_basic^2));
        datFC(fi).pP_basic = 2*(1-tcdf(abs(tB),dfB));
        % Store raw vectors for panel D scatter (residualization done inline)
        okD = ~isnan(fcV)&~isnan(fcMean.DAT_mean) & all(~isnan(covsDATfc),2);
        datFC(fi).dat_vals = fcMean.DAT_mean(okD);
        datFC(fi).fc_vals  = fcV(okD);
    end

    % FDR across supraIdx components for DAT→FC
    fcP = [datFC(supraIdx).pP];
    fdrDATfc = fdr_bh(fcP);
    for k = 1:nSupraFC
        datFC(supraIdx(k)).pFDR = fdrDATfc(k);
    end

    for ci = 1:nSupraFC
        fi = supraIdx(ci);
        fprintf('    %-35s  %+8.3f  %.4f  %.4f %s  n=%d   [basic: r=%+5.3f, p=%.4f]\n', ...
            fcLabels{fi}, datFC(fi).rP, datFC(fi).pP, datFC(fi).pFDR, ...
            pToStars(datFC(fi).pFDR), datFC(fi).n, ...
            datFC(fi).rP_basic, datFC(fi).pP_basic);
    end

    % Save DAT results (all components in supraIdx)
    datRows = {};
    for ci = 1:nSupraFC
        fi = supraIdx(ci);
        datRows{end+1} = {fcLabels{fi}, ...
            compR(fi).rP, compR(fi).pP, compR(fi).rP_DAT, compR(fi).pP_DAT, datSens(fi).pPd_FDR, ...
            datFC(fi).rP, datFC(fi).pP, datFC(fi).pFDR}; %#ok
    end
    writetable(cell2table(vertcat(datRows{:}), ...
        'VariableNames', {'FC_Component','r_p','p_unc','r_p_DAT','p_DAT','p_DAT_FDR', ...
                          'DAT_FC_r','DAT_FC_p','DAT_FC_pFDR'}), ...
        fullfile(outputDir, 'Table_Fig4_DAT_Results.csv'));
else
    fprintf('\n[5b] DAT analyses skipped (no DAT data).\n');
    emptySens = struct('rP0',NaN,'pP0',NaN,'rPd',NaN,'pPd',NaN,'pPd_FDR',NaN);
    emptyDatFC = struct('rP',NaN,'pP',NaN,'pFDR',NaN,'n',NaN, ...
                        'dat_vals',[],'fc_vals',[]);
    datSens = repmat(emptySens, 1, nAllFC);
    datFC   = repmat(emptyDatFC, 1, nAllFC);
end

%% ========================================================================
%  6.  MEDIATION — Apathy → FC → Global Cog Decline
%      Run for every primary component; idxMediation flags the one used
%      for Panel B visualization. Covariates = FD + TIV only.
% =========================================================================
fprintf('\n[6] Mediation: Apathy -> FC -> Global Decline (bootstrap %d reps)...\n', nBoot);

medResults = repmat(struct('fi',NaN,'label','','n',NaN,'a',NaN,'a_p',NaN, ...
    'b',NaN,'b_p',NaN,'cp',NaN,'cp_p',NaN,'c',NaN,'c_p',NaN, ...
    'ab',NaN,'ab_ci',[NaN NaN],'ab_sig',false,'prop',NaN, ...
    'ab_boot',[], ...
    'ab_ci_bca',[NaN NaN],'ab_sig_bca',false,'bca_z0',NaN,'bca_ahat',NaN), ...
    1, nAllFC);

medIdxList = unique([primaryIdx(:); idxMediation]);  % always include idxMediation

for fi = medIdxList(:)'
    medFC = fcData(:, fi);
    medSlope = slopeAll(:,1);   % global
    medAp = fcMean.apathy_composite;
    % Covariates = FD + TIV only. No per-outcome baselineG: the BLUP slope
    % already absorbs level-side covariates (age, sex, educ, etc.), so
    % adding baseline_g here would over-control. See feedback memory.
    okM = ~isnan(medFC)&~isnan(medSlope)&~isnan(medAp)& ...
          ~isnan(fcMean.mean_FD)&~isnan(fcMean.TIV);
    nM = sum(okM);
    fc_m=medFC(okM); sl_m=medSlope(okM); ap_m=medAp(okM);
    cov_m=[fcMean.mean_FD(okM), fcMean.TIV(okM)];

    % Paths
    mdl_a = fitlm([ap_m,cov_m], fc_m);
    a_ = mdl_a.Coefficients.Estimate(2); a_p = mdl_a.Coefficients.pValue(2);
    mdl_b = fitlm([fc_m,ap_m,cov_m], sl_m);
    b_ = mdl_b.Coefficients.Estimate(2); b_p = mdl_b.Coefficients.pValue(2);
    cp_ = mdl_b.Coefficients.Estimate(3); cp_p = mdl_b.Coefficients.pValue(3);
    mdl_c = fitlm([ap_m,cov_m], sl_m);
    c_ = mdl_c.Coefficients.Estimate(2); c_p = mdl_c.Coefficients.pValue(2);
    ab_ = a_ * b_;

    % Bootstrap
    ab_boot = zeros(nBoot,1);
    for bi = 1:nBoot
        ix = randi(nM,nM,1);
        ba = [ones(nM,1),ap_m(ix),cov_m(ix,:)]\fc_m(ix);
        bb = [ones(nM,1),fc_m(ix),ap_m(ix),cov_m(ix,:)]\sl_m(ix);
        ab_boot(bi) = ba(2)*bb(2);
    end
    ab_ci  = prctile(ab_boot,[2.5 97.5]);
    ab_sig = (ab_ci(1) > 0) || (ab_ci(2) < 0);
    prop_  = ternary(abs(c_)>1e-6, ab_/c_*100, NaN);

    medResults(fi).fi = fi;
    medResults(fi).label = fcLabels{fi};
    medResults(fi).n = nM;
    medResults(fi).a = a_; medResults(fi).a_p = a_p;
    medResults(fi).b = b_; medResults(fi).b_p = b_p;
    medResults(fi).cp = cp_; medResults(fi).cp_p = cp_p;
    medResults(fi).c = c_; medResults(fi).c_p = c_p;
    medResults(fi).ab = ab_; medResults(fi).ab_ci = ab_ci;
    medResults(fi).ab_sig = ab_sig; medResults(fi).prop = prop_;
    medResults(fi).ab_boot = ab_boot;
end

% --- Report all mediations, flag primary ---
fprintf('  %-45s  %8s  %8s  %8s  %8s  %15s  %s\n', ...
    '', 'a', 'b', 'c''', 'c', 'ab [CI]', 'mediated');
fprintf('  %s\n', repmat('-',1,110));
for fi = medIdxList(:)'
    m = medResults(fi);
    marker = ternary(fi == idxMediation, ' <- PRIMARY', '');
    fprintf('  %-45s  %+8.4f  %+8.4f  %+8.4f  %+8.4f  %+.3f[%+.3f,%+.3f] %s  %.1f%%%s\n', ...
        m.label, m.a, m.b, m.cp, m.c, ...
        m.ab, m.ab_ci(1), m.ab_ci(2), ...
        ternary(m.ab_sig,'SIG','n.s.'), m.prop, marker);
end

% Convenience: expose primary mediation as scalars (used by Panel B)
mp = medResults(idxMediation);
a=mp.a; a_p=mp.a_p; b=mp.b; b_p=mp.b_p;
cp=mp.cp; cp_p=mp.cp_p; c=mp.c; c_p=mp.c_p;
ab=mp.ab; ab_ci=mp.ab_ci; ab_sig=mp.ab_sig; prop=mp.prop;
ab_boot=mp.ab_boot; nM=mp.n;
labelMed = fcShort{idxMediation};

% Save full table
medRows = {};
for fi = medIdxList(:)'
    m = medResults(fi);
    medRows{end+1} = {m.label, m.n, m.a, m.a_p, m.b, m.b_p, m.cp, m.cp_p, ...
        m.c, m.c_p, m.ab, m.ab_ci(1), m.ab_ci(2), m.prop, ...
        double(m.ab_sig), double(fi == idxMediation)}; %#ok
end
writetable(cell2table(vertcat(medRows{:}), ...
    'VariableNames',{'FC_Component','N','a','a_p','b','b_p','cp','cp_p', ...
    'c','c_p','ab','CI_lo','CI_hi','PctMediated','IndirectSig','Primary'}), ...
    fullfile(outputDir,'Table_Fig4_Mediation.csv'));

%% ========================================================================
%  6_BCa. BCa VERIFICATION OF §6 INDIRECT-EFFECT CIs
%      Bias-corrected & accelerated bootstrap on the same ab_boot draws
%      and same data as §6, with jackknife-derived acceleration. Standard
%      refinement for indirect effects, whose sampling distribution is
%      asymmetric (Efron 1987; MacKinnon, Lockwood & Williams 2004). The
%      percentile CIs in §6 are kept as the primary report; BCa is a
%      robustness check against bias and skew. (Section §6c later in the
%      file is a different sensitivity — the Panel-d subset cohort.)
% =========================================================================
fprintf('\n[6_BCa] BCa verification of §6 percentile CIs (jackknife acceleration)...\n');

for fi = medIdxList(:)'
    % Re-derive the §6 working set (same filter, same covariates)
    medFC = fcData(:, fi);
    medSlope = slopeAll(:,1);
    medAp = fcMean.apathy_composite;
    okM = ~isnan(medFC)&~isnan(medSlope)&~isnan(medAp)& ...
          ~isnan(fcMean.mean_FD)&~isnan(fcMean.TIV);
    nM = sum(okM);
    fc_m=medFC(okM); sl_m=medSlope(okM); ap_m=medAp(okM);
    cov_m=[fcMean.mean_FD(okM), fcMean.TIV(okM)];

    % Jackknife: leave-one-out indirect effect, same paths as §6.
    ab_jack = zeros(nM,1);
    for ji = 1:nM
        keep = true(nM,1); keep(ji) = false;
        ba = [ones(nM-1,1), ap_m(keep), cov_m(keep,:)] \ fc_m(keep);
        bb = [ones(nM-1,1), fc_m(keep), ap_m(keep), cov_m(keep,:)] \ sl_m(keep);
        ab_jack(ji) = ba(2)*bb(2);
    end

    [ci_bca, z0, ahat] = bcaInterval( ...
        medResults(fi).ab, medResults(fi).ab_boot, ab_jack, 0.05);
    medResults(fi).ab_ci_bca  = ci_bca;
    medResults(fi).ab_sig_bca = (ci_bca(1) > 0) || (ci_bca(2) < 0);
    medResults(fi).bca_z0     = z0;
    medResults(fi).bca_ahat   = ahat;
end

% Side-by-side report (Percentile vs BCa)
fprintf('  %-45s  %9s  %25s  %25s\n', '', 'ab', 'Percentile 95% CI', 'BCa 95% CI');
fprintf('  %s\n', repmat('-',1,112));
for fi = medIdxList(:)'
    m = medResults(fi);
    marker = ternary(fi == idxMediation, ' <- PRIMARY', '');
    fprintf('  %-45s  %+9.4f  [%+.4f, %+.4f] %-3s  [%+.4f, %+.4f] %-3s%s\n', ...
        m.label, m.ab, ...
        m.ab_ci(1),     m.ab_ci(2),     ternary(m.ab_sig,    'SIG','n.s.'), ...
        m.ab_ci_bca(1), m.ab_ci_bca(2), ternary(m.ab_sig_bca,'SIG','n.s.'), ...
        marker);
end

% Expose BCa scalars for the primary mediator (parallel to §6's `ab_ci`)
ab_ci_bca  = medResults(idxMediation).ab_ci_bca;
ab_sig_bca = medResults(idxMediation).ab_sig_bca;

% Save BCa table (parallel to Table_Fig4_Mediation.csv, with both CIs side-by-side)
bcaRows = {};
for fi = medIdxList(:)'
    m = medResults(fi);
    bcaRows{end+1} = {m.label, m.n, m.ab, ...
        m.ab_ci(1),     m.ab_ci(2),     double(m.ab_sig), ...
        m.ab_ci_bca(1), m.ab_ci_bca(2), double(m.ab_sig_bca), ...
        m.bca_z0, m.bca_ahat, double(fi == idxMediation)}; %#ok
end
writetable(cell2table(vertcat(bcaRows{:}), ...
    'VariableNames',{'FC_Component','N','ab', ...
                     'Pct_CI_lo','Pct_CI_hi','Pct_Sig', ...
                     'BCa_CI_lo','BCa_CI_hi','BCa_Sig', ...
                     'BCa_z0','BCa_ahat','Primary'}), ...
    fullfile(outputDir,'Table_Fig4_Mediation_BCa.csv'));
fprintf('  -> Table_Fig4_Mediation_BCa.csv saved.\n');

%% ========================================================================
%  6b. AGE-STRATIFIED MEDIATION (needed by Figure 4 Panel B)
%      Compute bootstrap mediation within younger and older strata.
%      Covariates: FD + TIV (matches the primary mediation).
% =========================================================================
fprintf('\n[6b] Age-stratified mediation (median split)...\n');
ageMed_fig = nanmedian(fcMean.age);
strata_fig = {'Younger', fcMean.age <= ageMed_fig; ...
              'Older',   fcMean.age >  ageMed_fig};
ageStratBoot_fig = cell(2,1);
ageStratAB_fig   = [NaN NaN];
ageStratCI_fig   = [NaN NaN; NaN NaN];
ageStratN_fig    = [0 0];
for si = 1:2
    sel = strata_fig{si,2};
    mFC = fcData(:, idxMediation);
    mSl = slopeAll(:, 1);
    mAp = fcMean.apathy_composite;
    % Primary: FD + TIV only (no baselineG per user spec, see feedback memory).
    mCov = [fcMean.mean_FD, fcMean.TIV];
    ok = sel & ~isnan(mFC)&~isnan(mSl)&~isnan(mAp)&all(~isnan(mCov),2);
    n_ = sum(ok);
    fc_=mFC(ok); sl_=mSl(ok); ap_=mAp(ok); cov_=mCov(ok,:);
    ba = [ones(n_,1), ap_, cov_] \ fc_;  a_ = ba(2);
    bb = [ones(n_,1), fc_, ap_, cov_] \ sl_;  b_ = bb(2);
    ab_ = a_*b_;
    ab_b = zeros(nBoot,1);
    for bi = 1:nBoot
        ix = randi(n_,n_,1);
        baB = [ones(n_,1), ap_(ix), cov_(ix,:)] \ fc_(ix);
        bbB = [ones(n_,1), fc_(ix), ap_(ix), cov_(ix,:)] \ sl_(ix);
        ab_b(bi) = baB(2)*bbB(2);
    end
    ageStratBoot_fig{si} = ab_b;
    ageStratAB_fig(si)   = ab_;
    ageStratCI_fig(si,:) = prctile(ab_b, [2.5 97.5]);
    ageStratN_fig(si)    = n_;
    fprintf('    %s (n=%d): ab=%+.4f [%+.4f, %+.4f]\n', ...
        strata_fig{si,1}, n_, ab_, ageStratCI_fig(si,1), ageStratCI_fig(si,2));
end
deltaAB_fig = ageStratBoot_fig{2} - ageStratBoot_fig{1};
deltaCI_fig = prctile(deltaAB_fig, [2.5 97.5]);
deltaSig_fig = (deltaCI_fig(1) > 0) || (deltaCI_fig(2) < 0);
fprintf('    Δab (older - younger): %+.4f [%+.4f, %+.4f] %s\n', ...
    mean(deltaAB_fig), deltaCI_fig(1), deltaCI_fig(2), ...
    ternary(deltaSig_fig,'SIG','n.s.'));

%% ========================================================================
%  6b-DAT. AGE-STRATIFIED MEDIATION — DAT-controlled SENSITIVITY
%      Tests whether the age-moderated mediation is residual to baseline
%      dopaminergic state (DAT_mean partialled in both path a and path b).
%      Motivation: §5b (i) showed full-sample Path b for Ap(-)1 collapses
%      from +0.135 to +0.073 (n.s.) when DAT_mean is added — i.e., the
%      apathy-FC-cog signal is not fully residual to dopaminergic status
%      at baseline. This block does the parallel test stratified by age.
%      Output: per-stratum ab + CI under DAT control, plus Δab + CI; and
%      a side-by-side comparison with the §6b primary (FD+TIV only).
% =========================================================================
fprintf('\n[6b-DAT] Age-stratified mediation with DAT_mean controlled (sensitivity)...\n');
ageStratBoot_DAT = cell(2,1);
ageStratAB_DAT   = [NaN NaN];
ageStratCI_DAT   = [NaN NaN; NaN NaN];
ageStratN_DAT    = [0 0];
for si = 1:2
    sel = strata_fig{si,2};
    mFC = fcData(:, idxMediation);
    mSl = slopeAll(:, 1);
    mAp = fcMean.apathy_composite;
    mCov_DAT = [fcMean.mean_FD, fcMean.TIV, fcMean.DAT_mean];
    ok = sel & ~isnan(mFC) & ~isnan(mSl) & ~isnan(mAp) & all(~isnan(mCov_DAT),2);
    n_ = sum(ok);
    fc_=mFC(ok); sl_=mSl(ok); ap_=mAp(ok); cov_=mCov_DAT(ok,:);
    ba = [ones(n_,1), ap_, cov_] \ fc_;          a_ = ba(2);
    bb = [ones(n_,1), fc_, ap_, cov_] \ sl_;     b_ = bb(2);
    ab_ = a_ * b_;
    ab_b = zeros(nBoot,1);
    for bi = 1:nBoot
        ix = randi(n_, n_, 1);
        baB = [ones(n_,1), ap_(ix), cov_(ix,:)] \ fc_(ix);
        bbB = [ones(n_,1), fc_(ix), ap_(ix), cov_(ix,:)] \ sl_(ix);
        ab_b(bi) = baB(2) * bbB(2);
    end
    ageStratBoot_DAT{si} = ab_b;
    ageStratAB_DAT(si)   = ab_;
    ageStratCI_DAT(si,:) = prctile(ab_b, [2.5 97.5]);
    ageStratN_DAT(si)    = n_;
    sig_str = ternary((ageStratCI_DAT(si,1) > 0) || (ageStratCI_DAT(si,2) < 0), 'SIG', 'n.s.');
    fprintf('    %s (n=%d, +DAT): ab=%+.4f [%+.4f, %+.4f]  %s\n', ...
        strata_fig{si,1}, n_, ab_, ageStratCI_DAT(si,1), ageStratCI_DAT(si,2), sig_str);
end
deltaAB_DAT  = ageStratBoot_DAT{2} - ageStratBoot_DAT{1};
deltaCI_DAT  = prctile(deltaAB_DAT, [2.5 97.5]);
deltaSig_DAT = (deltaCI_DAT(1) > 0) || (deltaCI_DAT(2) < 0);
fprintf('    Δab (older - younger, +DAT): %+.4f [%+.4f, %+.4f]  %s\n', ...
    mean(deltaAB_DAT), deltaCI_DAT(1), deltaCI_DAT(2), ...
    ternary(deltaSig_DAT, 'SIG', 'n.s.'));

% --- Side-by-side comparison: §6b PRIMARY (FD+TIV) vs SENSITIVITY (FD+TIV+DAT) ---
fprintf('\n  Comparison: §6b PRIMARY (FD+TIV)  vs  SENSITIVITY (FD+TIV+DAT_mean)\n');
fprintf('  %-14s  %-30s  %-30s\n', 'Stratum', 'Primary ab [95%% CI]', '+DAT ab [95%% CI]');
fprintf('  %s\n', repmat('-',1,80));
for si = 1:2
    fprintf('  %-14s  %+.4f [%+.4f, %+.4f]  %+.4f [%+.4f, %+.4f]\n', ...
        strata_fig{si,1}, ...
        ageStratAB_fig(si), ageStratCI_fig(si,1), ageStratCI_fig(si,2), ...
        ageStratAB_DAT(si), ageStratCI_DAT(si,1), ageStratCI_DAT(si,2));
end
fprintf('  %-14s  %+.4f [%+.4f, %+.4f]  %+.4f [%+.4f, %+.4f]\n', ...
    'Δ (O − Y)', ...
    mean(deltaAB_fig), deltaCI_fig(1), deltaCI_fig(2), ...
    mean(deltaAB_DAT), deltaCI_DAT(1), deltaCI_DAT(2));
fprintf('  Interpretation:\n');
fprintf('    older ab survives + Δ survives → mediation independent of DAT\n');
fprintf('    older ab n.s.    + Δ survives → moderation present but absolute mediation collapses\n');
fprintf('    older ab n.s.    + Δ n.s.     → headline largely DAT-driven\n');
fprintf('    older ab survives + Δ n.s.    → mediation present but moderation not robust\n');

%% ========================================================================
%  6c. PANEL-D SUBSET SENSITIVITY — does the FC × cog-slope relationship
%      hold in the subjects who have longitudinal hippocampus data?
%      Panel (d) is built on the n≈71 subjects with both BL and 24M MRI;
%      this asks whether the cognitive Path-b finding (full-sample n=310)
%      survives within that smaller imaging-followed cohort. If yes, the
%      structural and cognitive findings converge on the same patient
%      subset (strengthens the joint narrative). If no, panel (d)'s
%      structural cohort is a different population from the cognitive one.
% =========================================================================
fprintf('\n[6c] Panel-d subset sensitivity: Path b in n=71 longitudinal cohort...\n');
hippLongCSV_6c = fullfile(outputDir, 'Table_Hippocampus_Longitudinal_Data.csv');
if exist(hippLongCSV_6c, 'file')
    Th_6c = readtable(hippLongCSV_6c);
    [~, longIdx_6c] = ismember(fcPATNO, double(Th_6c.PATNO));
    panelD_mask = false(nFC, 1);
    valid_6c = longIdx_6c > 0;
    panelD_mask(valid_6c) = ~isnan(Th_6c.slope_left(longIdx_6c(valid_6c)));
    nPanelD = sum(panelD_mask);
    fprintf('    Panel-d subset (subjects with L hipp slope): n = %d / %d\n', nPanelD, nFC);

    % Path b: FC × global cog slope | FD + TIV. Two predictors:
    %   idxMediation = Ap(-)1 RH SM-DA  (the §6 cognitive mediator)
    %   idxStructural = Diff(-)1 Ap-dom SM-DA  (the panel-d FC predictor)
    fc_to_test = [idxMediation, idxStructural];

    fprintf('\n  %-26s  %-30s  %-30s\n', 'FC predictor', ...
        'Full sample ρ [p, n]', 'Panel-d subset ρ [p, n]');
    fprintf('  %s\n', repmat('-', 1, 90));

    sub6c_rows = {};
    for fi_ = fc_to_test
        fcV_6c = fcData(:, fi_);
        sV_6c  = slopeAll(:, 1);
        Z_6c   = [fcMean.mean_FD, fcMean.TIV];

        % Full sample
        ok_full = ~isnan(fcV_6c) & ~isnan(sV_6c) & all(~isnan(Z_6c),2);
        [rho_full, p_full] = partialcorr(fcV_6c(ok_full), sV_6c(ok_full), ...
                                         Z_6c(ok_full,:), 'Type','Spearman');
        n_full = sum(ok_full);

        % Panel-d subset
        ok_sub = ok_full & panelD_mask;
        if sum(ok_sub) >= 10
            [rho_sub, p_sub] = partialcorr(fcV_6c(ok_sub), sV_6c(ok_sub), ...
                                           Z_6c(ok_sub,:), 'Type','Spearman');
            n_sub = sum(ok_sub);
        else
            rho_sub = NaN; p_sub = NaN; n_sub = sum(ok_sub);
        end

        fprintf('  %-26s  %+.3f [p=%.4f, n=%d]    %+.3f [p=%.4f, n=%d]\n', ...
            fcShort{fi_}, rho_full, p_full, n_full, rho_sub, p_sub, n_sub);
        sub6c_rows{end+1} = {fcShort{fi_}, rho_full, p_full, n_full, ...
                             rho_sub, p_sub, n_sub}; %#ok
    end

    if ~isempty(sub6c_rows)
        writetable(cell2table(vertcat(sub6c_rows{:}), ...
            'VariableNames', {'FC_Component', ...
                              'rho_full', 'p_full', 'n_full', ...
                              'rho_panelD', 'p_panelD', 'n_panelD'}), ...
            fullfile(outputDir, 'Table_Fig4_PanelD_Subset_Sensitivity.csv'));
        fprintf('\n  -> Table_Fig4_PanelD_Subset_Sensitivity.csv saved.\n');
    end

    fprintf('\n  Interpretation:\n');
    fprintf('    rho_panelD comparable to rho_full → cog and structural cohorts converge\n');
    fprintf('    rho_panelD attenuated/n.s.       → panel-d subset is a distinct population\n');
    fprintf('                                        for the cognitive-FC relationship\n');
else
    fprintf('    %s not found; run Figure4_Structural.m §13 first to enable this sensitivity.\n', ...
        hippLongCSV_6c);
end

%% ========================================================================
%  7.  PUBLICATION FIGURE  (2×2 layout, 30 × 15 cm, tiledlayout(2,12))
%      Row 1 (2 × [1,6] = 12 wide):
%        (a) Simple-slopes density scatter — Ap(-)1 × cog, age-stratified
%        (b) Age-stratified mediation bootstrap (KDE + CI bars)
%      Row 2 (2 × [1,6] = 12 wide):
%        (c) Structural forest (top 15 ROIs, from CSV, lollipop)
%        (d) Hippocampus slope × FC density scatter (Diff(-)1, longitudinal)
%      Density-scatter cmap = magma-style warm sequential.
%      Removed in 2026-05-07 layout passes: bar-plot panel (was old "a"),
%      DAT-slope scatter (was old "e"), DAT cross-sectional (relegated to
%      supplementary — analysis still runs in §5b iii but isn't plotted in
%      the main figure). The subdomain heatmap lives only in supplementary §11.
% =========================================================================
fprintf('\n[7] Building publication figure (4 panels in 2x2 layout)...\n');

fig4 = figure('Units','centimeters','Position',[1 1 15 15],'Color',COL.bg);
set(fig4, ...
    'DefaultAxesFontName','Helvetica','DefaultTextFontName','Helvetica', ...
    'DefaultAxesFontWeight','normal', ...        % axis tick NUMBERS not bold (matches Figure 2)
    'DefaultTextFontWeight','normal', ...        % text() annotations NOT bold
    'DefaultAxesTitleFontWeight','bold');
tl = tiledlayout(2, 12, 'TileSpacing','compact', 'Padding','compact');

colY = [0.35 0.70 0.90];    % Okabe-Ito sky blue: younger   (REVISED R1.10)
colO = [0.80 0.40 0.00];    % Okabe-Ito vermillion: older   (REVISED R1.10)

% ----- (a) Simple-slopes density scatter: FC x Cog Slope by age ----------
%   Density-coloured per stratum (alpha ∝ local KDE).
ax_b = nexttile(1, [1,6]); set(ax_b,'Color',COL.bg); hold on;

fcV_a  = fcData(:, idxMediation);
sV_a   = slopeAll(:, 1);
ageV_a = fcMean.age;
Z_a    = [fcMean.mean_FD, fcMean.TIV];
ok_a   = ~isnan(fcV_a) & ~isnan(sV_a) & ~isnan(ageV_a) & all(~isnan(Z_a),2);

grpCols_a = {colY, colO};
grpLabs_a = {sprintf('Younger (\\leq%.0fy)', ageMed_fig), ...
             sprintf('Older (>%.0fy)', ageMed_fig)};
hLeg_a = gobjects(2,1);

for gi = 1:2
    if gi == 1, sel_a = ok_a & (ageV_a <= ageMed_fig);
    else,       sel_a = ok_a & (ageV_a >  ageMed_fig); end
    Zs_a = Z_a(sel_a,:);
    xR_a = fcV_a(sel_a) - Zs_a*(Zs_a\fcV_a(sel_a));
    yR_a = sV_a(sel_a)  - Zs_a*(Zs_a\sV_a(sel_a));
    thisCol_a = grpCols_a{gi};

    % Per-stratum density → alpha modulation
    try
        dens_a = mvksdensity([xR_a yR_a], [xR_a yR_a], ...
            'Bandwidth', [std(xR_a)*0.35, std(yR_a)*0.35]);
    catch
        dens_a = ones(numel(xR_a), 1);
    end
    densN_a = (dens_a - min(dens_a)) / max(max(dens_a) - min(dens_a), eps);
    alphaVec_a = 0.12 + 0.68 * densN_a;

    scatter(ax_b, xR_a, yR_a, 30, repmat(thisCol_a, numel(xR_a), 1), ...
        'filled', 'MarkerFaceAlpha', 'flat', 'AlphaData', alphaVec_a);

    % Robust (bisquare) regression line
    % REVISED (R1.10): younger dashed, older solid, so the strata remain
    % separable with no colour information (grayscale print).
    lineStyle_a = ternary(gi == 1, '--', '-');
    try
        b_a = robustfit(xR_a, yR_a);
        xF_a = linspace(min(xR_a), max(xR_a), 100);
        hLeg_a(gi) = plot(ax_b, xF_a, b_a(1)+b_a(2)*xF_a, lineStyle_a, ...
            'Color', thisCol_a, 'LineWidth', 2.2);
    catch
        hLeg_a(gi) = plot(ax_b, NaN, NaN, lineStyle_a, 'Color', thisCol_a);
    end

    % Per-stratum Spearman partial (FD+TIV+baselineG, matches §4 primary)
    Z_strat = [primaryCovs(sel_a,:), baselineG(sel_a)];
    [rSp_a, pSp_a] = partialcorr(fcV_a(sel_a), sV_a(sel_a), Z_strat, 'Type','Spearman','Rows','complete');
    rSp_strat{gi} = rSp_a; pSp_strat{gi} = pSp_a; nSp_strat{gi} = sum(sel_a); %#ok
end

% Stats block in bottom-left: overall first, then per-stratum
% REVISED (R1.9): all stats text neutral dark; stratum identity is carried
% by the name prefix (colour-coded fonts were illegible in grayscale).
text(ax_b, 0.03, 0.19, ...
    sprintf('Overall:  \\rho = %+.2f,  p_{FDR} %s', compR(idxMediation).rP, fmtP(compR(idxMediation).pFDR)), ...
    'Units','normalized','FontSize',6.5,'FontWeight','normal','VerticalAlignment','top', ...
    'Color',COL.txt,'BackgroundColor',[1 1 1 0.75], 'Margin', 1);
grpNames_a = {'Younger','Older'};
for gi = 1:2
    text(ax_b, 0.03, 0.19 - 0.065 - (gi-1)*0.055, ...
        sprintf('%s:  \\rho = %+.2f,  p %s', ...
        grpNames_a{gi}, rSp_strat{gi}, fmtP(pSp_strat{gi})), ...
        'Units','normalized','FontSize',6,'FontWeight','normal','VerticalAlignment','top', ...
        'Color',COL.txt,'BackgroundColor',[1 1 1 0.75], 'Margin', 1);
end

yline(ax_b, 0, '--', 'Color', COL.grid, 'LineWidth', 0.5);
xlabel('Apathy SM-DA FC (resid)', 'FontSize', 8);
ylabel('Global Cog Slope (resid)', 'FontSize', 8);

% Invisible thicker handles for legend-only appearance (shorter+thicker lines)
% REVISED (R1.10): younger proxy dashed to match its plotted line style.
hLegThickA(1) = plot(ax_b, NaN, NaN, '--', 'Color', colY, 'LineWidth', 4);
hLegThickA(2) = plot(ax_b, NaN, NaN, '-', 'Color', colO, 'LineWidth', 4);
grpLabs_a = {sprintf('Younger (n=%d)', nSp_strat{1}), ...
             sprintf('Older (n=%d)', nSp_strat{2})};
lgdA = legend(hLegThickA, grpLabs_a, 'FontSize', 6.5, 'Box', 'off', 'Location', 'northwest');
lgdA.ItemTokenSize = [10, 10];
set(lgdA, 'FontWeight', 'bold');
title('\bfa\rm  Age Moderates FC \times Decline', 'FontSize', 9, 'Color', COL.txt);
formatAxis(ax_b, COL);
ylim(ax_b, [-0.4 0.2]);

% ----- (b) Age-stratified mediation bootstrap -----------------------------
ax_c = nexttile(7, [1,6]); set(ax_c,'Color',COL.bg); hold on;

abBY = ageStratBoot_fig{1};
abBO = ageStratBoot_fig{2};
[fY, xiY] = ksdensity(abBY);
[fO, xiO] = ksdensity(abBO);

% REVISED (R1.10): younger density edge dashed (matches panel a line style).
fill(ax_c, [xiY fliplr(xiY)], [fY zeros(size(fY))], colY, ...
    'FaceAlpha', 0.35, 'EdgeColor', colY, 'LineWidth', 1.5, 'LineStyle', '--');
fill(ax_c, [xiO fliplr(xiO)], [fO zeros(size(fO))], colO, ...
    'FaceAlpha', 0.35, 'EdgeColor', colO, 'LineWidth', 1.5);

% REVISED (R1.9): zero reference capped at the density envelope (a full-axis
% xline previously ran through the CI bars, bracket, and star above).
maxDens_b = max([fY fO]);
plot(ax_c, [0 0], [0 1.02*maxDens_b], '-', 'Color', [0.15 0.15 0.15], 'LineWidth', 1.5);

% CI vertical bars (2.5th / 97.5th) + point estimate
% REVISED (R1.9): bars lifted from 0.85/0.75 x max density (which drew them
% on top of the density curves) to 1.10/1.24 x max density, into clear
% space above the curves. Nothing occludes the densities any more.
for grp_b = 1:2
    if grp_b == 1, bVec = abBY; gCol = colY; yFrac = 1.10;
    else,          bVec = abBO; gCol = colO; yFrac = 1.24; end
    ciLo_b = prctile(bVec, 2.5);
    ciHi_b = prctile(bVec, 97.5);
    ptEst_b = mean(bVec);
    yPos_b  = yFrac * maxDens_b;
    plot(ax_c, [ciLo_b ciHi_b], [yPos_b yPos_b], '-', 'Color', gCol, 'LineWidth', 2.5);
    plot(ax_c, [ciLo_b ciLo_b], yPos_b + [-0.02 0.02]*maxDens_b, '-', 'Color', gCol, 'LineWidth', 1.5);
    plot(ax_c, [ciHi_b ciHi_b], yPos_b + [-0.02 0.02]*maxDens_b, '-', 'Color', gCol, 'LineWidth', 1.5);
    plot(ax_c, ptEst_b, yPos_b, 'o', 'Color', gCol, 'MarkerFaceColor', gCol, 'MarkerSize', 5);
end

% Delta stats
deltaAB_b = ageStratBoot_fig{2} - ageStratBoot_fig{1};
deltaEst_b = mean(deltaAB_b);
deltaCI_b  = prctile(deltaAB_b, [2.5 97.5]);
pBoot_b    = 2 * min(mean(deltaAB_b > 0), mean(deltaAB_b < 0));

% Per-stratum bootstrap p (two-sided) — for stats text block
pBootY = 2 * min(mean(abBY > 0), mean(abBY < 0));
pBootO = 2 * min(mean(abBO > 0), mean(abBO < 0));

% Stats block: per-stratum p (bottom-left to leave room for bracket)
% REVISED (R1.9): neutral dark text; stratum identity is in the name prefix.
text(ax_c, 0.03, 0.10, ...
    sprintf('Younger:  p %s', fmtP(pBootY)), ...
    'Units','normalized','FontSize',6.5,'FontWeight','normal','VerticalAlignment','bottom', ...
    'Color',COL.txt,'BackgroundColor',[1 1 1 0.75], 'Margin', 1);
text(ax_c, 0.03, 0.03, ...
    sprintf('Older:  p %s', fmtP(pBootO)), ...
    'Units','normalized','FontSize',6.5,'FontWeight','normal','VerticalAlignment','bottom', ...
    'Color',COL.txt,'BackgroundColor',[1 1 1 0.75], 'Margin', 1);

% Horizontal comparison bracket spanning the two point estimates (x = effect size)
% REVISED (R1.9): bracket geometry follows the lifted CI bars (1.10/1.24)
% and sits fully above them.
yBracketY = 1.10 * maxDens_b;
yBracketO = 1.24 * maxDens_b;
ptY_b = mean(abBY);
ptO_b = mean(abBO);
yTopBr = max(yBracketY, yBracketO) + 0.14 * maxDens_b;   % above highest bar
stemDrop = 0.04 * maxDens_b;

% Left stem (down to Younger bar)
plot(ax_c, [ptY_b ptY_b], [yTopBr yBracketY + stemDrop], '-', ...
    'Color', [0.25 0.25 0.25], 'LineWidth', 0.9);
% Right stem (down to Older bar)
plot(ax_c, [ptO_b ptO_b], [yTopBr yBracketO + stemDrop], '-', ...
    'Color', [0.25 0.25 0.25], 'LineWidth', 0.9);
% Horizontal bar (spans effect-size difference)
plot(ax_c, [ptY_b ptO_b], [yTopBr yTopBr], '-', ...
    'Color', [0.25 0.25 0.25], 'LineWidth', 0.9);

% Star indicator above the bar (Age-difference significance)
starStr = pToStars(pBoot_b);
if isempty(starStr), starStr = 'n.s.'; end
text(ax_c, (ptY_b + ptO_b)/2, yTopBr + 0.015 * maxDens_b, starStr, ...
    'FontSize', 11, 'FontWeight', 'bold', 'Color', [0.15 0.15 0.15], ...
    'HorizontalAlignment','center','VerticalAlignment','bottom');

xlabel('Indirect effect (Apathy \rightarrow FC \rightarrow Cog Slope)', 'FontSize', 8);
ylabel('Bootstrap density', 'FontSize', 8);
title('\bfb\rm  Age-Moderated Mediation', 'FontSize', 9, 'Color', COL.txt);
formatAxis(ax_c, COL);
% REVISED (R1.9): headroom for the lifted CI bars, bracket, and star.
ylim(ax_c, [0, 1.60 * maxDens_b]);

% ----- (c) Structural forest — RH SM-DA × volumes ------------------------
%   Reads Table_Structural_Correlates.csv (correct covariates from
%   Figure4_Structural.m).  Lollipop: FDR=blue solid, unc=open, n.s.=gray.
ax_f = nexttile(13, [1,6]); set(ax_f,'Color',COL.bg); hold on;

cSig_f  = [0.20 0.40 0.65];
cUnc_f  = [0.45 0.60 0.80];
cNS_f   = [0.78 0.78 0.78];
cBG_f   = [0.97 0.97 0.97];

structCSV_f = fullfile(outputDir, 'Table_Structural_Correlates.csv');

if exist(structCSV_f, 'file') == 2
    Tstruct_f = readtable(structCSV_f);
    smda_label_f = fcShort{idxMediation};
    mask_f = strcmp(Tstruct_f.FC_Component, smda_label_f);
    % CSV may carry an M1/M2 sensitivity column post 2026-05-06; restrict
    % the figure panel to M1 Primary. Older CSVs without this column are
    % treated as M1-only.
    if ismember('Model', Tstruct_f.Properties.VariableNames)
        mask_f = mask_f & strcmp(Tstruct_f.Model, 'M1_Primary');
    end
    Tstruct_f = Tstruct_f(mask_f, :);

    if height(Tstruct_f) > 0
        rROI_f   = Tstruct_f.rho;
        pROI_f   = Tstruct_f.p_unc;
        pFDR_f   = Tstruct_f.p_FDR;
        nameROI_f = Tstruct_f.ROI;
        nROI_f   = numel(rROI_f);

        % Select top ROIs: all FDR<.05, all p<.01, then fill to 15 by |r|
        maxShow_f = 15;
        isFDR_all = ~isnan(pFDR_f) & pFDR_f < 0.05;
        isUnc_all = ~isnan(pROI_f) & pROI_f < 0.01 & ~isFDR_all;
        keepIdx_f = find(isFDR_all | isUnc_all);
        if numel(keepIdx_f) < maxShow_f
            remaining = setdiff((1:nROI_f)', keepIdx_f);
            [~, rOrd] = sort(abs(rROI_f(remaining)), 'descend');
            nFill = min(maxShow_f - numel(keepIdx_f), numel(remaining));
            keepIdx_f = [keepIdx_f; remaining(rOrd(1:nFill))];
        end

        % Sort kept ROIs by |rho| ascending (smallest at bottom, largest at top)
        [~, sK] = sort(abs(rROI_f(keepIdx_f)), 'ascend');
        showOrd_f = keepIdx_f(sK);
        nShow_f   = numel(showOrd_f);

        % Identify subcort/cort boundary within shown ROIs
        subNames_f = {'Hippocampus','Amygdala','Caudate','Putamen','Accumbens', ...
                      'Thalamus','Pallidum','Cereb Cortex','Brain Stem'};
        isSub_shown = false(nShow_f, 1);
        for si = 1:numel(subNames_f)
            isSub_shown = isSub_shown | strcmpi(nameROI_f(showOrd_f), subNames_f{si});
        end
        nShowCort_f = sum(~isSub_shown);

        for ri = 1:nShow_f
            oi  = showOrd_f(ri);
            y_f = ri;
            rho = rROI_f(oi);
            if isnan(rho), continue; end

            isFDR_f = ~isnan(pFDR_f(oi)) && pFDR_f(oi) < 0.05;
            isUnc_f = ~isnan(pROI_f(oi)) && pROI_f(oi) < 0.01;

            if mod(ri,2) == 0
                fill(ax_f, [-0.35 0.35 0.35 -0.35], ...
                    [y_f-0.45 y_f-0.45 y_f+0.45 y_f+0.45], ...
                    cBG_f, 'EdgeColor','none','FaceAlpha',0.5);
            end

            if isFDR_f
                plot(ax_f, [0 rho], [y_f y_f], '-', 'Color', cSig_f, 'LineWidth', 1.8);
                plot(ax_f, rho, y_f, 'o', 'MarkerSize', 7, 'MarkerFaceColor', cSig_f, ...
                    'MarkerEdgeColor', 'w', 'LineWidth', 0.8);
                ha_f = ternary(rho >= 0, 'left', 'right');
                text(ax_f, rho + sign(rho)*0.028, y_f, sprintf('%.2f', rho), ...
                    'FontSize', 6.5, 'FontWeight', 'bold', 'Color', cSig_f, ...
                    'VerticalAlignment', 'middle', 'HorizontalAlignment', ha_f);
            elseif isUnc_f
                plot(ax_f, [0 rho], [y_f y_f], '-', 'Color', cUnc_f, 'LineWidth', 1.2);
                plot(ax_f, rho, y_f, 'o', 'MarkerSize', 6, 'MarkerFaceColor', 'w', ...
                    'MarkerEdgeColor', cUnc_f, 'LineWidth', 1.4);
                ha_f = ternary(rho >= 0, 'left', 'right');
                text(ax_f, rho + sign(rho)*0.028, y_f, sprintf('%.2f', rho), ...
                    'FontSize', 6, 'Color', [0.35 0.35 0.35], ...
                    'VerticalAlignment','middle', 'HorizontalAlignment', ha_f);
            else
                plot(ax_f, [0 rho], [y_f y_f], '-', 'Color', [0.85 0.85 0.85], 'LineWidth', 0.6);
                plot(ax_f, rho, y_f, '.', 'MarkerSize', 6, 'Color', cNS_f);
            end
        end

        xline(ax_f, 0, '-', 'Color', [0.4 0.4 0.4], 'LineWidth', 0.8);

        % Y-tick labels: all bold, color + italic encode significance tier
        yLabels_f = cell(nShow_f, 1);
        for ri = 1:nShow_f
            oi  = showOrd_f(ri);
            lbl = char(nameROI_f{oi});
            if ~isnan(pFDR_f(oi)) && pFDR_f(oi) < 0.05
                yLabels_f{ri} = ['\bf\color[rgb]{0.20,0.40,0.65}' lbl];
            elseif ~isnan(pROI_f(oi)) && pROI_f(oi) < 0.01
                yLabels_f{ri} = ['\bf\it\color[rgb]{0.30,0.30,0.30}' lbl];
            else
                yLabels_f{ri} = ['\bf\color[rgb]{0.55,0.55,0.55}' lbl];
            end
        end

        set(ax_f, 'YTick', 1:nShow_f, 'YTickLabel', yLabels_f, ...
            'FontSize', 8, 'TickDir', 'out', 'Box', 'off', ...
            'YLim', [0.3, nShow_f+0.7], 'XLim', [-0.35 0.35]);
        xlabel('\rho_{Sp} (partial)', 'FontSize', 8);
        text(ax_f, 0.03, 0.03, sprintf('%d / %d ROIs shown', nShow_f, nROI_f), ...
            'Units','normalized','FontSize',6.5,'FontWeight','normal','Color',COL.txt, ...
            'VerticalAlignment','bottom','BackgroundColor',[1 1 1 0.75],'Margin',1);
    else
        text(ax_f, 0.5, 0.5, 'No RH SM-DA data in CSV', ...
            'Units','normalized','FontSize',7,'Color',[0.5 0.5 0.5], ...
            'HorizontalAlignment','center','VerticalAlignment','middle');
        set(ax_f, 'XTick', [], 'YTick', []);
    end
else
    text(ax_f, 0.5, 0.5, {'Run Figure4\_Structural.m first', ...
        'to generate Table\_Structural\_Correlates.csv'}, ...
        'Units','normalized','FontSize',7,'Color',[0.5 0.5 0.5], ...
        'HorizontalAlignment','center','VerticalAlignment','middle');
    set(ax_f, 'XTick', [], 'YTick', []);
end

title('\bfc\rm  Structural Correlates', 'FontSize', 9, 'Color', COL.txt);
formatAxis(ax_f, COL);

% ----- (d) Longitudinal hippocampus × Ap-dom SM-DA FC density scatter ----
%   Loads per-subject longitudinal hippocampus slopes (computed in
%   Figure4_Structural.m §13). Primary uses raw slopes; winsorised copies
%   are available for sensitivity only.
%   FC predictor: Diff(-)1 Ap-dominant SM-DA (idxStructural).
%   Note: panel (g) uses a DIFFERENT FC component than panels d/e (which
%   use Ap(-)1 RH SM-DA = idxMediation). Diff(-)1 is the FDR-significant
%   L hipp correlate under primary Spearman partial. The Ap(-)1 Pearson
%   association is outlier-driven (Cook's D ≈ 0.49 from a single PATNO,
%   ~10× the 4/n threshold); the Spearman is n.s. (p≈0.15). Diff(-)1, by
%   contrast, has Spearman ≥ Pearson (robust monotonic effect signature).
%   Covariates (M1 PRIMARY — matches the headline write-up stat):
%     age + sex + FD + EDUCYRS + BL_hipp_vol.
%   No explicit eTIV covariate: §13 uses HC-normative-corrected slopes
%   (within-subject eTIV drift implicitly removed via the Buckner-style
%   HC β subtraction at each timepoint). Cross-sectional §1-§12 use the
%   uniform eTIV-as-covariate method, but for slope outcomes that approach
%   leaks within-subject eTIV-drift noise — so §13 keeps HC-normative.
%   Disease-severity sensitivity (+ duration + LEDD + UPDRS-3) is
%   reported in §14 and Table_Hippocampus_LM_Targeted.csv (Model column).
%   Stats annotation is Spearman partial via partialcorr (paper convention).
%   Linear bisquare fit (matches rank-space Pearson = Spearman partial).
%   Toggle hippSide_fig: 'left' (primary) or 'bilateral' (diluted).
hippSide_fig = 'left';

ax_g = nexttile(19, [1,6]); set(ax_g,'Color',COL.bg); hold on;

idxFC_g = idxStructural;   % 3 = Diff(-)1 Ap-dom SM-DA (FDR-sig L hipp correlate)
longCSV_g = fullfile(outputDir, 'Table_Hippocampus_Longitudinal_Data.csv');

if exist(longCSV_g, 'file') == 2
    Tlong_g = readtable(longCSV_g);
    [~, longIdx_g] = ismember(fcPATNO, double(Tlong_g.PATNO));
    slopeV_g = nan(nFC, 1);
    vBLV_g   = nan(nFC, 1);
    if strcmp(hippSide_fig, 'bilateral')
        sCol_g = 'slope_bilateral'; bCol_g = 'vBL_bilateral'; sideLbl_g = 'Bilateral';
    else
        sCol_g = 'slope_left';      bCol_g = 'vBL_left';      sideLbl_g = 'Left';
    end
    maskM_g = longIdx_g > 0;
    slopeV_g(maskM_g) = Tlong_g.(sCol_g)(longIdx_g(maskM_g));
    vBLV_g(maskM_g)   = Tlong_g.(bCol_g)(longIdx_g(maskM_g));
    % Note: vBL/slope are HC-normative-corrected (eTIV implicitly handled
    % at the volume level); no explicit eTIV in Z_g. CSV may carry an
    % `eTIV_BL` column from the earlier uniform-method experiment — unused.

    fcV_g = fcData(:, idxFC_g);
    Z_g   = [fcMean.age, sexNum, fcMean.mean_FD, baselineEduc, vBLV_g];
    ok_g  = ~isnan(slopeV_g) & ~isnan(fcV_g) & all(~isnan(Z_g),2);

    % Rank-residualise to match the paper's partialcorr Spearman convention
    % used for panels (a, d, e) and now §13/§14 too.
    fcV_rk_g    = tiedrank(fcV_g(ok_g));
    slopeV_rk_g = tiedrank(slopeV_g(ok_g));
    Zrk_g = zeros(sum(ok_g), size(Z_g,2));
    for cc = 1:size(Z_g,2), Zrk_g(:,cc) = tiedrank(Z_g(ok_g,cc)); end
    Zrk_g_int = [ones(sum(ok_g),1), Zrk_g];   % intercept matches partialcorr
    fcR_g    = fcV_rk_g    - Zrk_g_int*(Zrk_g_int\fcV_rk_g);
    slopeR_g = slopeV_rk_g - Zrk_g_int*(Zrk_g_int\slopeV_rk_g);

    % Dark-blue sequential cmap (original)
    cmapSubtle_g = interp1(linspace(0,1,5), ...
        [0.92 0.92 0.94; 0.72 0.76 0.82; 0.48 0.56 0.68; 0.28 0.38 0.56; 0.14 0.22 0.42], ...
        linspace(0,1,256), 'pchip');
    densityScatter(ax_g, fcR_g, slopeR_g, 'ShowOLS', false, 'Cmap', cmapSubtle_g, ...
                   'SmootherType', 'robust', 'MarkerSize', 28, ...
                   'RobustColor', [0.80 0.40 0.00]);   % REVISED (R1.10): Okabe-Ito vermillion (was default red)

    % partialcorr = rank-everything Spearman partial (matches scatter + paper convention).
    [rho_g, p_g] = partialcorr(fcV_g, slopeV_g, Z_g, 'Type','Spearman','Rows','complete');

    % Diagnostic: compare raw-residual Pearson (what the LOESS line fits)
    % with the partialcorr Spearman (what the annotation shows) and the
    % partialcorr rank-residual Pearson (alternative Spearman-partial
    % convention used by panels d/e).
    try
        fc_ok   = fcV_g(ok_g);   sl_ok = slopeV_g(ok_g);   Z_ok = Z_g(ok_g,:);
        rP_raw  = partialcorr(fc_ok, sl_ok, Z_ok, 'Type','Pearson');
        fc_rk   = tiedrank(fc_ok); sl_rk = tiedrank(sl_ok);
        Z_rk    = zeros(sum(ok_g), size(Z_ok,2));
        for cc = 1:size(Z_ok,2), Z_rk(:,cc) = tiedrank(Z_ok(:,cc)); end
        fcR_rk  = fc_rk - Z_rk*(Z_rk\fc_rk);
        slR_rk  = sl_rk - Z_rk*(Z_rk\sl_rk);
        rP_rkrk = corr(fcR_rk, slR_rk, 'Type','Pearson');
        fprintf(['  [Panel g diag] Raw-Pearson=%+.3f  partialcorr-Spearman=%+.3f  ' ...
                 'RankResid-Pearson(=partialcorr)=%+.3f  n=%d\n'], ...
            rP_raw, rho_g, rP_rkrk, sum(ok_g));
    catch
    end

    text(ax_g, 0.03, 0.10, sprintf('\\rho = %+.2f,  p %s', rho_g, fmtP(p_g)), ...
        'Units','normalized','FontSize',6.5,'FontWeight','normal','VerticalAlignment','bottom', ...
        'Color',COL.txt,'BackgroundColor',[1 1 1 0.75], 'Margin', 1);
    text(ax_g, 0.03, 0.03, sprintf('n = %d (%s)', sum(ok_g), lower(sideLbl_g)), ...
        'Units','normalized','FontSize',6.5,'FontWeight','normal','VerticalAlignment','bottom', ...
        'Color',COL.txt,'BackgroundColor',[1 1 1 0.75], 'Margin', 1);

    xlabel('Apathy > ICB SM-DA FC (resid)', 'FontSize', 8);
    ylabel(sprintf('%s hippocampus \\Delta/yr (resid)', sideLbl_g), 'FontSize', 8);
else
    text(ax_g, 0.5, 0.5, {'Run Figure4\_Structural.m first', ...
        'to generate Table\_Hippocampus\_Longitudinal\_Data.csv'}, ...
        'Units','normalized','FontSize',7,'Color',[0.5 0.5 0.5], ...
        'HorizontalAlignment','center','VerticalAlignment','middle');
    set(ax_g, 'XTick', [], 'YTick', []);
end

title('\bfd\rm  FC predicts hippocampal atrophy', 'FontSize', 9, 'Color', COL.txt);
formatAxis(ax_g, COL);

% --- Axis LABELS + titles bold; tick NUMBERS left normal (matches Figure 2) ---
%     Runs last, so it also overrides the bold tick weight that formatAxis and
%     the axes default set earlier. Standalone text() annotations (stats, n)
%     keep their own weight.
%     NOTE: panel (c) ROI y-labels embed \bf in the label strings, so they stay
%     bold regardless (their bold/colour/italic encodes significance tier).
drawnow;
allAx = findall(fig4, 'Type', 'axes');
for kk = 1:numel(allAx)
    try
        allAx(kk).FontWeight = 'normal';   % axis tick numbers NOT bold
        if ~isempty(allAx(kk).XLabel), allAx(kk).XLabel.FontWeight = 'bold'; end
        if ~isempty(allAx(kk).YLabel), allAx(kk).YLabel.FontWeight = 'bold'; end
        if ~isempty(allAx(kk).Title),  allAx(kk).Title.FontWeight  = 'bold'; end
    catch
    end
end

%% ========================================================================
%  8.  EXPORT
% =========================================================================
exportFigure(fig4, fullfile(outputDir, 'Fig4_FC_Apathy_CogDecline'));
fprintf('\n  -> Fig4_FC_Apathy_CogDecline saved.\n');

%% ========================================================================
%  9.  CONSOLE SUMMARY
% =========================================================================
fprintf('\n%s\n', repmat('=',1,72));
fprintf('  FIGURE 4 SUMMARY (2x3 layout)\n');
fprintf('  Primary stats: %s partial r | %s\n', primaryMethod, strjoin(primaryCovLbl,', '));
fprintf('  Mediation: OLS bootstrap | FD + TIV (matches primary)\n');
fprintf('%s\n', repmat('-',1,72));
fprintf('  a. FC x Global Decline:  rho=%+.3f  p_FDR=%.3f %s\n', ...
    compR(idxMediation).rP, compR(idxMediation).pFDR, pToStars(compR(idxMediation).pFDR));
fprintf('  b. Mediation (%s): ab=%.4f [%.4f,%.4f] %s (%.0f%% mediated)\n', labelMed, ...
    ab, ab_ci(1), ab_ci(2), ternary(ab_sig,'SIG','n.s.'), prop);
fprintf('  c. Conditional trajectory: FC median-split longitudinal\n');
fprintf('  d. 3D surface: age x FC x cog slope\n');
if exist('datFC','var') && numel(datFC) >= idxMediation
    fprintf('  e. DAT -> FC:  rho=%+.3f  p_FDR=%.3f %s\n', ...
        datFC(idxMediation).rP, datFC(idxMediation).pFDR, pToStars(datFC(idxMediation).pFDR));
end
fprintf('  f. Domain trajectories: Exec/Mem/Visuo x Hi/Lo FC\n');
if exist('dynR','var') && size(dynR,2) >= 3
    fprintf('  g. FC x DAT slope:  rho=%+.3f  p=%.3f\n', ...
        dynR(idxMediation, 3), dynP(idxMediation, 3));
end
fprintf('  h. Structural correlates: %s\n', ...
    ternary(exist('normROI','var'), 'plotted', 'not available'));
fprintf('%s\n', repmat('=',1,72));

%% ========================================================================
%  10. LONGITUDINAL DYNAMICS: Apathy / ICB / DAT slopes × Baseline FC
%      Slopes from reduced LME: y ~ YEAR + covariates + (1+YEAR|PATNO)
%      using table T.  Correlate with baseline mean FC (partial | FD,TIV).
% =========================================================================
fprintf('\n[10] Longitudinal dynamics x baseline FC...\n');

assert(exist('TpdSub','var')==1, 'Table TpdSub required for longitudinal dynamics.');

% --- Build DAT mean column if needed ---
if ismember('MIA_CAUDATE_BILAT',TpdSub.Properties.VariableNames) && ...
        ismember('MIA_PUTAMEN_BILAT',TpdSub.Properties.VariableNames) && ...
        ~ismember('DAT_mean',TpdSub.Properties.VariableNames)
    TpdSub.DAT_mean = nanmean([double(TpdSub.MIA_CAUDATE_BILAT), double(TpdSub.MIA_PUTAMEN_BILAT)],2);
end

% --- Compute baseline cognition (time-invariant covariate) ---
gVarT = 'bifactor_g';
assert(ismember(gVarT, TpdSub.Properties.VariableNames), ...
    'TpdSub is missing column %s — bifactor factor scores must be merged before running Figure 4.', gVarT);
if ~ismember('g_baseline',TpdSub.Properties.VariableNames)
    TpdSub.g_baseline = nan(height(TpdSub),1);
    uPat = unique(double(TpdSub.PATNO));
    for pi = 1:numel(uPat)
        pIdx = double(TpdSub.PATNO)==uPat(pi);
        blIdx = pIdx & TpdSub.YEAR==0;
        if any(blIdx)
            TpdSub.g_baseline(pIdx) = nanmean(double(TpdSub.(gVarT)(blIdx)));
        end
    end
    fprintf('    Baseline g (%s): mapped for %d/%d subjects\n', ...
        gVarT, sum(~isnan(TpdSub.g_baseline(TpdSub.YEAR==0))), numel(uPat));
end

% --- Compute baseline apathy & ICB (time-invariant cross-covariates) ---
if ~exist('uPat','var'), uPat = unique(double(TpdSub.PATNO)); end
for blVar = {'apathy_composite','quip'}
    vn = blVar{1};  blName = [vn '_bl'];
    if ismember(vn,TpdSub.Properties.VariableNames) && ~ismember(blName,TpdSub.Properties.VariableNames)
        TpdSub.(blName) = nan(height(TpdSub),1);
        for pi = 1:numel(uPat)
            pIdx = double(TpdSub.PATNO)==uPat(pi);
            blIdx = pIdx & TpdSub.YEAR==0;
            if any(blIdx)
                TpdSub.(blName)(pIdx) = nanmean(double(TpdSub.(vn)(blIdx)));
            end
        end
        fprintf('    Baseline %-20s: %d/%d subjects\n', blName, ...
            sum(~isnan(TpdSub.(blName)(TpdSub.YEAR==0))), numel(uPat));
    end
end

dynVars   = {'apathy_composite','quip','DAT_mean'};
dynLabels = {'Apathy','ICB (QUIP)','DAT Binding'};
% Cross-phenotype baselines (quip_bl, apathy_composite_bl, g_baseline)
% INTENTIONALLY DROPPED — keeping them in slope extraction would
% double-residualise the same constructs that the NBS-matched partial
% correlation already adjusts for, and over-strip the FC-relevant
% trait-level signal. Single-stage adjustment via NBS-matched covariates
% in the partial correlation step is sufficient.
dynExtra  = {'','',''};
keep = cellfun(@(v) ismember(v,TpdSub.Properties.VariableNames), dynVars);
dynVars = dynVars(keep);  dynLabels = dynLabels(keep);  dynExtra = dynExtra(keep);
nDyn = numel(dynVars);
fprintf('    Variables: %s\n', strjoin(dynLabels,', '));

% --- Covariates in T (g_baseline DROPPED — see comment above) ---
covCandsT = {'EDUCYRS','age','SEX','duration','depression_factor','LEDD','updrs3_score'};
covPresentT = unique(covCandsT(cellfun(@(v) ismember(v,TpdSub.Properties.VariableNames), covCandsT)), 'stable');
covStrT = ''; if ~isempty(covPresentT), covStrT = [' + ' strjoin(covPresentT,' + ')]; end

% --- Z-score numeric covariates (working copy, safe for re-runs) ---
TpdZ = TpdSub;
allCovToZ = unique([covPresentT, dynExtra(~cellfun(@isempty,dynExtra))], 'stable');
nZscored = 0;
for ci = 1:numel(allCovToZ)
    vn = allCovToZ{ci};
    if ~ismember(vn,TpdZ.Properties.VariableNames), continue; end
    if iscategorical(TpdZ.(vn)), continue; end
    vals = double(TpdZ.(vn));
    ok_z = ~isnan(vals);
    vals(ok_z) = zscore(vals(ok_z));
    TpdZ.(vn) = vals;
    nZscored = nZscored + 1;
end
fprintf('    Z-scored %d/%d covariates\n', nZscored, numel(allCovToZ));

TPATNO = double(TpdSub.PATNO);
okBaseT = ~isnan(double(TpdSub.YEAR));
for ci = 1:numel(covPresentT)
    vn = covPresentT{ci};
    if iscategorical(TpdSub.(vn)), okBaseT = okBaseT & ~isundefined(TpdSub.(vn));
    else, okBaseT = okBaseT & ~isnan(double(safeCol(TpdSub,vn))); end
end

% --- Extract slopes (LME + OLS fallback). Winsor copy saved as sensitivity. ---
dynSlopes = nan(nFC, nDyn);
for vi = 1:nDyn
    okD = okBaseT & ~isnan(double(safeCol(TpdSub, dynVars{vi})));
    % Add cross-covariate if available
    extraStr = '';
    if ~isempty(dynExtra{vi}) && ismember(dynExtra{vi}, TpdSub.Properties.VariableNames)
        okD = okD & ~isnan(double(safeCol(TpdSub, dynExtra{vi})));
        extraStr = [' + ' dynExtra{vi}];
    end
    fmla = sprintf('%s ~ YEAR%s%s + (1 + YEAR | PATNO)', dynVars{vi}, covStrT, extraStr);
    fprintf('    Formula: %s  (n=%d)\n', fmla, sum(okD));
    fprintf('    Covariates found: %s\n', strjoin(covPresentT,', '));
    try
        mdlD = fitlme(TpdZ(okD,:), fmla, 'FitMethod','REML');
        [sl,nM] = extractBLUPsFromModel(mdlD, fcPATNO, nFC);
        dynSlopes(:,vi) = sl;
        fprintf('    %-15s  AIC=%7.1f  matched=%d\n', dynLabels{vi}, mdlD.ModelCriterion.AIC, nM);
    catch ME
        fprintf('    %-15s  LME FAILED: %s\n', dynLabels{vi}, ME.message);
    end
end
% OLS fallback
for si = 1:nFC
    sIdx = TPATNO == fcPATNO(si); if sum(sIdx)<2, continue; end
    tV = double(TpdSub.YEAR(sIdx));
    for vi = 1:nDyn
        if ~isnan(dynSlopes(si,vi)), continue; end
        yV = double(safeCol(TpdSub(sIdx,:), dynVars{vi})); ok = ~isnan(tV)&~isnan(yV);
        if sum(ok)>=2, b=[ones(sum(ok),1),tV(ok)]\yV(ok); dynSlopes(si,vi)=b(2); end
    end
end
% Sensitivity-only: winsorised copy at 1/99 (dynSlopes itself stays raw)
dynSlopesWinsor = dynSlopes;
for vi = 1:nDyn
    v=dynSlopesWinsor(:,vi); ok=~isnan(v); if sum(ok)<10, continue; end
    p1=prctile(v(ok),1); p99=prctile(v(ok),99);
    v(ok&v<p1)=p1; v(ok&v>p99)=p99; dynSlopesWinsor(:,vi)=v;
end

% --- Dynamics covariate set: matches primary (FD + TIV).
%     The richer set [age, sex, apathy_bl, QUIP_bl, DAT_bl, dep_bl] is now
%     a sensitivity-only adjustment; primary stays at FD+TIV per the
%     paper-wide spec at L442.
if hasDAT
    dat_bl = fcMean.DAT_mean;
else
    dat_bl = nan(height(fcMean),1);
end
dynCovs   = [fcMean.mean_FD, fcMean.TIV];
dynCovLbl = {'FD','TIV'};

% --- Outcome-specific baseline (computed but NOT used in primary) -------
%     Empirical Corr (BLUP slope, baseline value):  DAT: -0.78  apathy: +0.23
%     Per user instruction (2026-05-06), baseline is NOT appended to the
%     primary covariate set for §10 dynamics. The variable is kept here for
%     potential sensitivity analyses but is unused in the primary loop below.
dynBaseline = nan(nFC, nDyn);
dynBaseLbl  = cell(1, nDyn);
for vi = 1:nDyn
    switch lower(dynVars{vi})
        case 'apathy_composite'
            dynBaseline(:,vi) = fcMean.apathy_composite;
            dynBaseLbl{vi} = 'apathy_bl';
        case 'quip'
            dynBaseline(:,vi) = baselineQUIP;
            dynBaseLbl{vi} = 'QUIP_bl';
        case 'dat_mean'
            dynBaseline(:,vi) = dat_bl;
            dynBaseLbl{vi} = 'DAT_bl';
        otherwise
            dynBaseLbl{vi} = '';
    end
end
fprintf('    Dynamics partial-correlation covariates: %s (no per-outcome baseline)\n', ...
    strjoin(dynCovLbl,', '));
fprintf('    Subjects complete on FD/TIV: %d / %d\n', sum(all(~isnan(dynCovs),2)), nFC);

% --- Correlate with baseline FC (primary: Spearman | dynCovs only, no outcome baseline) ---
dynR = nan(nAllFC, nDyn);        dynP = nan(nAllFC, nDyn);
dynR_basic = nan(nAllFC, nDyn);  dynP_basic = nan(nAllFC, nDyn);
for fi = 1:nAllFC
    fcV = fcData(:, fi);
    for vi = 1:nDyn
        % Primary: FD + TIV only (no per-outcome baseline; user spec 2026-05-06)
        covs_v = dynCovs;
        [rP_prim, pP_prim] = partialcorr(fcV, dynSlopes(:,vi), covs_v, 'Type', primaryMethod, 'Rows', 'complete');
        [~,~,rP_basic,pP_basic] = fcSlopeCorr(fcV, dynSlopes(:,vi), ...
                                              fcMean.mean_FD, fcMean.TIV);
        dynR(fi,vi) = rP_prim;        dynP(fi,vi) = pP_prim;
        dynR_basic(fi,vi) = rP_basic; dynP_basic(fi,vi) = pP_basic;
    end
end

% Hypothesis-driven FDR: apathy FC components x dynamic slopes only.
% Non-apathy (diff/shared) FC remain uncorrected exploratory.
dynPFDR = nan(nAllFC, nDyn);
P_ap = dynP(apathyIdx, :);
fdr_ap = fdr_bh(P_ap(:));
dynPFDR(apathyIdx, :) = reshape(fdr_ap, size(P_ap));
nPrimaryTests = numel(P_ap);

% --- Print results ---
% Build a human-readable cov string per outcome
covsByOut = repmat({strjoin(dynCovLbl,'+')}, 1, nDyn);
fprintf('\n  Primary %s partial r | dynamics covariates per outcome:\n', primaryMethod);
for vi = 1:nDyn
    fprintf('    %-15s : %s\n', dynLabels{vi}, covsByOut{vi});
end
fprintf('  (FDR across %d apathy-network tests; non-apathy uncorrected.)\n', nPrimaryTests);
fprintf('  %-25s', ''); for vi=1:nDyn, fprintf('  %14s', dynLabels{vi}); end
fprintf('\n  %s\n', repmat('-',1,25+15*nDyn));
for fi = 1:nAllFC
    fprintf('  %-25s', fcShort{fi});
    for vi = 1:nDyn
        star='';
        if ~isnan(dynPFDR(fi,vi)) && dynPFDR(fi,vi)<0.05, star='+';
        elseif dynP(fi,vi)<0.05, star='*'; end
        fprintf('  %+.3f%-3s      ', dynR(fi,vi), star);
    end
    fprintf('\n');
end

fprintf('\n  FDR-significant (apathy-family, %d tests):\n', nPrimaryTests);
hasSig = false;
for fi=1:nAllFC, for vi=1:nDyn
    if ~isnan(dynPFDR(fi,vi)) && dynPFDR(fi,vi)<0.05
        fprintf('    %-25s x %-15s r_p=%+.3f  p_FDR=%.4f\n', ...
            fcShort{fi}, dynLabels{vi}, dynR(fi,vi), dynPFDR(fi,vi));
        hasSig = true;
    end
end, end
if ~hasSig, fprintf('    (none)\n'); end

% --- Save ---
rows = {};
for fi=1:nAllFC, for vi=1:nDyn
    rows{end+1} = {fcLabels{fi},dynLabels{vi},dynR(fi,vi),dynP(fi,vi),dynPFDR(fi,vi)}; %#ok
end, end
% writetable(cell2table(vertcat(rows{:}), ...
%     'VariableNames',{'FC_Component','Slope_Variable','r_partial','p_unc','p_FDR'}), ...
%     fullfile(outputDir,'Table_Fig4_LongDynamics_x_FC.csv'));
% fprintf('    -> Table_Fig4_LongDynamics_x_FC.csv saved.\n');

%% ========================================================================
%  10b. LONGITUDINAL DYNAMICS — age-moderation analysis
%      Mirrors Section 14 sensitivity structure for the dynamics tests.
%        (A) Age-adjusted partial r (FD + TIV + age)
%        (B) Spearman partial r (same covariates)
%        (C) Age-stratified FC x slope
%        (D) FC x age-group interaction on dynamic slope
%      Hypothesis-driven: focuses on apathyIdx x 3 dynamic variables.
% =========================================================================
fprintf('\n[10b] Age-moderation of baseline FC x longitudinal slopes...\n');

% Use the same median split as Figure 4 Panel A/B
if ~exist('ageMed_fig','var'), ageMed_fig = nanmedian(fcMean.age); end
ageVecDyn  = fcMean.age;
strataDyn  = {'Younger (<=med)', ageVecDyn <= ageMed_fig; ...
              'Older  (>med)',   ageVecDyn >  ageMed_fig};
ageGrpDyn  = double(ageVecDyn > ageMed_fig);   % 1=older
Z_base_dyn = [fcMean.mean_FD, fcMean.TIV];
Z_ageAdj   = [fcMean.mean_FD, fcMean.TIV, fcMean.age];

% ----- (A) + (B): Partial r | dynCovs, Pearson and Spearman --------------
%      (dynCovs = FD + TIV + age + sex + apathy_bl + QUIP_bl + DAT_bl + dep_bl)
fprintf('\n  (A/B) Partial r | dynCovs (%s)   (Pearson | Spearman)\n', strjoin(dynCovLbl,', '));
fprintf('  %-22s', '');
for vi=1:nDyn, fprintf('  %20s', dynLabels{vi}); end
fprintf('\n  %s\n', repmat('-',1,22 + 22*nDyn));

dynR_age  = nan(nAllFC, nDyn);  dynP_age  = nan(nAllFC, nDyn);
dynR_sp   = nan(nAllFC, nDyn);  dynP_sp   = nan(nAllFC, nDyn);
for fi = 1:nAllFC
    fcV = fcData(:, fi);
    fprintf('  %-22s', fcShort{fi});
    for vi = 1:nDyn
        [rP, pP] = partialcorr(fcV, dynSlopes(:,vi), dynCovs, 'Type','Pearson','Rows','complete');
        [rS, pS] = partialcorr(fcV, dynSlopes(:,vi), dynCovs, 'Type','Spearman','Rows','complete');
        dynR_age(fi,vi) = rP; dynP_age(fi,vi) = pP;
        dynR_sp(fi,vi)  = rS; dynP_sp(fi,vi)  = pS;
        starP = '';
        if ~isnan(pP) && pP < 0.05, starP = '*'; end
        starS = '';
        if ~isnan(pS) && pS < 0.05, starS = '*'; end
        fprintf('  %+.2f%s / %+.2f%s ', rP, starP, rS, starS);
    end
    fprintf('\n');
end

% FDR on apathy family for the age-adjusted and Spearman matrices
dynPFDR_age = nan(nAllFC, nDyn);
dynPFDR_sp  = nan(nAllFC, nDyn);
Pap_age = dynP_age(apathyIdx, :);  dynPFDR_age(apathyIdx, :) = reshape(fdr_bh(Pap_age(:)), size(Pap_age));
Pap_sp  = dynP_sp(apathyIdx, :);   dynPFDR_sp(apathyIdx, :)  = reshape(fdr_bh(Pap_sp(:)),  size(Pap_sp));

fprintf('\n  FDR-significant under age-adjusted Pearson (apathy family):\n');
anySig = false;
for fi = apathyIdx(:)'
    for vi = 1:nDyn
        if ~isnan(dynPFDR_age(fi,vi)) && dynPFDR_age(fi,vi) < 0.05
            fprintf('    %-22s x %-15s  r=%+.3f  p_FDR=%.4f\n', ...
                fcShort{fi}, dynLabels{vi}, dynR_age(fi,vi), dynPFDR_age(fi,vi));
            anySig = true;
        end
    end
end
if ~anySig, fprintf('    (none)\n'); end

fprintf('  FDR-significant under age-adjusted Spearman (apathy family):\n');
anySig = false;
for fi = apathyIdx(:)'
    for vi = 1:nDyn
        if ~isnan(dynPFDR_sp(fi,vi)) && dynPFDR_sp(fi,vi) < 0.05
            fprintf('    %-22s x %-15s  r=%+.3f  p_FDR=%.4f\n', ...
                fcShort{fi}, dynLabels{vi}, dynR_sp(fi,vi), dynPFDR_sp(fi,vi));
            anySig = true;
        end
    end
end
if ~anySig, fprintf('    (none)\n'); end

% ----- (C) Age-stratified FC x slope (Pearson and Spearman) -------------
%      Within-stratum covariates: dynCovs minus age (we are stratifying on age).
dynCovs_noAge = dynCovs(:, ~strcmp(dynCovLbl,'age'));
dynCovLbl_noAge = dynCovLbl(~strcmp(dynCovLbl,'age'));
fprintf('\n  (C) Age-stratified Pearson / Spearman partial r | %s:\n', strjoin(dynCovLbl_noAge,', '));
dynRstrat  = cell(size(strataDyn,1),1);
dynPstrat  = cell(size(strataDyn,1),1);
dynRstratS = cell(size(strataDyn,1),1);
dynPstratS = cell(size(strataDyn,1),1);
for si = 1:size(strataDyn,1)
    sel = strataDyn{si,2};
    R = nan(nAllFC, nDyn); P = nan(nAllFC, nDyn);
    Rs = R; Ps = P;
    for fi = 1:nAllFC
        fcV = fcData(:, fi);
        for vi = 1:nDyn
            x = fcV(sel); y = dynSlopes(sel, vi);
            Z = dynCovs_noAge(sel,:);
            [rP, pP] = partialcorr(x, y, Z, 'Type','Pearson','Rows','complete');
            [rS, pS] = partialcorr(x, y, Z, 'Type','Spearman','Rows','complete');
            R(fi,vi) = rP; P(fi,vi) = pP;
            Rs(fi,vi) = rS; Ps(fi,vi) = pS;
        end
    end
    dynRstrat{si}  = R;  dynPstrat{si}  = P;
    dynRstratS{si} = Rs; dynPstratS{si} = Ps;

    fprintf('  --- %s (n=%d) ---\n', strataDyn{si,1}, sum(sel));
    fprintf('  %-22s', '');
    for vi=1:nDyn, fprintf('  %20s', dynLabels{vi}); end
    fprintf('\n  %s\n', repmat('-',1,22 + 22*nDyn));
    for fi = apathyIdx(:)'
        fprintf('  %-22s', fcShort{fi});
        for vi = 1:nDyn
            starP = '';
            if ~isnan(P(fi,vi))  && P(fi,vi)  < 0.05, starP = '*'; end
            starS = '';
            if ~isnan(Ps(fi,vi)) && Ps(fi,vi) < 0.05, starS = '*'; end
            fprintf('  %+.2f%s / %+.2f%s ', R(fi,vi), starP, Rs(fi,vi), starS);
        end
        fprintf('\n');
    end
end

% ----- (D) FC x age-group interaction on dynamic slopes -----------------
%      Now controlling for baseline apathy, baseline QUIP, baseline DAT, dep
%      (in addition to FD, TIV, sex). Age enters via the binary AgeGrp.
fprintf('\n  (D) FC x age-group interaction on dynamic slopes (apathy comps)\n');
fprintf('       | covariates: FD, TIV, sex, apathy_bl, QUIP_bl, DAT_bl, dep_bl + AgeGrp\n');
fprintf('  %-22s  %-15s  %-18s  %-8s\n', ...
    'FC Component', 'Slope', 'FC_x_Age β', 'p-value');
fprintf('  %s\n', repmat('-',1,70));
intRows = {};
for ci = 1:nApFC
    fi = apathyIdx(ci);
    fcV = fcData(:, fi);
    for vi = 1:nDyn
        y = dynSlopes(:, vi);
        X = [fcV, fcV.*ageGrpDyn, double(ageGrpDyn), ...
             fcMean.mean_FD, fcMean.TIV, sexNum, ...
             fcMean.apathy_composite, baselineQUIP, dat_bl, baselineDEP];
        okJ = ~isnan(y) & all(~isnan(X),2);
        if sum(okJ) < 30
            fprintf('  %-22s  %-15s  n<30, skipped\n', fcShort{fi}, dynLabels{vi});
            continue;
        end
        mdlInt_dyn = fitlm(X(okJ,:), y(okJ), ...
            'VarNames', {'FC','FC_x_Age','AgeGrp','FD','TIV','sex', ...
                         'apathy_bl','QUIP_bl','DAT_bl','dep_bl','slope'});
        intRow = strcmp(mdlInt_dyn.CoefficientNames, 'FC_x_Age');
        bInt = mdlInt_dyn.Coefficients.Estimate(intRow);
        pInt = mdlInt_dyn.Coefficients.pValue(intRow);
        fprintf('  %-22s  %-15s  %+7.3f (±%.3f)   %.4f %s\n', ...
            fcShort{fi}, dynLabels{vi}, bInt, ...
            mdlInt_dyn.Coefficients.SE(intRow), pInt, pToStars(pInt));
        intRows{end+1} = {fcShort{fi}, dynLabels{vi}, bInt, ...
            mdlInt_dyn.Coefficients.SE(intRow), pInt}; %#ok
    end
end
if ~isempty(intRows)
    writetable(cell2table(vertcat(intRows{:}), ...
        'VariableNames', {'FC_Component','Slope','FCxAge_beta','FCxAge_SE','FCxAge_p'}), ...
        fullfile(outputDir,'Table_Fig4_Dyn_FCxAge_Interactions.csv'));
    fprintf('  -> Table_Fig4_Dyn_FCxAge_Interactions.csv saved.\n');
end

% Save age-adjusted dynamics table
dynAgeRows = {};
for fi = 1:nAllFC
    for vi = 1:nDyn
        dynAgeRows{end+1} = {fcLabels{fi}, dynLabels{vi}, ...
            dynR_age(fi,vi), dynP_age(fi,vi), dynPFDR_age(fi,vi), ...
            dynR_sp(fi,vi),  dynP_sp(fi,vi),  dynPFDR_sp(fi,vi)}; %#ok
    end
end
writetable(cell2table(vertcat(dynAgeRows{:}), ...
    'VariableNames', {'FC_Component','Slope','Pearson_r','Pearson_p','Pearson_pFDR', ...
                      'Spearman_r','Spearman_p','Spearman_pFDR'}), ...
    fullfile(outputDir,'Table_Fig4_Dyn_AgeAdjusted.csv'));
fprintf('  -> Table_Fig4_Dyn_AgeAdjusted.csv saved.\n');

%% ========================================================================
%  11. SUPPLEMENTARY FIGURE S1 — FC Controls & Confirmations (2x2)
%      (A) Secondary FC (idxB) x Global      (B) Path a: Apathy -> primary FC (idxA)
%      (C) DAT sensitivity (apathy comps)    (D) Full domain heatmap (4 domains)
% =========================================================================
fprintf('\n[11] Supplementary Figure S1...\n');

figS1 = figure('Units','centimeters','Position',[1 1 26 8],'Color',COL.bg);
set(figS1,'DefaultAxesFontName','Helvetica','DefaultTextFontName','Helvetica');
% Layout switched from 2×2 to 1×3 after the domain-heatmap tile (former
% S1D) moved to Figure4_Supp_Domain.m. Three FC-control panels remain.
tiledlayout(1, 3, 'TileSpacing','compact', 'Padding','compact');

% Color for ±DaT comparison markers in panel S1C (was previously inherited
% from the deleted main-figure panel a; defined locally now).
datCol = [0.55 0.30 0.55];   % muted purple

% ----- S1A) Broad apathy hyperconn (idxA) x Global Slope ----------------
%     The non-mediation apathy FC signature (displaced from main figure).
ax = nexttile(1); set(ax,'Color',COL.bg); hold on;
if idxA <= nAllFC && ~isempty(compR(idxA).x)
    densityScatter(ax, compR(idxA).x, compR(idxA).y);
    yline(ax, 0, '--', 'Color', COL.grid, 'LineWidth', 0.5);
    text(ax, 0.03, 0.97, { ...
        sprintf('r_{partial} = %.2f', compR(idxA).rP), ...
        sprintf('p_{FDR} = %.3f %s', compR(idxA).pFDR, pToStars(compR(idxA).pFDR)), ...
        sprintf('n = %d', compR(idxA).n)}, ...
        'Units','normalized','FontSize',6,'VerticalAlignment','top', ...
        'Color',COL.txt,'BackgroundColor',[1 1 1 0.75]);
end
xlabel(labelA_axis, 'FontSize', 7);
ylabel('Global Cog Slope (yr^{-1})', 'FontSize', 7);
title(sprintf('\\bfA\\rm  %s \\times Global Decline', labelA), 'FontSize', 8, 'Color', COL.txt);
formatAxis(ax, COL);

% ----- S1B) Path a: Apathy -> Primary FC (idxA) --------------------------
ax = nexttile(2); set(ax,'Color',COL.bg); hold on;
if ~isempty(pathA(idxA).x)
    densityScatter(ax, pathA(idxA).x, pathA(idxA).y);
    text(ax, 0.03, 0.97, { ...
        sprintf('r_{partial} = %.2f', pathA(idxA).rP), ...
        sprintf('p = %.4f %s', pathA(idxA).pP, pToStars(pathA(idxA).pP)), ...
        sprintf('n = %d', pathA(idxA).n)}, ...
        'Units','normalized','FontSize',6,'VerticalAlignment','top', ...
        'Color',COL.txt,'BackgroundColor',[1 1 1 0.75]);
end
xlabel('Apathy Composite', 'FontSize', 7);
ylabel(labelA_axis, 'FontSize', 7);
title(sprintf('\\bfB\\rm  Path a: Apathy \\rightarrow %s', labelA), 'FontSize', 8, 'Color', COL.txt);
formatAxis(ax, COL);

% ----- S1C) DAT Sensitivity: paired dots +/-DAT --------------------------
%     datSens is now indexed by fi; only apathyIdx entries are populated.
ax = nexttile(3); set(ax,'Color',COL.bg); hold on;
if hasDAT
    for ci = 1:nApFC
        fi = apathyIdx(ci);
        yp = nApFC - ci + 1;
        r0 = datSens(fi).rP0;  rd = datSens(fi).rPd;
        isSig0 = ~isnan(datSens(fi).pP0) && datSens(fi).pP0 < 0.05;
        colF = scCol(min(fi,size(scCol,1)),:);
        if isSig0
            plot(ax, r0, yp+0.12, 'o', 'Color',colF, 'MarkerFaceColor',colF, 'MarkerSize',5, 'LineWidth',0.5);
        else
            plot(ax, r0, yp+0.12, 'o', 'Color',colF, 'MarkerFaceColor',COL.bg, 'MarkerSize',5, 'LineWidth',1);
        end
        isSigD = ~isnan(datSens(fi).pPd) && datSens(fi).pPd < 0.05;
        if isSigD
            plot(ax, rd, yp-0.12, 'd', 'Color',datCol, 'MarkerFaceColor',datCol, 'MarkerSize',5, 'LineWidth',0.5);
        else
            plot(ax, rd, yp-0.12, 'd', 'Color',datCol, 'MarkerFaceColor',COL.bg, 'MarkerSize',5, 'LineWidth',1);
        end
        plot(ax, [r0 rd], [yp+0.12 yp-0.12], '-', 'Color',[0.75 0.75 0.75], 'LineWidth',0.5);
    end
    xline(ax, 0, '-', 'Color', [0.78 0.78 0.78], 'LineWidth', 0.5);
    set(ax, 'YTick',1:nApFC, 'YTickLabel',flip(fcShort(apathyIdx)), ...
        'YLim',[0.4 nApFC+0.6], 'FontSize',6, 'TickDir','out', 'Box','off');
    xlabel('r_{partial} (FC \times Global Slope)', 'FontSize', 7);
    hO = plot(NaN,NaN,'o','Color',[0.4 0.4 0.4],'MarkerFaceColor',[0.4 0.4 0.4],'MarkerSize',4);
    hD = plot(NaN,NaN,'d','Color',datCol,'MarkerFaceColor',datCol,'MarkerSize',4);
    lg = legend([hO hD], {'FD + TIV','+ DaT'}, 'FontSize',5,'Box','off','Location','southeast');
    lg.ItemTokenSize = [8 8];
end
title('\bfC\rm  Sensitivity: \pm DaT Control', 'FontSize', 8, 'Color', COL.txt);
formatAxis(ax, COL);

% (Former S1D — Full Domain Heatmap — is now produced by
%  Figure4_Supp_Domain.m as a standalone 2×2 figure with 4 LME × FCpartial
%  variants and family-wise FDR. See Fig4_Supp_DomainSpecificity.pdf.)

exportFigure(figS1, fullfile(outputDir, 'SuppFig_S1_FC_Controls'));
fprintf('  -> SuppFig_S1_FC_Controls saved.\n');

%% ========================================================================
%  12. SUPPLEMENTARY FIGURE S2 — Longitudinal Dynamics Heatmap
% =========================================================================
fprintf('\n[12] Supplementary Figure S2...\n');

if exist('dynR','var') && ~isempty(dynR)
    figS2 = figure('Units','centimeters','Position',[1 1 14 12],'Color',COL.bg);
    set(figS2,'DefaultAxesFontName','Helvetica','DefaultTextFontName','Helvetica');
    ax = axes('Parent',figS2); set(ax,'Color',COL.bg);

    rHd = dynR; rHd(isnan(rHd))=0;
    mabd = max(abs(rHd(:))); if mabd<0.1, mabd=0.3; end
    imagesc(ax, rHd, [-1 1]*mabd*1.1);
    if ~exist('cmap','var') || isempty(cmap)   % guard
        cmap = [[linspace(0,1,128)'; ones(128,1)], [linspace(0.35,1,128)'; linspace(1,0.35,128)'], [ones(128,1); linspace(1,0,128)']];
    end
    colormap(ax, cmap); cb=colorbar(ax); cb.Label.String='r_{partial}'; cb.FontSize=6;

    hold(ax,'on');
    for fi=1:nAllFC
        for vi=1:nDyn
            rP=dynR(fi,vi); pP=dynP(fi,vi); pF=dynPFDR(fi,vi);
            if isnan(rP), continue; end
            tc=ternary(abs(rP)>mabd*0.55,[1 1 1],COL.txt);
            str=sprintf('%.2f',rP);
            if ~isnan(pF)&&pF<0.05
                str=[str '+']; %#ok
                rectangle(ax,'Position',[vi-0.5,fi-0.5,1,1],'EdgeColor',[0.1 0.1 0.1],'LineWidth',1.8);
            elseif pP<0.05
                str=[str '*']; %#ok
            end
            text(ax,vi,fi,str,'HorizontalAlignment','center','FontSize',5.5,'Color',tc);
        end
    end

    nAp_hm=sum(fcContrast==1); nICB_hm=sum(fcContrast==2); nInt_hm=sum(fcContrast==3);
    yline(ax, nAp_hm+0.5, '-', 'Color',COL.txt, 'LineWidth',0.8);
    yline(ax, nAp_hm+nICB_hm+0.5, '-', 'Color',COL.txt, 'LineWidth',0.8);
    text(ax, nDyn+0.55, (1+nAp_hm)/2, 'Ap', 'FontSize',5,'Color',COL.apathy,'FontWeight','bold','HorizontalAlignment','left');
    text(ax, nDyn+0.55, nAp_hm+(1+nICB_hm)/2, 'ICB', 'FontSize',5,'Color',COL.icb,'FontWeight','bold','HorizontalAlignment','left');
    text(ax, nDyn+0.55, nAp_hm+nICB_hm+(1+nInt_hm)/2, 'Int', 'FontSize',5,'Color',COL.overlap,'FontWeight','bold','HorizontalAlignment','left');
    set(ax,'XTick',1:nDyn,'XTickLabel',dynLabels,'XTickLabelRotation',30,...
        'YTick',1:nAllFC,'YTickLabel',fcShort,'FontSize',6,'TickDir','out','Box','off');
    title('Baseline FC \times Longitudinal Slopes','FontSize',9,'Color',COL.txt);

    exportFigure(figS2, fullfile(outputDir, 'SuppFig_S2_LongDynamics'));
    fprintf('  -> SuppFig_S2_LongDynamics saved.\n');
else
    fprintf('  Skipped (no longitudinal dynamics data).\n');
end

%% ========================================================================
%  13. INSPECTION FIGURE — Density-coloured scatter + OLS + robust regression
%      NOT for publication. Diagnostic view of residualized partial scatter
%      for every FC × domain combination.
%
%      Each scatter shows (x, y) residualized against FD, TIV — i.e. the
%      partial-correlation geometry, so the slope you see matches r_partial.
%      - Point colour   = 2D Gaussian KDE density (parula)
%      - Blue dashed    = OLS fit
%      - Red solid      = robust bisquare fit (robustfit)
%      - Spearman r     = non-parametric sanity check
%      Stars: + = FDR<.05, * = p_unc<.05.
% =========================================================================
fprintf('\n[13] Building inspection figure...\n');

inspectIdx = 1:nAllFC;   % change to apathyIdx or primaryIdx to restrict
nInsp = numel(inspectIdx);

figW = max(18, 6*nDom + 3);
figH = max(14, 4.5*nInsp + 3);
figInsp = figure('Units','centimeters', ...
    'Position',[1 1 figW figH], ...
    'Color',COL.bg, 'Name','Inspection: FC x Cognitive Slopes');
set(figInsp,'DefaultAxesFontName','Helvetica','DefaultTextFontName','Helvetica');
tI = tiledlayout(nInsp, nDom, 'TileSpacing','compact','Padding','compact');

for ri = 1:nInsp
    fi = inspectIdx(ri);
    fcV = fcData(:, fi);
    for di = 1:nDom
        ax = nexttile; set(ax,'Color',COL.bg); hold on;

        sV = slopeAll(:, di);
        Z  = [fcMean.mean_FD, fcMean.TIV];
        ok = ~isnan(fcV) & ~isnan(sV) & all(~isnan(Z),2);
        n  = sum(ok);

        if n < 10
            text(ax, 0.5, 0.5, sprintf('n=%d', n), ...
                'Units','normalized','HorizontalAlignment','center', ...
                'FontSize',7,'Color',[0.5 0.5 0.5]);
            set(ax,'XTick',[],'YTick',[]);
            if ri == 1, title(ax, cogLabels{di}, 'FontSize', 8, 'Color', COL.txt); end
            if di == 1, ylabel(ax, fcShort{fi}, 'FontSize', 6, 'FontWeight','bold'); end
            continue;
        end

        % --- Residualize against FD, TIV (partial-correlation geometry) ---
        Zo = Z(ok,:);
        x  = fcV(ok) - Zo*(Zo\fcV(ok));
        y  = sV(ok)  - Zo*(Zo\sV(ok));

        % --- 2D density estimate at each point (Gaussian KDE) ---
        dens = nan(size(x));
        try
            dens = mvksdensity([x y], [x y], ...
                'Bandwidth', [std(x)*0.35, std(y)*0.35]);
        catch
            % Fallback: bin-based density
            try
                [N, xE, yE] = histcounts2(x, y, 20, 20);
                xb = max(min(discretize(x,xE), size(N,1)), 1);
                yb = max(min(discretize(y,yE), size(N,2)), 1);
                for i = 1:numel(x), dens(i) = N(xb(i),yb(i)); end
            catch
                dens = ones(size(x));
            end
        end

        scatter(ax, x, y, 18, dens, 'filled', ...
            'MarkerFaceAlpha',0.75, 'MarkerEdgeColor','none');
        colormap(ax, parula);

        % --- OLS and robust regression lines ---
        xF = linspace(min(x), max(x), 100);
        b_ols = [ones(n,1), x] \ y;
        plot(ax, xF, b_ols(1)+b_ols(2)*xF, '--', ...
            'Color', [0.15 0.35 0.70], 'LineWidth', 1.3);

        b_rob = [NaN; NaN];
        try
            b_rob = robustfit(x, y);  % bisquare weight (default)
            plot(ax, xF, b_rob(1)+b_rob(2)*xF, '-', ...
                'Color', [0.80 0.20 0.20], 'LineWidth', 1.6);
        catch
        end

        rSp = corr(x, y, 'Type','Spearman');

        yline(ax, 0, ':', 'Color', COL.grid, 'LineWidth', 0.5);
        xline(ax, 0, ':', 'Color', COL.grid, 'LineWidth', 0.5);

        % --- Stats overlay -------------------------------------------------
        % Global (di=1) reuses the §4 main-analysis stats (rP, pP, pFDR).
        % Subdomains are diagnostic only (the formal FDR family for
        % subdomains lives in Figure4_Supp_Domain.m) — compute a
        % Spearman partial inline, omit FDR.
        if di == 1
            rP   = compR(fi).rP;
            pP   = compR(fi).pP;
            pFDR = compR(fi).pFDR;
        else
            sV_diag = slopeAll(:, di);
            Z_diag  = [fcMean.mean_FD, fcMean.TIV];
            okI = ~isnan(fcV) & ~isnan(sV_diag) & all(~isnan(Z_diag), 2);
            if sum(okI) >= 10
                [rP, pP] = partialcorr(fcV, sV_diag, Z_diag, ...
                    'Type', primaryMethod, 'Rows', 'complete');
            else
                rP = NaN; pP = NaN;
            end
            pFDR = NaN;
        end
        star = '';
        if ~isnan(pFDR) && pFDR < 0.05, star = '+';
        elseif ~isnan(pP) && pP < 0.05, star = '*'; end

        text(ax, 0.03, 0.97, { ...
            sprintf('r_p=%+.2f (p=%.3f%s)', rP, pP, star), ...
            sprintf('r_{Sp}=%+.2f', rSp), ...
            sprintf('\\beta_{rob}=%+.3f', b_rob(2)), ...
            sprintf('n=%d', n)}, ...
            'Units','normalized','FontSize',5.5,'VerticalAlignment','top', ...
            'Color',COL.txt, 'BackgroundColor',[1 1 1 0.65], 'Margin', 1);

        % --- Labels (edges only) ---
        if ri == 1,    title(ax, cogLabels{di}, 'FontSize', 8, 'Color', COL.txt); end
        if ri == nInsp, xlabel(ax, sprintf('FC resid [%s]', fcSource{fi}), 'FontSize', 6); end
        if di == 1,    ylabel(ax, {fcShort{fi},'slope (resid)'}, ...
                              'FontSize', 6, 'FontWeight','bold'); end

        set(ax,'FontSize',6,'Box','off','TickDir','out');
    end
end

title(tI, sprintf(['Inspection: FC \\times Cognitive Slope (residualized | FD, TIV)   ' ...
    'blue--=OLS  red=robust  +=FDR<.05  *=p_{unc}<.05']), ...
    'FontSize', 9, 'Color', COL.txt);

exportFigure(figInsp, fullfile(outputDir, 'Inspection_FC_Slopes'));
fprintf('  -> Inspection_FC_Slopes saved.\n');

% ----- PATH A inspection: Apathy x FC for every component with supra -----
fprintf('\n  Building path-a inspection figure (Apathy x FC)...\n');
figPathA = figure('Units','centimeters','Position',[1 1 8*nSupraFC+2 9], ...
    'Color',COL.bg,'Name','Inspection: Apathy x FC (path a)');
set(figPathA,'DefaultAxesFontName','Helvetica','DefaultTextFontName','Helvetica');
tPA = tiledlayout(1, nSupraFC, 'TileSpacing','compact','Padding','compact');

for ci = 1:nSupraFC
    fi = supraIdx(ci);
    ax = nexttile; set(ax,'Color',COL.bg); hold on;

    x = fcMean.apathy_composite;
    y = fcSupra.(supraFCvars{ci});  % suprathresh — matches NBS edges
    Z  = [fcMean.mean_FD, fcMean.TIV];
    ok = ~isnan(x) & ~isnan(y) & all(~isnan(Z),2);
    if sum(ok) < 10
        text(ax, 0.5, 0.5, 'n<10', 'Units','normalized', ...
            'HorizontalAlignment','center','FontSize',8);
        continue;
    end

    Zo = Z(ok,:);
    xR = x(ok) - Zo*(Zo\x(ok));
    yR = y(ok) - Zo*(Zo\y(ok));

    densityScatter(ax, xR, yR);

    rP_p = partialcorr(x, y, Z, 'Type','Pearson','Rows','complete');
    rP_s = partialcorr(x, y, Z, 'Type','Spearman','Rows','complete');
    text(ax, 0.03, 0.97, { ...
        sprintf('r_{Pearson} = %+.3f', rP_p), ...
        sprintf('r_{Spearman} = %+.3f', rP_s), ...
        sprintf('n = %d', sum(ok))}, ...
        'Units','normalized','FontSize',7,'VerticalAlignment','top', ...
        'Color',COL.txt,'BackgroundColor',[1 1 1 0.75]);
    xlabel(ax, 'Apathy composite (resid)','FontSize',8);
    ylabel(ax, sprintf('%s FC (supra, resid)', fcShort{fi}), 'FontSize', 8);
    title(ax, sprintf('Path a: Apathy \\rightarrow %s', fcShort{fi}), ...
        'FontSize', 9, 'Color', COL.txt);
    formatAxis(ax, COL);
end
title(tPA, 'Path a: Apathy x FC (residualized | FD, TIV)   blue--=OLS  red=robust', ...
    'FontSize', 10, 'Color', COL.txt);
exportFigure(figPathA, fullfile(outputDir, 'Inspection_PathA_ApathyFC'));
fprintf('  -> Inspection_PathA_ApathyFC saved.\n');

%% ========================================================================
%  14. SENSITIVITY ANALYSIS — robustness of primary findings
%      Primary spec: Spearman partial, RAW slopes, FD+TIV
%      Four conditions per test:
%        (1) Baseline     : Pearson partial, RAW slopes,       FD+TIV
%        (2) Winsor       : Pearson partial, winsorised slopes, FD+TIV  (sensitivity)
%        (3) Spearman     : Spearman partial, RAW slopes,      FD+TIV  (== primary)
%        (4) ExtCovars    : Pearson partial, RAW slopes,
%                           FD + TIV + duration + LEDD + UPDRS-3
% =========================================================================
fprintf('\n[14] Sensitivity analysis...\n');

Z_base = [fcMean.mean_FD, fcMean.TIV];
Z_ext  = [fcMean.mean_FD, fcMean.TIV, ...
          baselineDur, baselineLEDD, baselineUPDRS3];

% NBS-matched covariates: age, sex, TIV, FD, QUIP, apathy×QUIP, depression.
% apathy×QUIP is supplied separately where apathy is the predictor so we can
% form the interaction per-row with the right apathy vector.
sexNum = sexFC;  % numeric sex (M=1, F=0); computed in Section 1
Z_nbs_noInt = [fcMean.mean_FD, fcMean.TIV, fcMean.age, sexNum, ...
               baselineQUIP, baselineDEP];

fprintf('    Z_ext n-complete    = %d\n', sum(all(~isnan(Z_ext),2)));
fprintf('    Z_base n-complete   = %d\n', sum(all(~isnan(Z_base),2)));
fprintf('    Z_nbs n-complete    = %d (pre-interaction; adds apathy x QUIP)\n', ...
    sum(all(~isnan(Z_nbs_noInt),2)));

condLabels = {'Baseline','Winsor','Spearman','ExtCovars'};
sensRows = {};

% ----- Path b: FC -> Global slope (apathy components) --------------------
fprintf('\n  Path b: FC -> Global slope\n');
fprintf('  %-28s  %-16s  %-16s  %-16s  %-16s\n', 'Test', ...
    condLabels{1}, condLabels{2}, condLabels{3}, condLabels{4});
fprintf('  %s\n', repmat('-',1,100));
for ci = 1:nSupraFC
    fi = supraIdx(ci);
    x = fcData(:, fi);
    y_raw  = slopeAll(:, 1);          % primary (raw, NoWinsor)
    y_wins = slopeAllWinsor(:, 1);    % sensitivity (winsorised at 1/99)
    [r1, p1] = partialcorr(x, y_raw,  Z_base, 'Type','Pearson','Rows','complete');
    n1 = sum(~isnan(x) & ~isnan(y_raw)  & all(~isnan(Z_base),2));    [r2, p2] = partialcorr(x, y_wins, Z_base, 'Type','Pearson','Rows','complete');
    n2 = sum(~isnan(x) & ~isnan(y_wins) & all(~isnan(Z_base),2));    [r3, p3] = partialcorr(x, y_raw,  Z_base, 'Type','Spearman','Rows','complete');
    n3 = sum(~isnan(x) & ~isnan(y_raw)  & all(~isnan(Z_base),2));    [r4, p4] = partialcorr(x, y_raw,  Z_ext,  'Type','Pearson','Rows','complete');
    n4 = sum(~isnan(x) & ~isnan(y_raw)  & all(~isnan(Z_ext),2));
    fprintf('  %-28s  %+.3f(p=%.3f)  %+.3f(p=%.3f)  %+.3f(p=%.3f)  %+.3f(p=%.3f)\n', ...
        fcShort{fi}, r1,p1, r2,p2, r3,p3, r4,p4);
    sensRows{end+1} = {'FCxGlobal', fcShort{fi}, ...
        r1,p1,n1, r2,p2,n2, r3,p3,n3, r4,p4,n4};
end

% ----- DAT -> FC (supraIdx) ---------------------------------------------
if hasDAT
    fprintf('\n  DAT -> FC (all supraIdx components)\n');
    fprintf('  %-28s  %-16s  %-16s  %-16s  %-16s\n', 'Test', ...
        condLabels{1}, condLabels{2}, condLabels{3}, condLabels{4});
    fprintf('  %s\n', repmat('-',1,100));
    for ci = 1:nSupraFC
        fi = supraIdx(ci);
        x = fcMean.DAT_mean;
        y = fcData(:, fi);
        % Extended covars for DAT->FC add age too (as in section 5b(iii))
        Zb_dat = [fcMean.mean_FD, fcMean.TIV, fcMean.age];
        Ze_dat = [fcMean.mean_FD, fcMean.TIV, fcMean.age, ...
                  baselineDur, baselineLEDD, baselineUPDRS3];
        [r1, p1] = partialcorr(x, y, Zb_dat, 'Type','Pearson','Rows','complete');
        n1 = sum(~isnan(x) & ~isnan(y) & all(~isnan(Zb_dat),2));        % Winsor column is same for DAT->FC (no slope on x or y)
        [r2,p2,n2] = deal(r1, p1, n1);
        [r3, p3] = partialcorr(x, y, Zb_dat, 'Type','Spearman','Rows','complete');
        n3 = sum(~isnan(x) & ~isnan(y) & all(~isnan(Zb_dat),2));        [r4, p4] = partialcorr(x, y, Ze_dat, 'Type','Pearson','Rows','complete');
        n4 = sum(~isnan(x) & ~isnan(y) & all(~isnan(Ze_dat),2));
        fprintf('  %-28s  %+.3f(p=%.3f)  %-16s  %+.3f(p=%.3f)  %+.3f(p=%.3f)\n', ...
            fcShort{fi}, r1,p1, '(n/a)', r3,p3, r4,p4);
        sensRows{end+1} = {'DATxFC', fcShort{fi}, ...
            r1,p1,n1, r2,p2,n2, r3,p3,n3, r4,p4,n4};
    end
end

% ----- Path a: Apathy -> FC (suprathresh, all supraIdx components) -------
fprintf('\n  Path a: Apathy -> FC (suprathresh, supraIdx)\n');
fprintf('  %-28s  %-16s  %-16s  %-16s  %-16s\n', 'Test', ...
    condLabels{1}, condLabels{2}, condLabels{3}, condLabels{4});
fprintf('  %s\n', repmat('-',1,100));
for ci = 1:nSupraFC
    fi = supraIdx(ci);
    x = fcMean.apathy_composite;
    y = fcSupra.(supraFCvars{ci});
    [r1, p1] = partialcorr(x, y, Z_base, 'Type','Pearson','Rows','complete');
    n1 = sum(~isnan(x) & ~isnan(y) & all(~isnan(Z_base),2));    [r2,p2,n2] = deal(r1, p1, n1);  % no slope here, no winsor to swap
    [r3, p3] = partialcorr(x, y, Z_base, 'Type','Spearman','Rows','complete');
    n3 = sum(~isnan(x) & ~isnan(y) & all(~isnan(Z_base),2));    [r4, p4] = partialcorr(x, y, Z_ext, 'Type','Pearson','Rows','complete');
    n4 = sum(~isnan(x) & ~isnan(y) & all(~isnan(Z_ext),2));
    fprintf('  %-28s  %+.3f(p=%.3f)  %-16s  %+.3f(p=%.3f)  %+.3f(p=%.3f)\n', ...
        fcShort{fi}, r1,p1, '(n/a)', r3,p3, r4,p4);
    sensRows{end+1} = {'ApathyToFC', fcShort{fi}, ...
        r1,p1,n1, r2,p2,n2, r3,p3,n3, r4,p4,n4};
end

% ----- Path a: NBS-matched covariate sensitivity -------------------------
%     Replicates NBS edge-level model at the aggregate level. Adds QUIP,
%     apathy x QUIP, depression, age, sex to FD+TIV.
%     Reports both Pearson and Spearman partial r.
fprintf('\n  Path a (NBS-matched covariates: age, sex, TIV, FD, QUIP, apathy*QUIP, dep):\n');
fprintf('  %-28s  %-18s  %-18s\n', 'Test', 'Pearson partial', 'Spearman partial');
fprintf('  %s\n', repmat('-',1,70));
nbsRows = {};
for ci = 1:nSupraFC
    fi = supraIdx(ci);
    x = fcMean.apathy_composite;
    y = fcSupra.(supraFCvars{ci});
    interAQ = fcMean.apathy_composite .* baselineQUIP;
    Z_nbs = [Z_nbs_noInt, interAQ];
    [rN, pN] = partialcorr(x, y, Z_nbs, 'Type','Pearson','Rows','complete');
    nN = sum(~isnan(x) & ~isnan(y) & all(~isnan(Z_nbs),2));    [rS, pS] = partialcorr(x, y, Z_nbs, 'Type','Spearman','Rows','complete');
    nS = sum(~isnan(x) & ~isnan(y) & all(~isnan(Z_nbs),2));
    fprintf('  %-28s  %+.3f (p=%.4f) %-4s  %+.3f (p=%.4f) %s   n=%d/%d\n', ...
        fcShort{fi}, rN, pN, pToStars(pN), rS, pS, pToStars(pS), nN, nS);
    nbsRows{end+1} = {fcShort{fi}, rN, pN, nN, rS, pS, nS};
end
if ~isempty(nbsRows)
    writetable(cell2table(vertcat(nbsRows{:}), ...
        'VariableNames', {'Component','Pearson_r','Pearson_p','Pearson_n', ...
                          'Spearman_r','Spearman_p','Spearman_n'}), ...
        fullfile(outputDir,'Table_Fig4_PathA_NBSmatched.csv'));
    fprintf('  -> Table_Fig4_PathA_NBSmatched.csv saved.\n');
end

% ----- Path b: NBS-matched and covariate-isolation tests -----------------
%     Same covariate sets as mediation sensitivity, applied to
%     FC → Global slope directly (to isolate which path loses signal).
fprintf('\n  Path b (FC -> Global slope) — covariate isolation:\n');
fprintf('  %-26s  %-10s  %-18s  %-18s\n', ...
    'Component', 'CovSet', 'Pearson partial', 'Spearman partial');
fprintf('  %s\n', repmat('-',1,78));
pathBIsolRows = {};
sexNumLocal = sexNum;
for ci = 1:nSupraFC
    fi = supraIdx(ci);
    x  = fcData(:, fi);
    y  = slopeAll(:, 1);
    interAQ_all = fcMean.apathy_composite .* baselineQUIP;
    covSetsDef = { ...
        'DepOnly',     [fcMean.mean_FD, fcMean.TIV, baselineDEP]; ...
        'QuipOnly',    [fcMean.mean_FD, fcMean.TIV, baselineQUIP, interAQ_all]; ...
        'AgeSexOnly',  [fcMean.mean_FD, fcMean.TIV, fcMean.age, sexNumLocal]; ...
        'NBSmatched',  [fcMean.mean_FD, fcMean.TIV, fcMean.age, sexNumLocal, ...
                        baselineQUIP, baselineDEP, interAQ_all]};
    for cs = 1:size(covSetsDef,1)
        csLbl = covSetsDef{cs,1};
        Z_cs  = covSetsDef{cs,2};
        [rP, pP] = partialcorr(x, y, Z_cs, 'Type','Pearson','Rows','complete');
        nP = sum(~isnan(x) & ~isnan(y) & all(~isnan(Z_cs),2));        [rS, pS] = partialcorr(x, y, Z_cs, 'Type','Spearman','Rows','complete');
        nS = sum(~isnan(x) & ~isnan(y) & all(~isnan(Z_cs),2));
        fprintf('  %-26s  %-10s  %+.3f (p=%.4f) %-4s  %+.3f (p=%.4f) %s\n', ...
            fcShort{fi}, csLbl, rP, pP, pToStars(pP), rS, pS, pToStars(pS));
        pathBIsolRows{end+1} = {fcShort{fi}, csLbl, rP, pP, nP, rS, pS, nS}; %#ok
    end
end
if ~isempty(pathBIsolRows)
    writetable(cell2table(vertcat(pathBIsolRows{:}), ...
        'VariableNames', {'Component','CovSet', ...
                          'Pearson_r','Pearson_p','Pearson_n', ...
                          'Spearman_r','Spearman_p','Spearman_n'}), ...
        fullfile(outputDir,'Table_Fig4_PathB_Isolation.csv'));
    fprintf('  -> Table_Fig4_PathB_Isolation.csv saved.\n');
end

% Save sensitivity table
sensT = cell2table(vertcat(sensRows{:}), 'VariableNames', { ...
    'Test','Component', ...
    'Baseline_r','Baseline_p','Baseline_n', ...
    'Winsor_r','Winsor_p','Winsor_n', ...
    'Spearman_r','Spearman_p','Spearman_n', ...
    'ExtCovars_r','ExtCovars_p','ExtCovars_n'});
writetable(sensT, fullfile(outputDir,'Table_Fig4_Sensitivity.csv'));
fprintf('\n  -> Table_Fig4_Sensitivity.csv saved.\n');

% ----- Mediation sensitivity --------------------------------------------
%     Re-run the §6 mediation for idxMediation across covariate sets and
%     transformations. The 'Baseline' row is intentionally identical to
%     §6's primary mediation: same outcome (slopeAll(:,1)), same
%     predictor (apathy_composite), same mediator (fcData(:,idxMediation)),
%     same covariates (FD + TIV), same okM filter (no per-outcome
%     baselineG — see feedback memory). Every other row perturbs one or
%     more of {slope transform, covariate set, rank vs linear} relative
%     to that baseline.
%
%     Each row reports BOTH percentile and BCa 95% CIs. BCa matches the
%     §6_BCa convention — standard for asymmetric sampling distributions
%     of the indirect effect (Efron 1987; MacKinnon et al. 2004).
fprintf('\n  Mediation sensitivity (idxMediation = fi=%d, %s):\n', ...
    idxMediation, labelMed);
fprintf('  %-22s  %-6s  %10s  %-23s  %-23s  %s\n', ...
    'Condition', 'N', 'ab', 'Percentile 95% CI', 'BCa 95% CI', 'prop');
fprintf('  %s\n', repmat('-',1,108));

% Covariate sets:
%   'base'      = FD + TIV (current primary)
%   'ext'       = FD + TIV + duration + LEDD + UPDRS-3 (disease severity —
%                 includes likely mediators; CAUTION: may over-control)
%   'nbs'       = FD + TIV + age + sex + QUIP + apathy*QUIP + depression
%                 (matches NBS covariates — confounders only, no mediators)
%   'depOnly'   = FD + TIV + depression  (isolate depression contribution)
%   'quipOnly'  = FD + TIV + QUIP + apathy*QUIP  (isolate QUIP contribution)
%   'ageSexOnly'= FD + TIV + age + sex  (isolate demographic contribution)
%
% wins = true uses the winsorised (1/99) slope copy (sensitivity);
% wins = false uses the raw slope (primary).
% useRank = true ranks all variables (apathy, FC, slope, covariates) before
% running the mediation — Spearman-equivalent linear mediation. This
% captures monotonic but non-linear relationships that OLS misses.
medSensCfg = { ...
    'Baseline',         false, 'base',       false; ...   % PRIMARY (raw + FD+TIV)
    'Winsor',           true,  'base',       false; ...   % sensitivity: winsor + FD+TIV
    'DepOnly',          false, 'depOnly',    false; ...
    'QuipOnly',         false, 'quipOnly',   false; ...
    'AgeSexOnly',       false, 'ageSexOnly', false; ...
    'NBSmatched',       false, 'nbs',        false; ...
    'ExtCovars',        false, 'ext',        false; ...
    'NBS+Winsor',       true,  'nbs',        false; ...
    'Ext+Winsor',       true,  'ext',        false; ...
    'Rank+Base',        false, 'base',       true;  ...
    'Rank+AgeSex',      false, 'ageSexOnly', true;  ...
    'Rank+NBS',         false, 'nbs',        true;  ...
    'Rank+Ext',         false, 'ext',        true};

medSensRows = {};
for ki = 1:size(medSensCfg,1)
    lbl     = medSensCfg{ki,1};
    wins    = medSensCfg{ki,2};
    cvSet   = medSensCfg{ki,3};
    useRank = medSensCfg{ki,4};

    mFC = fcData(:, idxMediation);
    mSl = ternary(wins, slopeAllWinsor(:,1), slopeAll(:,1));
    mAp = fcMean.apathy_composite;
    switch cvSet
        case 'base'
            mCov = [fcMean.mean_FD, fcMean.TIV];
        case 'ext'
            mCov = [fcMean.mean_FD, fcMean.TIV, ...
                    baselineDur, baselineLEDD, baselineUPDRS3];
        case 'nbs'
            interAQ = fcMean.apathy_composite .* baselineQUIP;
            mCov = [fcMean.mean_FD, fcMean.TIV, fcMean.age, sexNum, ...
                    baselineQUIP, baselineDEP, interAQ];
        case 'depOnly'
            mCov = [fcMean.mean_FD, fcMean.TIV, baselineDEP];
        case 'quipOnly'
            interAQ = fcMean.apathy_composite .* baselineQUIP;
            mCov = [fcMean.mean_FD, fcMean.TIV, baselineQUIP, interAQ];
        case 'ageSexOnly'
            mCov = [fcMean.mean_FD, fcMean.TIV, fcMean.age, sexNum];
        otherwise
            error('unknown cvSet');
    end

    % Rank-transform (Spearman-equivalent linear mediation)
    if useRank
        okR = ~isnan(mFC)&~isnan(mSl)&~isnan(mAp)&all(~isnan(mCov),2);
        mFC_r = nan(size(mFC));  mFC_r(okR) = tiedrank(mFC(okR));   mFC = mFC_r;
        mSl_r = nan(size(mSl));  mSl_r(okR) = tiedrank(mSl(okR));   mSl = mSl_r;
        mAp_r = nan(size(mAp));  mAp_r(okR) = tiedrank(mAp(okR));   mAp = mAp_r;
        for c = 1:size(mCov,2)
            col_r = nan(size(mCov,1),1);
            col_r(okR) = tiedrank(mCov(okR,c));
            mCov(:,c) = col_r;
        end
    end
    ok = ~isnan(mFC)&~isnan(mSl)&~isnan(mAp)&all(~isnan(mCov),2);
    n_ = sum(ok);
    if n_ < 20
        fprintf('  %-22s  n<20, skipped\n', lbl);
        continue;
    end
    fc_=mFC(ok); sl_=mSl(ok); ap_=mAp(ok); cov_=mCov(ok,:);

    % Observed a*b
    ba_obs = [ones(n_,1), ap_, cov_] \ fc_;  a_ = ba_obs(2);
    bb_obs = [ones(n_,1), fc_, ap_, cov_] \ sl_;  b_ = bb_obs(2);
    cp_ = bb_obs(3);
    bc_obs = [ones(n_,1), ap_, cov_] \ sl_;  c_ = bc_obs(2);
    ab_ = a_ * b_;
    prop_ = ternary(abs(c_)>1e-6, ab_/c_*100, NaN);

    % Bootstrap (percentile CI)
    ab_b = zeros(nBoot,1);
    for bi = 1:nBoot
        ix = randi(n_,n_,1);
        ba = [ones(n_,1), ap_(ix), cov_(ix,:)] \ fc_(ix);
        bb = [ones(n_,1), fc_(ix), ap_(ix), cov_(ix,:)] \ sl_(ix);
        ab_b(bi) = ba(2)*bb(2);
    end
    ci_ = prctile(ab_b, [2.5 97.5]);
    sig_ = (ci_(1) > 0) || (ci_(2) < 0);

    % Jackknife + BCa CI (same convention as §6_BCa)
    ab_j = zeros(n_,1);
    for ji = 1:n_
        keep = true(n_,1); keep(ji) = false;
        baJ = [ones(n_-1,1), ap_(keep), cov_(keep,:)] \ fc_(keep);
        bbJ = [ones(n_-1,1), fc_(keep), ap_(keep), cov_(keep,:)] \ sl_(keep);
        ab_j(ji) = baJ(2)*bbJ(2);
    end
    [ciB, z0_, ahat_] = bcaInterval(ab_, ab_b, ab_j, 0.05);
    sigB = (ciB(1) > 0) || (ciB(2) < 0);

    fprintf('  %-22s  %-6d  %+9.4f  [%+.4f, %+.4f] %-3s  [%+.4f, %+.4f] %-3s  %.1f%%\n', ...
        lbl, n_, ab_, ...
        ci_(1), ci_(2), ternary(sig_,'SIG','n.s.'), ...
        ciB(1), ciB(2), ternary(sigB,'SIG','n.s.'), prop_);
    medSensRows{end+1} = {lbl, n_, a_, b_, cp_, c_, ab_, ...
        ci_(1), ci_(2), double(sig_), ...
        ciB(1), ciB(2), double(sigB), z0_, ahat_, prop_}; %#ok
end
if ~isempty(medSensRows)
    writetable(cell2table(vertcat(medSensRows{:}), 'VariableNames', { ...
        'Condition','N','a','b','cp','c','ab', ...
        'Pct_CI_lo','Pct_CI_hi','Pct_Sig', ...
        'BCa_CI_lo','BCa_CI_hi','BCa_Sig','BCa_z0','BCa_ahat', ...
        'PctMediated'}), ...
        fullfile(outputDir,'Table_Fig4_Mediation_Sensitivity.csv'));
    fprintf('  -> Table_Fig4_Mediation_Sensitivity.csv saved.\n');
end

% ----- Age-stratified mediation ------------------------------------------
%     If mediation holds within younger and older subgroups separately,
%     age is just a marginal-distribution driver and the effect is genuine.
%     If it vanishes in both strata, age is driving everything.
fprintf('\n  Age-stratified mediation (median split):\n');
ageVec = fcMean.age;
ageMed = nanmedian(ageVec);
fprintf('  Age median split at %.1f years\n', ageMed);
fprintf('  %-22s  %-8s  %-30s  %s\n', ...
    'Stratum', 'N', 'ab [95%% CI]', 'mediated');
fprintf('  %s\n', repmat('-',1,85));

strata = {'Younger (age<=median)', ageVec <= ageMed; ...
          'Older  (age>median)', ageVec >  ageMed};
ageStratRows = {};
ageStratBoot = cell(size(strata,1), 1);  % store bootstrap distributions
ageStratN    = zeros(size(strata,1), 1);
for si = 1:size(strata,1)
    lbl = strata{si,1};
    sel = strata{si,2};
    mFC = fcData(:, idxMediation);
    mSl = slopeAll(:, 1);
    mAp = fcMean.apathy_composite;
    mCov = [fcMean.mean_FD, fcMean.TIV];
    ok  = sel & ~isnan(mFC)&~isnan(mSl)&~isnan(mAp)&all(~isnan(mCov),2);
    n_  = sum(ok);
    if n_ < 30
        fprintf('  %-22s  %-8d  n<30, skipped\n', lbl, n_);
        continue;
    end
    fc_=mFC(ok); sl_=mSl(ok); ap_=mAp(ok); cov_=mCov(ok,:);

    ba_obs = [ones(n_,1), ap_, cov_] \ fc_;     a_ = ba_obs(2);
    bb_obs = [ones(n_,1), fc_, ap_, cov_] \ sl_; b_ = bb_obs(2);
    cp_ = bb_obs(3);
    bc_obs = [ones(n_,1), ap_, cov_] \ sl_;     c_ = bc_obs(2);
    ab_ = a_ * b_;
    prop_ = ternary(abs(c_)>1e-6, ab_/c_*100, NaN);

    ab_b = zeros(nBoot,1);
    for bi = 1:nBoot
        ix = randi(n_,n_,1);
        ba = [ones(n_,1), ap_(ix), cov_(ix,:)] \ fc_(ix);
        bb = [ones(n_,1), fc_(ix), ap_(ix), cov_(ix,:)] \ sl_(ix);
        ab_b(bi) = ba(2)*bb(2);
    end
    ci_ = prctile(ab_b, [2.5 97.5]);
    sig_ = (ci_(1) > 0) || (ci_(2) < 0);
    fprintf('  %-22s  %-8d  %+.4f [%+.4f, %+.4f] %-3s  %.1f%%\n', ...
        lbl, n_, ab_, ci_(1), ci_(2), ternary(sig_,'SIG','n.s.'), prop_);
    ageStratRows{end+1} = {lbl, n_, a_, b_, cp_, c_, ab_, ci_(1), ci_(2), ...
        double(sig_), prop_}; %#ok
    ageStratBoot{si} = ab_b;
    ageStratN(si)    = n_;
end
if ~isempty(ageStratRows)
    writetable(cell2table(vertcat(ageStratRows{:}), 'VariableNames', { ...
        'Stratum','N','a','b','cp','c','ab','CI_lo','CI_hi','Sig','PctMediated'}), ...
        fullfile(outputDir,'Table_Fig4_Mediation_AgeStratified.csv'));
    fprintf('  -> Table_Fig4_Mediation_AgeStratified.csv saved.\n');
end

% ----- Moderated-mediation test (formal age x indirect-effect test) ------
%     H0: indirect effect is the same in both strata.
%     Bootstrap the difference: delta_ab = ab_older - ab_younger.
%     Since strata are independent, pair bootstrap iterations arbitrarily.
if ~isempty(ageStratBoot{1}) && ~isempty(ageStratBoot{2})
    fprintf('\n  Moderated-mediation test (bootstrap difference in indirect effect):\n');

    abY = ageStratBoot{1};       % younger
    abO = ageStratBoot{2};       % older
    nY  = ageStratN(1);
    nO  = ageStratN(2);
    delta = abO - abY;           % older minus younger
    deltaCI = prctile(delta, [2.5, 97.5]);
    deltaSig = (deltaCI(1) > 0) || (deltaCI(2) < 0);
    % One-sided p-value (proportion of bootstraps where direction flips)
    % If observed is negative (older has stronger mediation), count
    % bootstraps where older >= younger.
    meanDelta = mean(delta);
    pEmpirical = mean((delta >= 0) == (meanDelta < 0));

    fprintf('  %-28s  %-10s  %-22s  %-10s\n', ...
        'Contrast', 'n(Y|O)', 'delta_ab [95%% CI]', 'p_boot');
    fprintf('  %s\n', repmat('-',1,78));
    fprintf('  %-28s  %-10s  %+.4f [%+.4f, %+.4f] %-3s  %.4f\n', ...
        'Older - Younger', sprintf('%d|%d', nY, nO), meanDelta, ...
        deltaCI(1), deltaCI(2), ...
        ternary(deltaSig,'SIG','n.s.'), pEmpirical);

    writetable(table({'Older - Younger'}, [meanDelta], [deltaCI(1)], ...
        [deltaCI(2)], [double(deltaSig)], [pEmpirical], ...
        'VariableNames',{'Contrast','Delta_ab','CI_lo','CI_hi','Sig','p_boot'}), ...
        fullfile(outputDir,'Table_Fig4_ModeratedMediation.csv'));
    fprintf('  -> Table_Fig4_ModeratedMediation.csv saved.\n');

    % --- Age × apathy × FC interaction as a complementary test -----------
    %     Fit a single model on the full sample with apathy x age_group
    %     interacting in the b path. If the interaction term is significant,
    %     this corroborates the bootstrap difference.
    fprintf('\n  Complementary: apathy x age-group interaction on b path\n');
    ageGrp = double(ageVec > ageMed);    % 1 = older, 0 = younger
    mFC = fcData(:, idxMediation);
    mSl = slopeAll(:,1);
    mAp = fcMean.apathy_composite;
    okJ = ~isnan(mFC)&~isnan(mSl)&~isnan(mAp)&~isnan(fcMean.mean_FD)& ...
          ~isnan(fcMean.TIV)&~isnan(ageVec);
    % Slope ~ FC + FC*ageGrp + apathy + ageGrp + FD + TIV
    X = [mFC, mFC.*ageGrp, mAp, ageGrp, fcMean.mean_FD, fcMean.TIV];
    mdlInt = fitlm(X(okJ,:), mSl(okJ), ...
        'VarNames', {'FC','FC_x_Age','Apathy','AgeGrp','FD','TIV','Slope'});
    fprintf('  FC × AgeGrp interaction on slope:\n');
    disp(mdlInt.Coefficients);
end

% ----- (D1) Path a by age stratum ----------------------------------------
%     Does the apathy → FC association itself differ by age?
%     If path a is stable across strata, moderation is located in path b;
%     if path a is weaker in younger, the cascade begins at path a.
fprintf('\n  Path a (Apathy -> FC) by age stratum (FD+TIV adjusted):\n');
fprintf('  %-26s  %-10s  %-16s  %-16s\n', ...
    'Component', 'Stratum', 'Pearson partial', 'Spearman partial');
fprintf('  %s\n', repmat('-',1,76));
pathA_age_rows = {};
for ci = 1:nSupraFC
    fi = supraIdx(ci);
    x_all = fcMean.apathy_composite;
    y_all = fcSupra.(supraFCvars{ci});
    Z_all = [fcMean.mean_FD, fcMean.TIV];
    for si = 1:size(strata,1)
        sel = strata{si,2};
        xS = x_all(sel); yS = y_all(sel); ZS = Z_all(sel,:);
        [rP, pP] = partialcorr(xS, yS, ZS, 'Type','Pearson','Rows','complete');
        nP = sum(~isnan(xS) & ~isnan(yS) & all(~isnan(ZS),2));        [rS, pS] = partialcorr(xS, yS, ZS, 'Type','Spearman','Rows','complete');
        nS = sum(~isnan(xS) & ~isnan(yS) & all(~isnan(ZS),2));
        stratumShort = strata{si,1}(1:min(8,end));
        fprintf('  %-26s  %-10s  %+.3f (p=%.4f) %-4s  %+.3f (p=%.4f) %s\n', ...
            fcShort{fi}, stratumShort, rP, pP, pToStars(pP), rS, pS, pToStars(pS));
        pathA_age_rows{end+1} = {fcShort{fi}, strata{si,1}, rP, pP, nP, rS, pS, nS}; %#ok
    end
end
if ~isempty(pathA_age_rows)
    writetable(cell2table(vertcat(pathA_age_rows{:}), ...
        'VariableNames', {'Component','Stratum','Pearson_r','Pearson_p','Pearson_n', ...
                          'Spearman_r','Spearman_p','Spearman_n'}), ...
        fullfile(outputDir,'Table_Fig4_PathA_byAge.csv'));
    fprintf('  -> Table_Fig4_PathA_byAge.csv saved.\n');
end

% ----- (D2) Age-stratified mediation with NBS-matched covariates --------
%     Add depression, QUIP, apathy x QUIP, sex back in within each stratum
%     (omit age since we are already stratifying on age).
fprintf('\n  Age-stratified mediation with NBS-matched covariates:\n');
fprintf('  (covariates: FD, TIV, sex, QUIP, apathy*QUIP, depression)\n');
fprintf('  %-22s  %-8s  %-30s  %s\n', ...
    'Stratum', 'N', 'ab [95%% CI]', 'mediated');
fprintf('  %s\n', repmat('-',1,78));
ageStratNBS_rows = {};
ageStratNBS_boot = cell(size(strata,1),1);
for si = 1:size(strata,1)
    sel = strata{si,2};
    mFC = fcData(:, idxMediation);
    mSl = slopeAll(:,1);
    mAp = fcMean.apathy_composite;
    interAQ = fcMean.apathy_composite .* baselineQUIP;
    mCov = [fcMean.mean_FD, fcMean.TIV, sexNum, baselineQUIP, baselineDEP, interAQ];
    ok = sel & ~isnan(mFC)&~isnan(mSl)&~isnan(mAp)&all(~isnan(mCov),2);
    n_ = sum(ok);
    if n_ < 30
        fprintf('  %-22s  %-8d  n<30, skipped\n', strata{si,1}, n_);
        continue;
    end
    fc_=mFC(ok); sl_=mSl(ok); ap_=mAp(ok); cov_=mCov(ok,:);
    ba = [ones(n_,1), ap_, cov_] \ fc_;   a_ = ba(2);
    bb = [ones(n_,1), fc_, ap_, cov_] \ sl_;  b_ = bb(2);
    ab_ = a_ * b_;
    bc = [ones(n_,1), ap_, cov_] \ sl_;  c_ = bc(2);
    prop_ = ternary(abs(c_)>1e-6, ab_/c_*100, NaN);

    ab_b = zeros(nBoot,1);
    for bi = 1:nBoot
        ix = randi(n_,n_,1);
        baB = [ones(n_,1), ap_(ix), cov_(ix,:)] \ fc_(ix);
        bbB = [ones(n_,1), fc_(ix), ap_(ix), cov_(ix,:)] \ sl_(ix);
        ab_b(bi) = baB(2)*bbB(2);
    end
    ci_ = prctile(ab_b, [2.5 97.5]);
    sig_ = (ci_(1) > 0) || (ci_(2) < 0);
    fprintf('  %-22s  %-8d  %+.4f [%+.4f, %+.4f] %-3s  %.1f%%\n', ...
        strata{si,1}, n_, ab_, ci_(1), ci_(2), ternary(sig_,'SIG','n.s.'), prop_);
    ageStratNBS_rows{end+1} = {strata{si,1}, n_, ab_, ci_(1), ci_(2), ...
        double(sig_), prop_}; %#ok
    ageStratNBS_boot{si} = ab_b;
end

% Bootstrap difference of NBS-adjusted stratified mediations
if ~isempty(ageStratNBS_boot{1}) && ~isempty(ageStratNBS_boot{2})
    deltaN  = ageStratNBS_boot{2} - ageStratNBS_boot{1};
    deltaNc = prctile(deltaN, [2.5, 97.5]);
    deltaNs = (deltaNc(1) > 0) || (deltaNc(2) < 0);
    fprintf('  Moderated-mediation Δ (older - younger), NBS-adjusted:\n');
    fprintf('    %+.4f [%+.4f, %+.4f]  %s\n', mean(deltaN), deltaNc(1), deltaNc(2), ...
        ternary(deltaNs, 'SIG', 'n.s.'));
end
if ~isempty(ageStratNBS_rows)
    writetable(cell2table(vertcat(ageStratNBS_rows{:}), ...
        'VariableNames', {'Stratum','N','ab','CI_lo','CI_hi','Sig','PctMediated'}), ...
        fullfile(outputDir,'Table_Fig4_MediationByAge_NBSadjusted.csv'));
    fprintf('  -> Table_Fig4_MediationByAge_NBSadjusted.csv saved.\n');
end

% ----- (D2b) Age-stratified mediation with NBS + disease duration -------
%     Tests whether the age moderation is actually a "duration moderation
%     in disguise" (since age and duration correlate in PD).
fprintf('\n  Age-stratified mediation with NBS + duration covariates:\n');
fprintf('  (adds baseline disease duration to NBS covariate set)\n');
fprintf('  %-22s  %-8s  %-30s  %s\n', 'Stratum', 'N', 'ab [95%% CI]', 'mediated');
fprintf('  %s\n', repmat('-',1,78));
ageStratDur_boot = cell(size(strata,1),1);
for si = 1:size(strata,1)
    sel = strata{si,2};
    mFC = fcData(:, idxMediation);
    mSl = slopeAll(:,1);
    mAp = fcMean.apathy_composite;
    interAQ = fcMean.apathy_composite .* baselineQUIP;
    mCov = [fcMean.mean_FD, fcMean.TIV, sexNum, baselineQUIP, baselineDEP, ...
            interAQ, baselineDur];
    ok = sel & ~isnan(mFC)&~isnan(mSl)&~isnan(mAp)&all(~isnan(mCov),2);
    n_ = sum(ok);
    if n_ < 30
        fprintf('  %-22s  %-8d  n<30, skipped\n', strata{si,1}, n_);
        continue;
    end
    fc_=mFC(ok); sl_=mSl(ok); ap_=mAp(ok); cov_=mCov(ok,:);
    ba = [ones(n_,1), ap_, cov_] \ fc_;   a_ = ba(2);
    bb = [ones(n_,1), fc_, ap_, cov_] \ sl_;  b_ = bb(2);
    ab_ = a_ * b_;
    bc = [ones(n_,1), ap_, cov_] \ sl_;  c_ = bc(2);
    prop_ = ternary(abs(c_)>1e-6, ab_/c_*100, NaN);

    ab_b = zeros(nBoot,1);
    for bi = 1:nBoot
        ix = randi(n_,n_,1);
        baB = [ones(n_,1), ap_(ix), cov_(ix,:)] \ fc_(ix);
        bbB = [ones(n_,1), fc_(ix), ap_(ix), cov_(ix,:)] \ sl_(ix);
        ab_b(bi) = baB(2)*bbB(2);
    end
    ci_ = prctile(ab_b, [2.5 97.5]);
    sig_ = (ci_(1) > 0) || (ci_(2) < 0);
    fprintf('  %-22s  %-8d  %+.4f [%+.4f, %+.4f] %-3s  %.1f%%\n', ...
        strata{si,1}, n_, ab_, ci_(1), ci_(2), ternary(sig_,'SIG','n.s.'), prop_);
    ageStratDur_boot{si} = ab_b;
end
if ~isempty(ageStratDur_boot{1}) && ~isempty(ageStratDur_boot{2})
    deltaD  = ageStratDur_boot{2} - ageStratDur_boot{1};
    deltaDc = prctile(deltaD, [2.5, 97.5]);
    deltaDs = (deltaDc(1) > 0) || (deltaDc(2) < 0);
    fprintf('  Moderated-mediation Δ (older - younger), NBS + duration:\n');
    fprintf('    %+.4f [%+.4f, %+.4f]  %s\n', mean(deltaD), deltaDc(1), deltaDc(2), ...
        ternary(deltaDs, 'SIG', 'n.s.'));
end

% ----- (D3) Alternative age cutoffs -------------------------------------
%     Test whether the moderation pattern (older mediation SIG, younger n.s.)
%     is reproduced at multiple age thresholds. Reports older-stratum
%     mediation at cutoffs 55, 60, 65, 70, plus the median split already done.
fprintf('\n  Mediation at alternative age cutoffs (older stratum only):\n');
fprintf('  %-10s  %-10s  %-30s  %s\n', 'Cutoff (y)', 'N(older)', 'ab [95%% CI]', 'mediated');
fprintf('  %s\n', repmat('-',1,75));
altCutoffs = [55, 60, 65, 70];
altCut_rows = {};
for cut = altCutoffs
    sel = fcMean.age > cut;
    mFC = fcData(:, idxMediation);
    mSl = slopeAll(:,1);
    mAp = fcMean.apathy_composite;
    mCov = [fcMean.mean_FD, fcMean.TIV];
    ok = sel & ~isnan(mFC)&~isnan(mSl)&~isnan(mAp)&all(~isnan(mCov),2);
    n_ = sum(ok);
    if n_ < 30
        fprintf('  %-10.0f  %-10d  n<30, skipped\n', cut, n_);
        continue;
    end
    fc_=mFC(ok); sl_=mSl(ok); ap_=mAp(ok); cov_=mCov(ok,:);
    ba = [ones(n_,1), ap_, cov_] \ fc_;  a_ = ba(2);
    bb = [ones(n_,1), fc_, ap_, cov_] \ sl_;  b_ = bb(2);
    ab_ = a_ * b_;
    bc = [ones(n_,1), ap_, cov_] \ sl_;  c_ = bc(2);
    prop_ = ternary(abs(c_)>1e-6, ab_/c_*100, NaN);

    ab_b = zeros(nBoot,1);
    for bi = 1:nBoot
        ix = randi(n_,n_,1);
        baB = [ones(n_,1), ap_(ix), cov_(ix,:)] \ fc_(ix);
        bbB = [ones(n_,1), fc_(ix), ap_(ix), cov_(ix,:)] \ sl_(ix);
        ab_b(bi) = baB(2)*bbB(2);
    end
    ci_ = prctile(ab_b, [2.5 97.5]);
    sig_ = (ci_(1) > 0) || (ci_(2) < 0);
    fprintf('  %-10.0f  %-10d  %+.4f [%+.4f, %+.4f] %-3s  %.1f%%\n', ...
        cut, n_, ab_, ci_(1), ci_(2), ternary(sig_,'SIG','n.s.'), prop_);
    altCut_rows{end+1} = {cut, n_, ab_, ci_(1), ci_(2), double(sig_), prop_}; %#ok
end
if ~isempty(altCut_rows)
    writetable(cell2table(vertcat(altCut_rows{:}), ...
        'VariableNames', {'Cutoff_y','N_older','ab','CI_lo','CI_hi','Sig','PctMediated'}), ...
        fullfile(outputDir,'Table_Fig4_AltCutoffs.csv'));
    fprintf('  -> Table_Fig4_AltCutoffs.csv saved.\n');
end

%% ========================================================================
%  LOCAL FUNCTIONS
% =========================================================================

function yS = gaussSmooth(x, y, sigma)
    yS = nan(size(y));
    for i = 1:numel(x)
        w = exp(-0.5 * ((x - x(i)) / sigma).^2);
        w = w / sum(w);
        yS(i) = sum(w .* y);
    end
end

function [r,p,rP,pP,n,nP,x,y,mdl] = fcSlopeCorr(fcV,sV,meanFD,TIV)
    ok=~isnan(fcV)&~isnan(sV); x=fcV(ok);y=sV(ok);n=sum(ok);
    r=NaN;p=NaN;mdl=[];
    if n>=5,[r,p]=corr(x,y,'Type','Pearson'); mdl=fitlm(x,y); end
    okP=ok&~isnan(meanFD)&~isnan(TIV); nP=sum(okP); rP=NaN;pP=NaN;
    if nP>=10
        rP=partialcorr(fcV(okP),sV(okP),[meanFD(okP),TIV(okP)],'Rows','complete');
        df=nP-4; tS=rP*sqrt(df/(1-rP^2)); pP=2*(1-tcdf(abs(tS),df));
    end
end

function [slopes,nM] = extractBLUPsFromModel(mdl,fcPATNO,nFC)
    slopes=nan(nFC,1); nM=0;
    [reV,reN]=randomEffects(mdl);
    reT=unique(reN.Name); tRE='';
    for tn=reT', if ~strcmp(tn{1},'(Intercept)'), tRE=tn{1}; break; end, end
    if isempty(tRE), return; end
    isS=strcmp(reN.Name,tRE); sRE=reV(isS); lvl=reN.Level(isS);
    if iscategorical(lvl), rP=double(string(lvl));
    elseif iscell(lvl)||isstring(lvl), rP=str2double(string(lvl));
    else, rP=double(lvl); end
    feS=mdl.Coefficients.Estimate(strcmp(mdl.Coefficients.Name,tRE));
    for si=1:nFC
        row=find(rP==fcPATNO(si),1);
        if ~isempty(row), slopes(si)=feS+sRE(row); nM=nM+1; end
    end
end

function pFDR=fdr_bh(pV)
    pV=pV(:);m=numel(pV);[ps,si]=sort(pV);pa=nan(m,1);pa(m)=ps(m);
    for i=m-1:-1:1, pa(i)=min(pa(i+1),ps(i)*m/i); end
    pa=min(pa,1);pFDR=nan(m,1);pFDR(si)=pa;
end

function [ci, z0, ahat] = bcaInterval(theta_hat, theta_boot, theta_jack, alpha)
% Bias-corrected & accelerated bootstrap CI (Efron 1987; DiCiccio & Efron 1996).
%   theta_hat  : observed point estimate (scalar)
%   theta_boot : bootstrap replicates of theta (B x 1)
%   theta_jack : jackknife replicates of theta, leave-one-out (n x 1)
%   alpha      : two-sided level (e.g., 0.05 for 95% CI)
    theta_boot = theta_boot(~isnan(theta_boot));
    theta_jack = theta_jack(~isnan(theta_jack));
    B = numel(theta_boot);
    if B < 100 || numel(theta_jack) < 3
        ci = [NaN NaN]; z0 = NaN; ahat = NaN; return;
    end
    % Bias correction: z0 = Phi^-1( P(theta_boot < theta_hat) ).
    prop_lt = mean(theta_boot < theta_hat) + 0.5*mean(theta_boot == theta_hat);
    prop_lt = min(max(prop_lt, 1/B), 1 - 1/B);  % clamp to avoid +/- Inf
    z0 = norminv(prop_lt);
    % Acceleration from jackknife skewness.
    jmean = mean(theta_jack);
    num = sum((jmean - theta_jack).^3);
    den = 6 * (sum((jmean - theta_jack).^2))^1.5;
    if den == 0, ahat = 0; else, ahat = num / den; end
    % Adjusted quantiles.
    zL = norminv(alpha/2);
    zU = norminv(1 - alpha/2);
    aL = normcdf(z0 + (z0 + zL)/(1 - ahat*(z0 + zL)));
    aU = normcdf(z0 + (z0 + zU)/(1 - ahat*(z0 + zU)));
    ci = prctile(theta_boot, 100*[aL, aU]);
end

function v=safeNumeric(v_in)
    if isnumeric(v_in),v=double(v_in);
    elseif iscell(v_in)||isstring(v_in),v=str2double(string(v_in));
    elseif iscategorical(v_in),v=double(v_in);
    elseif islogical(v_in),v=double(v_in);
    else,v=nan(numel(v_in),1);end
end

function v=safeCol(T,name)
    if ismember(name,T.Properties.VariableNames),v=safeNumeric(T.(name));
    else,v=nan(height(T),1);end
end

function s=pToStars(p)
    if isnan(p),s='';elseif p<0.001,s='***';elseif p<0.01,s='**';elseif p<0.05,s='*';else,s='';end
end

function out=ternary(c,a,b), if c,out=a;else,out=b;end, end

function s=fmtP(p)
    % Format p-value: show '< 0.001' when p rounds to 0.000
    if isnan(p), s='NaN';
    elseif p < 0.001, s='< 0.001';
    else, s=sprintf('%.3f', p); end
end

function formatAxis(ax,COL)
    set(ax,'FontSize',8,'Box','off','TickDir','out','XGrid','off','YGrid','off',...
        'GridColor',COL.grid,'GridAlpha',0.4,'FontWeight','bold');
    try
        if ~isempty(ax.XLabel), set(ax.XLabel, 'FontWeight', 'bold'); end
        if ~isempty(ax.YLabel), set(ax.YLabel, 'FontWeight', 'bold'); end
        if ~isempty(ax.Title),  set(ax.Title,  'FontWeight', 'bold'); end
    catch
    end
end

function exportFigure(fig,fp)
    set(fig,'PaperPositionMode','auto','Renderer','painters');
    print(fig,[fp '.pdf'],'-dpdf','-bestfit');
    print(fig,[fp '.png'],'-dpng','-r600');
    close(fig);
end

function densityScatter(ax, x, y, varargin)
    % DENSITYSCATTER  Density-coloured scatter with OLS + non-parametric smoother.
    %   Colours points by local 2D KDE density.
    %   Overlays an OLS line (dashed) and a LOWESS smoother (solid) —
    %   the LOWESS curve visualises the monotonic trend consistent with
    %   a Spearman partial correlation, and does not assume linearity.
    %
    %   Options (name-value):
    %     MarkerSize  (default 24)
    %     OLSColor    (default [0.30 0.30 0.30])
    %     RobustColor (default [0.85 0.20 0.20])         — smoother colour
    %     Cmap        (default densityColormap())
    %     Alpha       (default 0.75)
    %     ShowOLS     (default true)                      — dashed linear line
    %     ShowRobust  (default true)                      — LOWESS smoother
    %     SmootherType (default 'lowess')                 — 'lowess' | 'loess' |
    %                                                       'rlowess' | 'robust'
    %                  ('robust' reverts to linear bisquare fit)
    %     SmoothSpan  (default 0.33 of N, min 15)
    p = inputParser;
    p.addParameter('MarkerSize', 24);
    p.addParameter('OLSColor',    [0.30 0.30 0.30]);
    p.addParameter('RobustColor', [0.85 0.20 0.20]);
    p.addParameter('Cmap',        densityColormap());
    p.addParameter('Alpha',       0.75);
    p.addParameter('ShowOLS',     true);
    p.addParameter('ShowRobust',  true);
    p.addParameter('SmootherType','rloess');  % robust local quadratic — default
    p.addParameter('SmoothSpan',  []);          % fraction of N (0-1); default adaptive
    p.parse(varargin{:});
    opt = p.Results;

    ok = ~isnan(x) & ~isnan(y);
    x = x(ok); y = y(ok);
    n = numel(x);
    if n < 5, return; end

    % 2D density at each point (Gaussian KDE; fall back to bin density)
    dens = nan(n,1);
    try
        dens = mvksdensity([x y], [x y], ...
            'Bandwidth', [std(x)*0.35, std(y)*0.35]);
    catch
        try
            [N,xE,yE] = histcounts2(x, y, 25, 25);
            xb = max(min(discretize(x,xE), size(N,1)), 1);
            yb = max(min(discretize(y,yE), size(N,2)), 1);
            for i=1:n, dens(i) = N(xb(i),yb(i)); end
        catch
            dens = ones(n,1);
        end
    end

    hold(ax, 'on');
    scatter(ax, x, y, opt.MarkerSize, dens, 'filled', ...
        'MarkerFaceAlpha', opt.Alpha, 'MarkerEdgeColor', 'none');
    colormap(ax, opt.Cmap);

    xF = linspace(min(x), max(x), 200);
    if opt.ShowOLS
        b_ols = [ones(n,1), x] \ y;
        plot(ax, xF, b_ols(1)+b_ols(2)*xF, '--', ...
            'Color', opt.OLSColor, 'LineWidth', 1.2);
    end
    if opt.ShowRobust
        if strcmpi(opt.SmootherType, 'robust')
            % Legacy behaviour: linear bisquare robustfit
            try
                b_rob = robustfit(x, y);
                plot(ax, xF, b_rob(1)+b_rob(2)*xF, '-', ...
                    'Color', opt.RobustColor, 'LineWidth', 1.8);
            catch
            end
        else
            % Non-parametric smoother (default rLOESS = robust local quadratic)
            %   Span adapts to sample size: smaller for large n (captures trend),
            %   larger for small n (stability). 'rloess' down-weights outliers.
            if ~isempty(opt.SmoothSpan) && opt.SmoothSpan < 1
                winSize = max(round(opt.SmoothSpan * n), 15);
            elseif ~isempty(opt.SmoothSpan)
                winSize = opt.SmoothSpan;
            else
                % Adaptive: 20-40% of N, floor 20
                if n < 50,        winSize = max(round(0.45*n), 15);
                elseif n < 150,   winSize = round(0.30*n);
                else              winSize = round(0.22*n);
                end
            end
            [xs, si] = sort(x);
            ys = y(si);
            ok_fit = true;
            try
                ySmooth = smoothdata(ys, opt.SmootherType, winSize);
            catch
                try
                    % Fallback: built-in smooth() from Curve Fitting Toolbox
                    ySmooth = smooth(xs, ys, winSize/n, 'rloess'); %#ok
                catch
                    % Final fallback: linear robustfit so the panel isn't
                    % left without any reference line.
                    try
                        b_rob = robustfit(x, y);
                        plot(ax, xF, b_rob(1)+b_rob(2)*xF, '-', ...
                            'Color', opt.RobustColor, 'LineWidth', 1.8);
                    catch
                    end
                    ok_fit = false;
                end
            end
            if ok_fit
                plot(ax, xs, ySmooth, '-', 'Color', opt.RobustColor, ...
                    'LineWidth', 1.8);
            end
        end
    end
end

function cm = densityColormap()
    % Sequential colormap for density plots: navy → indigo → magenta → orange → gold
    key = [ ...
        0.05 0.05 0.15;
        0.22 0.12 0.45;
        0.55 0.15 0.50;
        0.90 0.35 0.20;
        1.00 0.82 0.35];
    n   = 256;
    pos = linspace(0, 1, size(key,1));
    t   = linspace(0, 1, n)';
    cm  = zeros(n, 3);
    for ci = 1:3
        cm(:,ci) = interp1(pos, key(:,ci), t, 'pchip');
    end
end

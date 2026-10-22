% Role: Connectivity-to-outcome chain, step 2 of 3 for Figure 4. Spearman partial correlations of
%     each connectivity component with 42 bilaterally averaged FreeSurfer volumes (9 subcortical,
%     33 cortical), adjusted for age, sex, mean FD and eTIV (M1) and in addition for education,
%     disease duration, LEDD and MDS-UPDRS Part III (M2), with BH FDR across regions within
%     component. Then the 24-month hippocampal atrophy rate (bilateral, left, right; volumes
%     eTIV-corrected at both time points with a healthy-control normative model) against baseline
%     connectivity, adjusted for age, sex, mean FD and baseline volume (M1) and in addition for
%     education, disease duration, LEDD and MDS-UPDRS Part III (M2), with BH FDR across components
%     within side. Writes Table_Structural_Correlates.csv, Table_Hippocampus_Slope_x_FC.csv and
%     Table_Hippocampus_Longitudinal_Data.csv to outputDir.
% Run: after fc_cognitive_slope_mediation.m, on the workspace it leaves.
% Manuscript item: Figure 4b,c values (drawn by figure4_connectivity_outcomes.m); Results
%     (structural correlates, hippocampal atrophy); the 42-region supplementary table.
fprintf('\n[Structural] FreeSurfer volumes x FC components...\n');

%% 1. Load structural measures
structFile = fullfile(getenv('PPMI_CLINICAL_ROOT'), 'MRI_Data/all_structural_measures.csv');
structT = readtable(structFile);
nStruct = height(structT);
structPATNO = double(structT.PATNO);
fprintf('  Structural file: %d subjects\n', nStruct);

% eTIV
structeTIV = double(structT.fs_EstimatedTotalIntraCranialVol);

%% 2. HC normative eTIV model for subcortical volumes
%    HC subjects available in FS7_ASEG download (n=201).
fsDir = fullfile(getenv('PPMI_CLINICAL_ROOT'), 'ppmi_data/download');
asegT = readtable(fullfile(fsDir,'FS7_ASEG_VOL_16Mar2026.csv'));

% Identify HC via curated subgroup
if ~exist('Tcurated','var') || ~istable(Tcurated)
    curFile = fullfile(getenv('PPMI_CLINICAL_ROOT'), 'ppmi_data/PPMI_Curated_Data_Cut_Public_20260223.xlsx');
    Tcurated = readtable(curFile, 'Sheet','20260105');
end
[~,uIdx] = unique(Tcurated.PATNO,'first');
curBL = Tcurated(uIdx,:);
if iscell(curBL.PATNO)||isstring(curBL.PATNO), curPATNO = str2double(string(curBL.PATNO));
else, curPATNO = double(curBL.PATNO); end
curSubgr = string(curBL.subgroup);

asegPATNO = double(asegT.PATNO);
isHC = false(height(asegT),1);
for si=1:height(asegT)
    r = find(curPATNO == asegPATNO(si), 1);
    if ~isempty(r) && curSubgr(r) == "Healthy Control", isHC(si) = true; end
end
hc_eTIV = double(asegT.EstimatedTotalIntraCranialVol(isHC));
hc_mean_eTIV = nanmean(hc_eTIV);
fprintf('  HC for normative model: n=%d (mean eTIV = %.0f)\n', sum(isHC), hc_mean_eTIV);

%% 3. Define ROIs — subcortical volumes (bilateral avg)
subPairs = { ...
    'Hippocampus',  'fs_vol_Left_Hippocampus','fs_vol_Right_Hippocampus'; ...
    'Amygdala',     'fs_vol_Left_Amygdala','fs_vol_Right_Amygdala'; ...
    'Caudate',      'fs_vol_Left_Caudate','fs_vol_Right_Caudate'; ...
    'Putamen',      'fs_vol_Left_Putamen','fs_vol_Right_Putamen'; ...
    'Accumbens',    'fs_vol_Left_Accumbens_area','fs_vol_Right_Accumbens_area'; ...
    'Thalamus',     'fs_vol_Left_Thalamus','fs_vol_Right_Thalamus'; ...
    'Pallidum',     'fs_vol_Left_Pallidum','fs_vol_Right_Pallidum'; ...
    'Cereb Cortex', 'fs_vol_Left_Cerebellum_Cortex','fs_vol_Right_Cerebellum_Cortex'};
nSubPairs = size(subPairs,1);
nSub = nSubPairs + 1;  % + brainstem
subNames = [subPairs(:,1); {'Brain Stem'}];

% HC normative betas for subcortical ROIs
% Map structural-file column names to ASEG column names (strip 'fs_vol_')
hc_beta_sub = nan(nSub, 1);
for ri=1:nSubPairs
    asegL = strrep(subPairs{ri,2},'fs_vol_','');
    asegR = strrep(subPairs{ri,3},'fs_vol_','');
    hcV = (double(asegT.(asegL)(isHC)) + double(asegT.(asegR)(isHC))) / 2;
    ok = ~isnan(hcV) & ~isnan(hc_eTIV);
    if sum(ok)>15, b=polyfit(hc_eTIV(ok),hcV(ok),1); hc_beta_sub(ri)=b(1);
    else, hc_beta_sub(ri)=0; end
end
% Brainstem
hcBS = double(asegT.Brain_Stem(isHC));
ok = ~isnan(hcBS) & ~isnan(hc_eTIV);
if sum(ok)>15, b=polyfit(hc_eTIV(ok),hcBS(ok),1); hc_beta_sub(nSub)=b(1);
else, hc_beta_sub(nSub)=0; end

%% 4. Define ROIs — cortical gray matter volumes (bilateral avg)
%   Note: 'bankssts' (banks of superior temporal sulcus) is excluded from
%   the Desikan-Killiany list — it's defined by a sulcal bank rather than
%   cytoarchitecture, is anatomically heterogeneous, and is dropped from
%   many modern aparc analyses. 33 cortical regions retained.
aparcReg = {'caudalanteriorcingulate','caudalmiddlefrontal','cuneus', ...
    'entorhinal','fusiform','inferiorparietal','inferiortemporal','isthmuscingulate', ...
    'lateraloccipital','lateralorbitofrontal','lingual','medialorbitofrontal', ...
    'middletemporal','parahippocampal','paracentral','parsopercularis','parsorbitalis', ...
    'parstriangularis','pericalcarine','postcentral','posteriorcingulate','precentral', ...
    'precuneus','rostralanteriorcingulate','rostralmiddlefrontal','superiorfrontal', ...
    'superiorparietal','superiortemporal','supramarginal','frontalpole','temporalpole', ...
    'transversetemporal','insula'};
nCort = numel(aparcReg);
nROI  = nSub + nCort;
allNames = [subNames; aparcReg(:)];
fprintf('  ROIs: %d subcort + %d cortical = %d\n', nSub, nCort, nROI);

%% 5. Match to FC sample and build ROI matrix
nFC = numel(fcPATNO);
structIdx = nan(nFC,1);
for si=1:nFC
    r = find(structPATNO == fcPATNO(si), 1);
    if ~isempty(r), structIdx(si) = r; end
end
nMatched = sum(~isnan(structIdx));
fprintf('  FC subjects matched to structural: %d / %d\n', nMatched, nFC);

% Extract raw values
rawROI = nan(nFC, nROI);
fc_eTIV = nan(nFC, 1);
for si=1:nFC
    if isnan(structIdx(si)), continue; end
    idx = structIdx(si);
    fc_eTIV(si) = structeTIV(idx);
    % Subcortical bilateral avg
    for ri=1:nSubPairs
        L = double(structT.(subPairs{ri,2})(idx));
        R = double(structT.(subPairs{ri,3})(idx));
        rawROI(si,ri) = (L + R) / 2;
    end
    rawROI(si,nSub) = double(structT.fs_vol_Brain_Stem(idx));
    % Cortical bilateral avg (gray matter volume)
    for ri=1:nCort
        colL = ['fs_lh_' aparcReg{ri} '_volume'];
        colR = ['fs_rh_' aparcReg{ri} '_volume'];
        L = double(structT.(colL)(idx));
        R = double(structT.(colR)(idx));
        rawROI(si, nSub+ri) = (L + R) / 2;
    end
end

% HC normative correction of the subcortical volumes (used for the hippocampal slopes in
% section 13; the cross-sectional analysis uses rawROI with eTIV as a covariate).
normROI = rawROI;
for ri=1:nSub
    for si=1:nFC
        if ~isnan(rawROI(si,ri)) && ~isnan(fc_eTIV(si)) && ~isnan(hc_beta_sub(ri))
            normROI(si,ri) = rawROI(si,ri) - hc_beta_sub(ri) * (fc_eTIV(si) - hc_mean_eTIV);
        end
    end
end
fprintf('  Primary volumes: rawROI (uncorrected); eTIV partialled in covariate set.\n');

%% 6. Predictors + covariates
predData  = nan(nFC, nAllFC);
predNames = cell(1, nAllFC);
for fi=1:nAllFC
    predData(:,fi) = fcData(:,fi);   % suprathreshold (from fc_cognitive_slope_mediation.m)
    predNames{fi}  = fcShort{fi};
end
nPred = nAllFC;

% Covariates, two models:
%   M1 (primary): age + sex + FD + eTIV, the same for subcortical and cortical ROIs
%   M2 (sensitivity): M1 + education + disease duration + LEDD + UPDRS-III
Z_sub  = [fcMean.age, sexFC, fcMean.mean_FD, fc_eTIV];
Z_cort = [fcMean.age, sexFC, fcMean.mean_FD, fc_eTIV];

% Defensive existence checks for M2 sensitivity covariates (fc_cognitive_slope_mediation.m
% defines these at L423-428; this guards against running structural alone).
if ~exist('baselineEduc','var'),   baselineEduc   = nan(nFC,1); end
if ~exist('baselineDur','var'),    baselineDur    = nan(nFC,1); end
if ~exist('baselineLEDD','var'),   baselineLEDD   = nan(nFC,1); end
if ~exist('baselineUPDRS3','var'), baselineUPDRS3 = nan(nFC,1); end

Z_sub_M2  = [Z_sub,  baselineEduc, baselineDur, baselineLEDD, baselineUPDRS3];
Z_cort_M2 = [Z_cort, baselineEduc, baselineDur, baselineLEDD, baselineUPDRS3];

fprintf('  M1 Primary covariates (uniform): age, sex, FD, eTIV\n');
fprintf('  M2 Sensitivity adds: EDUC + duration + LEDD + UPDRS-3\n');
fprintf('    M2 n-complete (subcort Z): %d / %d\n', sum(all(~isnan(Z_sub_M2),2)), nFC);
fprintf('    M2 n-complete (cort  Z):   %d / %d\n', sum(all(~isnan(Z_cort_M2),2)), nFC);
fprintf('  Total tests: %d ROIs x %d FC x 2 models = %d\n', nROI, nPred, nROI*nPred*2);

%% 7. Spearman partial correlation — M1 (primary) and M2 (sensitivity)
%   Uses rawROI uniformly for both subcortical and cortical (eTIV is in Z).
rMat    = nan(nROI, nPred);  pMat    = nan(nROI, nPred);   % M1 PRIMARY
rMat_M2 = nan(nROI, nPred);  pMat_M2 = nan(nROI, nPred);   % M2 SENSITIVITY
for ri=1:nROI
    y = rawROI(:,ri);   % uncorrected; eTIV partialled in Zuse below
    if ri <= nSub
        Zuse_M1 = Z_sub;   Zuse_M2 = Z_sub_M2;
    else
        Zuse_M1 = Z_cort;  Zuse_M2 = Z_cort_M2;
    end
    for pi=1:nPred
        x = predData(:,pi);
        % M1 PRIMARY
        ok1 = ~isnan(x) & ~isnan(y) & all(~isnan(Zuse_M1),2);
        if sum(ok1) >= 20
            [rho1, pv1] = partialcorr(x(ok1), y(ok1), Zuse_M1(ok1,:), 'Type','Spearman');
            rMat(ri,pi) = rho1; pMat(ri,pi) = pv1;
        end
        % M2 SENSITIVITY
        ok2 = ~isnan(x) & ~isnan(y) & all(~isnan(Zuse_M2),2);
        if sum(ok2) >= 20
            [rho2, pv2] = partialcorr(x(ok2), y(ok2), Zuse_M2(ok2,:), 'Type','Spearman');
            rMat_M2(ri,pi) = rho2; pMat_M2(ri,pi) = pv2;
        end
    end
end

%% 8. Within-component FDR
%    BH FDR across the 42 ROIs for each FC component separately, within each model (M1, M2).
pFDR    = nan(nROI, nPred);   % M1
pFDR_M2 = nan(nROI, nPred);   % M2
for pi=1:nPred
    % --- M1 ---
    pVec = pMat(:,pi); v = ~isnan(pVec);
    if any(v)
        [sP,sI] = sort(pVec(v)); m=numel(sP);
        adj=nan(m,1); adj(m)=sP(m);
        for k=m-1:-1:1, adj(k)=min(adj(k+1), sP(k)*m/k); end
        adj=min(adj,1);
        tmp=nan(sum(v),1); tmp(sI)=adj;
        out=nan(nROI,1); out(v)=tmp;
        pFDR(:,pi) = out;
    end
    % --- M2 ---
    pVec2 = pMat_M2(:,pi); v2 = ~isnan(pVec2);
    if any(v2)
        [sP,sI] = sort(pVec2(v2)); m=numel(sP);
        adj=nan(m,1); adj(m)=sP(m);
        for k=m-1:-1:1, adj(k)=min(adj(k+1), sP(k)*m/k); end
        adj=min(adj,1);
        tmp=nan(sum(v2),1); tmp(sI)=adj;
        out=nan(nROI,1); out(v2)=tmp;
        pFDR_M2(:,pi) = out;
    end
end

piSMDA = 2;  % RH SM-DA index in predData
nTests_smda = sum(~isnan(pMat(:,piSMDA)));
nSigFDR = sum(pFDR(:,piSMDA)<0.05);
nSigUnc = sum(pMat(:,piSMDA)<0.01);
fprintf('  [M1 Primary] Within-component FDR (RH SM-DA, %d tests): FDR<.05 = %d  |  p_unc<.01 = %d\n', ...
    nTests_smda, nSigFDR, nSigUnc);

% Also report other components briefly (M1)
for pi=1:nPred
    n_fdr = sum(pFDR(:,pi)<0.05);
    n_unc = sum(pMat(:,pi)<0.01);
    fprintf('    %-25s  FDR<.05: %d  p<.01: %d\n', predNames{pi}, n_fdr, n_unc);
end

% M2 sensitivity summary
nTests_smda_M2 = sum(~isnan(pMat_M2(:,piSMDA)));
nSigFDR_M2 = sum(pFDR_M2(:,piSMDA)<0.05);
nSigUnc_M2 = sum(pMat_M2(:,piSMDA)<0.01);
fprintf('  [M2 Sensitivity] Within-component FDR (RH SM-DA, %d tests): FDR<.05 = %d  |  p_unc<.01 = %d\n', ...
    nTests_smda_M2, nSigFDR_M2, nSigUnc_M2);
for pi=1:nPred
    n_fdr2 = sum(pFDR_M2(:,pi)<0.05);
    n_unc2 = sum(pMat_M2(:,pi)<0.01);
    fprintf('    %-25s  FDR<.05: %d  p<.01: %d\n', predNames{pi}, n_fdr2, n_unc2);
end

%% 9. Print RH SM-DA results — M1 primary, then M2 sensitivity side-by-side
fprintf('\n  RH SM-DA structural correlates (M1 Primary | M2 Sensitivity, within-component FDR):\n');
fprintf('  %-28s  %7s  %10s  %10s   |  %7s  %10s  %10s\n', ...
    'ROI','M1 rho','M1 p_unc','M1 p_FDR','M2 rho','M2 p_unc','M2 p_FDR');
fprintf('  %s\n', repmat('-',1,98));
for ri=1:nROI
    isFDR = ~isnan(pFDR(ri,piSMDA)) && pFDR(ri,piSMDA)<0.05;
    isUnc = ~isnan(pMat(ri,piSMDA)) && pMat(ri,piSMDA)<0.01;
    isFDR2 = ~isnan(pFDR_M2(ri,piSMDA)) && pFDR_M2(ri,piSMDA)<0.05;
    isUnc2 = ~isnan(pMat_M2(ri,piSMDA)) && pMat_M2(ri,piSMDA)<0.01;
    if isFDR || isUnc || isFDR2 || isUnc2
        tag1 = ''; if isFDR, tag1='*FDR'; elseif isUnc, tag1='*unc'; end
        tag2 = ''; if isFDR2, tag2='*FDR'; elseif isUnc2, tag2='*unc'; end
        fprintf('  %-28s  %+.3f  %.4f  %.4f %-4s |  %+.3f  %.4f  %.4f %-4s\n', ...
            allNames{ri}, ...
            rMat(ri,piSMDA),    pMat(ri,piSMDA),    pFDR(ri,piSMDA),    tag1, ...
            rMat_M2(ri,piSMDA), pMat_M2(ri,piSMDA), pFDR_M2(ri,piSMDA), tag2);
    end
end

%% 11. Save table (all components, within-component FDR, both models)
%   Long format: one row per (Model, ROI, FC) combination. Model column
%   distinguishes M1_Primary (age+sex+FD[+eTIV]) from M2_Sensitivity
%   (M1 + EDUC + duration + LEDD + UPDRS-3).
rows = {};
for ri=1:nROI, for pi=1:nPred
    rows{end+1} = {'M1_Primary',    allNames{ri}, predNames{pi}, ...
        rMat(ri,pi),    pMat(ri,pi),    pFDR(ri,pi)}; %#ok
    rows{end+1} = {'M2_Sensitivity', allNames{ri}, predNames{pi}, ...
        rMat_M2(ri,pi), pMat_M2(ri,pi), pFDR_M2(ri,pi)}; %#ok
end, end
writetable(cell2table(vertcat(rows{:}), ...
    'VariableNames',{'Model','ROI','FC_Component','rho','p_unc','p_FDR'}), ...
    fullfile(outputDir,'Table_Structural_Correlates.csv'));
fprintf('  Saved: Table_Structural_Correlates.csv (M1 + M2, %d rows)\n', numel(rows));

%% =========================================================================
%  13. LONGITUDINAL HIPPOCAMPUS SLOPE x FC — hypothesis-driven follow-up
%      Motivation: Section 9 found RH SM-DA x hippocampal volume rho = +0.19
%      (within-component FDR-sig). This section asks whether FC predicts
%      prospective hippocampal atrophy rate.
%      Source: all_structural_measures_ses-24M.csv (24-month follow-up).
%      Slope = (V_24M - V_BL) / 2 yr (single follow-up visit, so multi-
%      timepoint LME reduces to the same 2-point difference).
%      Method: HC-normative eTIV correction at BOTH timepoints (Buckner
%      2004; Sanfilipo 2004): hc_beta_sub(1) for bilateral, hc_beta_lat(1:2)
%      for L/R. Slope = corrected_V24 - corrected_VBL implicitly removes
%      within-subject eTIV drift between scans (FreeSurfer eTIV is not
%      perfectly stable across rescans). RAW-volume slope variants
%      (vBL_*_raw, v24_*_raw) computed in parallel for sensitivity.
%      The cross-sectional analysis (sections 1 to 11) uses eTIV as a covariate; for the
%      slopes the HC-normative residual is used because eTIV drifts between scans.
%      Spearman partial r | age + sex + FD + BL_vol (M1); M2 adds education, disease
%      duration, LEDD and UPDRS-III.
%      The bilateral hippocampus is the primary test; left and right are examined
%      separately, each with BH FDR across components within side.
% =========================================================================
fprintf('\n[13] Longitudinal hippocampus slope x FC (hypothesis-driven)...\n');

%% 13.0 Baseline education (M2 covariate)
%  baselineEduc is defined by fc_cognitive_slope_mediation.m; it is rebuilt from Tlong here if
%  it is missing from the workspace.
if ~exist('baselineEduc','var')
    baselineEduc = nan(nFC,1);
    if exist('Tlong','var') && istable(Tlong)
        educCol = '';
        for cn = {'educ','EDUCYRS','education','education_years','EDUC'}
            if ismember(cn{1}, Tlong.Properties.VariableNames), educCol = cn{1}; break; end
        end
        if ~isempty(educCol)
            eV = double(Tlong.(educCol));
            if ismember('YEAR', Tlong.Properties.VariableNames)
                yV = double(Tlong.YEAR);
            elseif exist('timeVar','var') && ismember(timeVar, Tlong.Properties.VariableNames)
                yV = double(Tlong.(timeVar));
            else
                yV = zeros(height(Tlong),1);
            end
            pV = double(Tlong.PATNO);
            for si = 1:nFC
                pidx = pV == fcPATNO(si);
                if ~any(pidx), continue; end
                [~,ord] = sort(yV(pidx));
                valS = eV(pidx); valS = valS(ord);
                k = find(~isnan(valS),1);
                if ~isempty(k), baselineEduc(si) = valS(k); end
            end
            fprintf('  baselineEduc extracted inline from Tlong.%s (n=%d non-NaN)\n', ...
                educCol, sum(~isnan(baselineEduc)));
        else
            fprintf('  [warn] baselineEduc: no education column found in Tlong; all NaN\n');
        end
    else
        fprintf('  [warn] Tlong unavailable; baselineEduc will be all NaN\n');
    end
end

% HC normative eTIV slopes for the lateralised subcortical volumes (used for the hippocampus)
latSubPairs = { ...
    'L Hippocampus',  'fs_vol_Left_Hippocampus';  'R Hippocampus',  'fs_vol_Right_Hippocampus'; ...
    'L Amygdala',     'fs_vol_Left_Amygdala';     'R Amygdala',     'fs_vol_Right_Amygdala'; ...
    'L Caudate',      'fs_vol_Left_Caudate';       'R Caudate',      'fs_vol_Right_Caudate'; ...
    'L Putamen',      'fs_vol_Left_Putamen';       'R Putamen',      'fs_vol_Right_Putamen'; ...
    'L Accumbens',    'fs_vol_Left_Accumbens_area';'R Accumbens',    'fs_vol_Right_Accumbens_area'; ...
    'L Thalamus',     'fs_vol_Left_Thalamus';      'R Thalamus',     'fs_vol_Right_Thalamus'; ...
    'L Pallidum',     'fs_vol_Left_Pallidum';       'R Pallidum',     'fs_vol_Right_Pallidum'; ...
    'L Cereb Cortex', 'fs_vol_Left_Cerebellum_Cortex'; 'R Cereb Cortex', 'fs_vol_Right_Cerebellum_Cortex'; ...
    'Brain Stem',     'fs_vol_Brain_Stem'};
nLatSub = size(latSubPairs,1);
hc_beta_lat = nan(nLatSub,1);
for ri = 1:nLatSub
    asegCol = strrep(latSubPairs{ri,2}, 'fs_vol_', '');
    if ismember(asegCol, asegT.Properties.VariableNames)
        hcV = double(asegT.(asegCol)(isHC));
        ok = ~isnan(hcV) & ~isnan(hc_eTIV);
        if sum(ok)>15, b = polyfit(hc_eTIV(ok), hcV(ok), 1); hc_beta_lat(ri) = b(1);
        else, hc_beta_lat(ri) = 0; end
    else
        hc_beta_lat(ri) = 0;
    end
end

%% 13.1 Load 24M follow-up
structFile_24M = fullfile(getenv('PPMI_IMAGING_ROOT'), 'MRI_Data/derivatives/all_structural_measures_ses-24M.csv');
structT_24M = readtable(structFile_24M);
structPATNO_24M = double(structT_24M.PATNO);
structeTIV_24M = double(structT_24M.fs_EstimatedTotalIntraCranialVol);
fprintf('  24M follow-up: %d subjects\n', height(structT_24M));

% Nominal interval: ses-24M => 24 months (no per-scan dates available)
dtYears = 2;

%% 13.2 Match to FC sample
struct24Idx = nan(nFC,1);
for si = 1:nFC
    r = find(structPATNO_24M == fcPATNO(si), 1);
    if ~isempty(r), struct24Idx(si) = r; end
end
hasBoth = ~isnan(struct24Idx) & ~isnan(structIdx);
fprintf('  FC subjects with both BL and 24M: %d / %d\n', sum(hasBoth), nFC);

%% 13.3 Extract hippocampus volumes at BL and 24M
%   Volumes are HC-normative eTIV-corrected at both time points, so the slope removes the
%   within-subject eTIV drift between the baseline and 24-month scans (the cross-sectional
%   analysis in sections 1 to 11 uses eTIV as a covariate instead).
%   Raw-volume variants (vBL_*_raw, v24_*_raw) are computed alongside.
hippRow = find(strcmp(allNames,'Hippocampus'), 1);  % row 1 in bilateral ROI matrix
vBL_bi     = normROI(:, hippRow);    % PRIMARY HC-corrected bilateral
vBL_bi_raw = rawROI(:, hippRow);     % SENSITIVITY raw
v24_bi = nan(nFC,1);  v24_bi_raw = nan(nFC,1);
vBL_L  = nan(nFC,1);  v24_L = nan(nFC,1);
vBL_R  = nan(nFC,1);  v24_R = nan(nFC,1);
vBL_L_raw = nan(nFC,1); v24_L_raw = nan(nFC,1);
vBL_R_raw = nan(nFC,1); v24_R_raw = nan(nFC,1);
fc_eTIV_24M = nan(nFC,1);

for si = 1:nFC
    if ~isnan(struct24Idx(si))
        idx = struct24Idx(si);
        fc_eTIV_24M(si) = structeTIV_24M(idx);
        L24 = double(structT_24M.fs_vol_Left_Hippocampus(idx));
        R24 = double(structT_24M.fs_vol_Right_Hippocampus(idx));
        eTIV24 = structeTIV_24M(idx);
        % PRIMARY (HC-normative)
        v24_bi(si) = (L24 + R24)/2 - hc_beta_sub(1) * (eTIV24 - hc_mean_eTIV);
        v24_L(si)  = L24           - hc_beta_lat(1) * (eTIV24 - hc_mean_eTIV);
        v24_R(si)  = R24           - hc_beta_lat(2) * (eTIV24 - hc_mean_eTIV);
        % SENSITIVITY (raw)
        v24_bi_raw(si) = (L24 + R24)/2;
        v24_L_raw(si)  = L24;
        v24_R_raw(si)  = R24;
    end
    if ~isnan(structIdx(si))
        idx = structIdx(si);
        L_BL = double(structT.fs_vol_Left_Hippocampus(idx));
        R_BL = double(structT.fs_vol_Right_Hippocampus(idx));
        % PRIMARY (HC-normative)
        vBL_L(si) = L_BL - hc_beta_lat(1) * (fc_eTIV(si) - hc_mean_eTIV);
        vBL_R(si) = R_BL - hc_beta_lat(2) * (fc_eTIV(si) - hc_mean_eTIV);
        % SENSITIVITY (raw)
        vBL_L_raw(si) = L_BL;
        vBL_R_raw(si) = R_BL;
    end
end

slope_bi = (v24_bi - vBL_bi) / dtYears;
slope_L  = (v24_L  - vBL_L)  / dtYears;
slope_R  = (v24_R  - vBL_R)  / dtYears;

% Sensitivity-only: winsorised copies at 1/99 (primary slopes stay raw)
slope_bi_winsor = slope_bi;
slope_L_winsor  = slope_L;
slope_R_winsor  = slope_R;
for sv = {'slope_bi_winsor','slope_L_winsor','slope_R_winsor'}
    v = eval(sv{1}); okV = ~isnan(v);
    if sum(okV) < 10, continue; end
    p1 = prctile(v(okV), 1); p99 = prctile(v(okV), 99);
    v(okV & v < p1)  = p1;
    v(okV & v > p99) = p99;
    eval(sprintf('%s = v;', sv{1}));
end
fprintf('  Primary slopes: raw (NoWinsor). Winsor copies (1/99) saved as sensitivity.\n');

%% 13.4 Partial correlation: two models, each hipp x 4 FC
%       M1 (primary): age + sex + FD + BL_vol
%       M2 (sensitivity): M1 + education + duration + LEDD + UPDRS-III
%       Slopes from HC-normative-corrected volumes (within-subject eTIV drift removed), so
%       eTIV is not entered again as a covariate. BL_vol controls regression to the mean.
hippLabels = {'Bilateral','Left','Right'};
slopeSet   = {slope_bi, slope_L, slope_R};
blSet      = {vBL_bi,   vBL_L,   vBL_R};
nHipp = numel(hippLabels);

modelLabels = {'M1_Primary','M2_Sensitivity'};
modelDesc   = {'age+sex+FD+BLvol', ...
               'age+sex+FD+EDUC+dur+LEDD+UPDRS3+BLvol'};

rHipp        = nan(nHipp, nPred, 2);
pHipp        = nan(nHipp, nPred, 2);
pFDR_hipp    = nan(nHipp, nPred, 2);
nHippPerCell = nan(nHipp, nPred, 2);

for mi = 1:2
    for hi = 1:nHipp
        y  = slopeSet{hi};
        bl = blSet{hi};
        if mi == 1
            Zuse = [fcMean.age, sexFC, fcMean.mean_FD, bl];
        else
            Zuse = [fcMean.age, sexFC, fcMean.mean_FD, baselineEduc, ...
                    baselineDur, baselineLEDD, baselineUPDRS3, bl];
        end
        for pi = 1:nPred
            x = predData(:,pi);
            % partialcorr (rank everything incl. Z) — matches the paper's
            % primary Spearman-partial convention used throughout Figure 4.
            [rho, pval] = partialcorr(x, y, Zuse, 'Type','Spearman','Rows','complete');
            nOK = sum(~isnan(x) & ~isnan(y) & all(~isnan(Zuse),2));
            nHippPerCell(hi,pi,mi) = nOK;
            if nOK < 20, continue; end
            rHipp(hi,pi,mi) = rho;
            pHipp(hi,pi,mi) = pval;
        end
    end
    % Within-side BH FDR (4 FC per hipp) within each model
    for hi = 1:nHipp
        pv = pHipp(hi,:,mi); v = ~isnan(pv);
        [sP,sI] = sort(pv(v)); m = numel(sP);
        adj = nan(1,m); adj(m) = sP(m);
        for k = m-1:-1:1, adj(k) = min(adj(k+1), sP(k)*m/k); end
        adj = min(adj,1);
        tmp = nan(1,m); tmp(sI) = adj;
        out = nan(1,nPred); out(v) = tmp;
        pFDR_hipp(hi,:,mi) = out;
    end
end

%% 13.5 Report (both models)
for mi = 1:2
    fprintf('\n  === [%s: %s] ===\n', modelLabels{mi}, modelDesc{mi});
    for hi = 1:nHipp
        n_hi = max(nHippPerCell(hi,:,mi));
        fprintf('\n  %s Hippocampus  (n = %d):\n', hippLabels{hi}, n_hi);
        fprintf('  %-22s  %7s  %10s  %10s\n', 'FC component','rho','p_unc','p_FDR');
        fprintf('  %s\n', repmat('-',1,56));
        for pi = 1:nPred
            if isnan(rHipp(hi,pi,mi))
                fprintf('  %-22s  %7s  %10s  %10s\n', predNames{pi}, 'NaN', 'n<20', '-');
                continue
            end
            tag = '';
            if pFDR_hipp(hi,pi,mi) < 0.05,  tag = ' *FDR';
            elseif pHipp(hi,pi,mi) < 0.05,  tag = ' *unc'; end
            fprintf('  %-22s  %+.3f  %.4f  %.4f%s\n', predNames{pi}, ...
                rHipp(hi,pi,mi), pHipp(hi,pi,mi), pFDR_hipp(hi,pi,mi), tag);
        end
    end
end

%% 13.6 Save focused CSV (24 rows: 3 hipp x 4 FC x 2 models)
rows_sl = {};
for mi = 1:2, for hi = 1:nHipp, for pi = 1:nPred
    rows_sl{end+1} = {modelLabels{mi}, hippLabels{hi}, predNames{pi}, ...
        nHippPerCell(hi,pi,mi), rHipp(hi,pi,mi), pHipp(hi,pi,mi), pFDR_hipp(hi,pi,mi)}; %#ok
end, end, end
writetable(cell2table(vertcat(rows_sl{:}), ...
    'VariableNames',{'Model','Hippocampus','FC_Component','n','rho','p_unc','p_FDR_within_side'}), ...
    fullfile(outputDir,'Table_Hippocampus_Slope_x_FC.csv'));
fprintf('\n  Saved: Table_Hippocampus_Slope_x_FC.csv (both models, 24 rows)\n');

% Per-subject longitudinal data for fc_cognitive_slope_mediation.m panel (g)
TlongHipp = table(fcPATNO, slope_bi, slope_L, slope_R, vBL_bi, vBL_L, vBL_R, fc_eTIV, ...
    'VariableNames', {'PATNO','slope_bilateral','slope_left','slope_right', ...
                      'vBL_bilateral','vBL_left','vBL_right','eTIV_BL'});
writetable(TlongHipp, fullfile(outputDir,'Table_Hippocampus_Longitudinal_Data.csv'));
fprintf('  Saved: Table_Hippocampus_Longitudinal_Data.csv (per-subject for panel g, includes BL eTIV)\n');

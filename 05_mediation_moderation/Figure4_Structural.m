%% =========================================================================
%  STRUCTURAL CORRELATES — FreeSurfer volumes x FC components
%  Source: all_structural_measures.csv (bilateral-averaged subcort vol +
%          cortical gray matter vol).
%  Method: UNIFORM eTIV-as-covariate — Spearman partial
%          r | age + sex + FD + eTIV  (subcort and cortical treated identically;
%          slope and lateralised analyses use the same covariate set).
%          HC-normative residual correction is computed in parallel and used
%          for the §7b sensitivity comparison only.
%  Spearman partial via partialcorr; within-component BH FDR.
%  Requires workspace from Figure4.m.
% =========================================================================
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

% Apply HC normative correction to subcortical volumes (SENSITIVITY only,
% post 2026-05-06; primary now uses rawROI + eTIV-in-Z uniformly).
normROI = rawROI;
for ri=1:nSub
    for si=1:nFC
        if ~isnan(rawROI(si,ri)) && ~isnan(fc_eTIV(si)) && ~isnan(hc_beta_sub(ri))
            normROI(si,ri) = rawROI(si,ri) - hc_beta_sub(ri) * (fc_eTIV(si) - hc_mean_eTIV);
        end
    end
end
fprintf('  Primary volumes: rawROI (uncorrected); eTIV partialled in covariate set.\n');
fprintf('  Sensitivity:     HC-normative-corrected subcortical (normROI) — see §7b.\n');

%% 6. Predictors + covariates
predData  = nan(nFC, nAllFC);
predNames = cell(1, nAllFC);
for fi=1:nAllFC
    predData(:,fi) = fcData(:,fi);   % suprathreshold (from Figure4.m)
    predNames{fi}  = fcShort{fi};
end
nPred = nAllFC;

% Covariates — two-model sensitivity design (uniform method post 2026-05-06):
%   M1 PRIMARY    : age + sex + FD + eTIV       (uniform across subcort + cortical)
%   M2 SENSITIVITY: M1 + EDUC + duration + LEDD + UPDRS-3
%
% Rationale: switched to eTIV-as-covariate as the primary method across the
% board (was HC-normative residual for subcortical only — see §7b for the
% empirical comparison). This eliminates the cortical/subcortical asymmetry
% and gives a single uniform method for all 43 ROIs. HC-normative correction
% is retained as a parallel sensitivity (rMat_subHC in §7b).
% age^2 was previously computed but unused (kept as ageSq for back-compat).
ageSq  = fcMean.age .^ 2;
Z_sub  = [fcMean.age, sexFC, fcMean.mean_FD, fc_eTIV];
Z_cort = [fcMean.age, sexFC, fcMean.mean_FD, fc_eTIV];

% Defensive existence checks for M2 sensitivity covariates (Figure4.m
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

%% 7b. SENSITIVITY: subcortical with HC-normative residual
%   Recomputes subcortical correlations using normROI (HC-normative correction)
%   while §7 primary uses rawROI + eTIV-in-Z. This is the methodological
%   robustness check for the May 2026 switch to uniform covariate method.
fprintf('\n  [Sensitivity §7b] Subcortical with HC-normative residual (vs primary covariate method)\n');
Z_sub_HC = [fcMean.age, sexFC, fcMean.mean_FD];   % NO eTIV — already in normROI

rMat_subHC = nan(nSub, nPred);  pMat_subHC = nan(nSub, nPred);
for ri = 1:nSub
    y = normROI(:,ri);   % HC-NORMATIVE-corrected subcortical volume
    for pi = 1:nPred
        x = predData(:,pi);
        ok = ~isnan(x) & ~isnan(y) & all(~isnan(Z_sub_HC),2);
        if sum(ok) < 20, continue; end
        [rho, pval] = partialcorr(x(ok), y(ok), Z_sub_HC(ok,:), 'Type','Spearman');
        rMat_subHC(ri,pi) = rho;  pMat_subHC(ri,pi) = pval;
    end
end

% Side-by-side comparison for each FC component, focused on subcortical only
fprintf('\n  Subcortical: PRIMARY (eTIV-as-covariate)  vs  SENSITIVITY (HC-normative)\n');
maxAbsDelta_all = 0;  maxCell = '';
for pi = 1:nPred
    fprintf('\n  --- %s ---\n', predNames{pi});
    fprintf('  %-15s  %8s  %8s   %8s  %8s   %8s\n', ...
        'ROI', 'rho_cov', 'p_cov', 'rho_HC', 'p_HC', 'Δrho');
    fprintf('  %s\n', repmat('-',1,68));
    maxAbsDelta_pi = 0;
    for ri = 1:nSub
        rho_cov = rMat(ri, pi);          % PRIMARY (covariate method, from §7)
        p_cov   = pMat(ri, pi);
        rho_hc  = rMat_subHC(ri, pi);
        p_hc    = pMat_subHC(ri, pi);
        if isnan(rho_cov) && isnan(rho_hc), continue; end
        dRho = rho_hc - rho_cov;
        if abs(dRho) > maxAbsDelta_pi, maxAbsDelta_pi = abs(dRho); end
        if abs(dRho) > maxAbsDelta_all
            maxAbsDelta_all = abs(dRho);
            maxCell = sprintf('%s × %s', subNames{ri}, predNames{pi});
        end
        fprintf('  %-15s  %+8.3f  %8.4f   %+8.3f  %8.4f   %+8.3f\n', ...
            subNames{ri}, rho_cov, p_cov, rho_hc, p_hc, dRho);
    end
    fprintf('  Max |Δρ| (%s): %.3f\n', predNames{pi}, maxAbsDelta_pi);
end
fprintf('\n  Across all subcortical × FC cells: max |Δρ| = %.3f  (%s)\n', ...
    maxAbsDelta_all, maxCell);
fprintf('  Interpretation: methods are equivalent if |Δρ| < ~0.02 across all cells.\n');
fprintf('  Larger divergence in any cell flags a per-region nonlinearity in the eTIV-volume\n');
fprintf('  relationship that the covariate method may under-correct for.\n');

%% 8. Within-component FDR — RH SM-DA as primary (fi=2)
%    FDR across 43 ROIs for each FC component separately. Done within
%    each model (M1 / M2) to keep families internally consistent. Main
%    analysis and forest plot focus on M1 RH SM-DA (the cognitive mediator).
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

%% 10. Forest plot — RH SM-DA component, compact single panel
cSig  = [0.20 0.40 0.65];   % blue for FDR-significant
cUnc  = [0.45 0.60 0.80];   % lighter blue for uncorrected
cNS   = [0.78 0.78 0.78];   % gray for n.s.
cBG   = [0.97 0.97 0.97];   % alternating row shade

rhoVec = rMat(:, piSMDA);
pVec   = pMat(:, piSMDA);
fdrVec = pFDR(:, piSMDA);

fig_f = figure('Color','w','Position',[100 80 420 max(400, 11*nROI)],'Visible','off');
ax = axes('Parent',fig_f,'Units','normalized','Position',[0.38 0.08 0.56 0.86]);
hold(ax,'on');

for ri = 1:nROI
    y = nROI - ri + 1;
    rho = rhoVec(ri);
    if isnan(rho), continue; end

    isFDR = ~isnan(fdrVec(ri)) && fdrVec(ri) < 0.05;
    isUnc = ~isnan(pVec(ri))   && pVec(ri) < 0.01;

    % Alternating row shading
    if mod(ri,2)==0
        fill(ax, [-0.35 0.35 0.35 -0.35], [y-0.45 y-0.45 y+0.45 y+0.45], ...
            cBG, 'EdgeColor','none','FaceAlpha',0.5);
    end

    % Horizontal line from 0 to rho (lollipop stem)
    if isFDR
        plot(ax, [0 rho], [y y], '-', 'Color', cSig, 'LineWidth', 1.8);
        plot(ax, rho, y, 'o', 'MarkerSize', 8, 'MarkerFaceColor', cSig, ...
            'MarkerEdgeColor', 'w', 'LineWidth', 0.8);
        if rho >= 0, ha = 'left'; else, ha = 'right'; end
        text(ax, rho + sign(rho)*0.015, y, sprintf('%.2f', rho), ...
            'FontSize', 6, 'FontWeight', 'bold', 'Color', cSig, ...
            'VerticalAlignment', 'middle', 'HorizontalAlignment', ha);
    elseif isUnc
        plot(ax, [0 rho], [y y], '-', 'Color', cUnc, 'LineWidth', 1.2);
        plot(ax, rho, y, 'o', 'MarkerSize', 6.5, 'MarkerFaceColor', 'w', ...
            'MarkerEdgeColor', cUnc, 'LineWidth', 1.4);
        text(ax, rho + sign(rho)*0.015, y, sprintf('%.2f', rho), ...
            'FontSize', 5.5, 'Color', [0.35 0.35 0.35], ...
            'VerticalAlignment', 'middle');
    else
        plot(ax, [0 rho], [y y], '-', 'Color', [0.85 0.85 0.85], 'LineWidth', 0.6);
        plot(ax, rho, y, '.', 'MarkerSize', 6, 'Color', cNS);
    end
end

xline(ax, 0, '-', 'Color', [0.4 0.4 0.4], 'LineWidth', 0.8);
% Subcortical / cortical separator
yline(ax, nROI - nSub + 0.5, '-', 'Color', [0.3 0.3 0.3], 'LineWidth', 0.8);
text(ax, -0.33, nROI - nSub/2 + 0.5, 'Subcortical', 'FontSize', 6, ...
    'FontWeight','bold', 'Color',[0.4 0.4 0.4], 'HorizontalAlignment','center', 'Rotation', 90);
text(ax, -0.33, (nROI - nSub)/2 + 0.5, 'Cortical', 'FontSize', 6, ...
    'FontWeight','bold', 'Color',[0.4 0.4 0.4], 'HorizontalAlignment','center', 'Rotation', 90);

% Y-axis labels with formatting
yLabels = cell(nROI,1);
for ri = 1:nROI
    ori = nROI - ri + 1;
    lbl = allNames{ori};
    isFDR = ~isnan(fdrVec(ori)) && fdrVec(ori) < 0.05;
    isUnc = ~isnan(pVec(ori))   && pVec(ori) < 0.01;
    if isFDR,     yLabels{ri} = ['\bf\color[rgb]{0.20,0.40,0.65}' lbl];
    elseif isUnc, yLabels{ri} = ['\it\color[rgb]{0.35,0.35,0.35}' lbl];
    else,         yLabels{ri} = ['\color[rgb]{0.55,0.55,0.55}' lbl];
    end
end

set(ax, 'YTick', 1:nROI, 'YTickLabel', yLabels, 'FontSize', 6.5, ...
    'TickDir', 'out', 'Box', 'off', 'YLim', [0.3 nROI+0.7], ...
    'XLim', [-0.35 0.35], 'FontName', 'Helvetica');
xlabel(ax, 'Spearman partial \rho', 'FontSize', 9);
title(ax, sprintf('RH Somatomotor-DorsAttn FC  \\times  Regional Volume'), ...
    'FontSize', 10, 'FontWeight', 'bold', 'Color', cSig);

% Legend
annotation(fig_f, 'textbox', [0.40 0.01 0.55 0.05], 'String', ...
    sprintf('\\bf\\color[rgb]{0.20,0.40,0.65}\\bullet FDR<.05 (%d)   \\rm\\color[rgb]{0.45,0.60,0.80}\\circ p<.01 (%d)   \\color[rgb]{0.6,0.6,0.6}\\cdot n.s.   |   %d ROIs tested', ...
    nSigFDR, nSigUnc, nTests_smda), ...
    'EdgeColor', 'none', 'FontSize', 6.5, 'HorizontalAlignment', 'center', ...
    'VerticalAlignment', 'bottom', 'FitBoxToText', 'on');

exportgraphics(fig_f, fullfile(outputDir,'Fig4_Structural_Forest.pdf'), 'ContentType','vector');
exportgraphics(fig_f, fullfile(outputDir,'Fig4_Structural_Forest.png'), 'Resolution',300);
fprintf('  Saved: Fig4_Structural_Forest.pdf / .png\n');
close(fig_f);

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

%% 12. LATERALISED ANALYSIS — all L/R volumes x all FC components
%    No bilateral averaging. Includes subcortical (fs_vol_Left/Right_*)
%    and cortical GM volumes (fs_lh/rh_*_volume). Method (post 2026-05-06):
%    UNIFORM eTIV-as-covariate (Z_sub = Z_cort = [age, sex, FD, eTIV]) on
%    raw lateralised volumes (latROI). HC-normative-corrected lateralised
%    subcortical (latROI_HC) is preserved for parallel sensitivity.
%    FDR reported both within-component and across all components.
fprintf('\n[12] Lateralised volumes x FC (no bilateral averaging)...\n');

% --- Build lateralised subcortical ROI list ---
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

% HC normative beta per lateralised subcortical ROI
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

% --- Build lateralised cortical ROI list (GM volume) ---
latCortNames = {};
latCortCols  = {};
for ri = 1:nCort
    for hemi = {'lh','rh'}
        h = hemi{1};
        if strcmp(h,'lh'), pre = 'L '; else, pre = 'R '; end
        latCortNames{end+1} = [pre aparcReg{ri}]; %#ok
        latCortCols{end+1}  = ['fs_' h '_' aparcReg{ri} '_volume']; %#ok
    end
end
nLatCort = numel(latCortNames);
nLatROI  = nLatSub + nLatCort;
latNames = [latSubPairs(:,1); latCortNames(:)];
fprintf('  Lateralised ROIs: %d subcort + %d cortical = %d\n', nLatSub, nLatCort, nLatROI);

% --- Extract lateralised values for FC sample ---
%   Primary uses raw volumes; eTIV is partialled via Z_sub / Z_cort which
%   both include fc_eTIV (uniform covariate method). HC-normative variants
%   are kept inline as latROI_HC for parallel sensitivity comparison.
latROI    = nan(nFC, nLatROI);
latROI_HC = nan(nFC, nLatROI);   % subcortical only; cortical = same as latROI
for si = 1:nFC
    if isnan(structIdx(si)), continue; end
    idx = structIdx(si);
    for ri = 1:nLatSub
        col = latSubPairs{ri,2};
        if ismember(col, structT.Properties.VariableNames)
            raw = double(structT.(col)(idx));
            latROI(si,ri)    = raw;                                                       % PRIMARY
            latROI_HC(si,ri) = raw - hc_beta_lat(ri) * (fc_eTIV(si) - hc_mean_eTIV);      % SENSITIVITY
        end
    end
    for ri = 1:nLatCort
        col = latCortCols{ri};
        if ismember(col, structT.Properties.VariableNames)
            latROI(si, nLatSub+ri)    = double(structT.(col)(idx));
            latROI_HC(si, nLatSub+ri) = double(structT.(col)(idx));   % cortical: same in both
        end
    end
end

% --- Spearman partial correlation ---
rLat = nan(nLatROI, nPred); pLat = nan(nLatROI, nPred);
for ri = 1:nLatROI
    y = latROI(:,ri);
    if ri <= nLatSub, Zuse = Z_sub; else, Zuse = Z_cort; end
    for pi = 1:nPred
        x = predData(:,pi);
        ok = ~isnan(x) & ~isnan(y) & all(~isnan(Zuse),2);
        if sum(ok)<20, continue; end
        [rho,pval] = partialcorr(x(ok), y(ok), Zuse(ok,:), 'Type','Spearman');
        rLat(ri,pi) = rho; pLat(ri,pi) = pval;
    end
end

% --- FDR: within-component ---
pFDR_lat_wc = nan(nLatROI, nPred);
for pi = 1:nPred
    pv = pLat(:,pi); v = ~isnan(pv);
    [sP,sI] = sort(pv(v)); m = numel(sP);
    adj = nan(m,1); adj(m) = sP(m);
    for k = m-1:-1:1, adj(k) = min(adj(k+1), sP(k)*m/k); end
    adj = min(adj,1);
    tmp = nan(m,1); tmp(sI) = adj;
    out = nan(nLatROI,1); out(v) = tmp;
    pFDR_lat_wc(:,pi) = out;
end

% --- FDR: across ALL components ---
pFDR_lat_all = nan(nLatROI, nPred);
pFlat = pLat(:); v = ~isnan(pFlat);
[sP,sI] = sort(pFlat(v)); m = numel(sP);
adj = nan(m,1); adj(m) = sP(m);
for k = m-1:-1:1, adj(k) = min(adj(k+1), sP(k)*m/k); end
adj = min(adj,1);
tmp = nan(m,1); tmp(sI) = adj;
pFDR_flat = nan(size(pFlat)); pFDR_flat(v) = tmp;
pFDR_lat_all = reshape(pFDR_flat, nLatROI, nPred);

nTotal = sum(~isnan(pLat(:)));
fprintf('  Total tests: %d\n', nTotal);
for pi = 1:nPred
    nWC = sum(pFDR_lat_wc(:,pi)<0.05);
    nAll = sum(pFDR_lat_all(:,pi)<0.05);
    fprintf('    %-25s  FDR-within: %d   FDR-across: %d   p<.01: %d\n', ...
        predNames{pi}, nWC, nAll, sum(pLat(:,pi)<0.01));
end

% --- Print significant results ---
fprintf('\n  Lateralised significant (FDR-within < .05 OR p_unc < .005):\n');
fprintf('  %-25s  %-22s  %7s  %8s  %8s  %8s\n', 'ROI','FC','rho','p_unc','FDR_wc','FDR_all');
fprintf('  %s\n', repmat('-',1,98));
for ri = 1:nLatROI, for pi = 1:nPred
    isFDR_wc  = ~isnan(pFDR_lat_wc(ri,pi))  && pFDR_lat_wc(ri,pi) < 0.05;
    isFDR_all = ~isnan(pFDR_lat_all(ri,pi)) && pFDR_lat_all(ri,pi) < 0.05;
    isUnc     = ~isnan(pLat(ri,pi)) && pLat(ri,pi) < 0.005;
    if isFDR_wc || isFDR_all || isUnc
        tag = '';
        if isFDR_all, tag = ' +++ALL';
        elseif isFDR_wc, tag = ' ++WC';
        elseif isUnc, tag = ' *unc'; end
        fprintf('  %-25s  %-22s  %+.3f  %.4f  %.4f  %.4f%s\n', ...
            latNames{ri}, predNames{pi}, rLat(ri,pi), pLat(ri,pi), ...
            pFDR_lat_wc(ri,pi), pFDR_lat_all(ri,pi), tag);
    end
end, end

% --- Save lateralised table ---
rows_lat = {};
for ri = 1:nLatROI, for pi = 1:nPred
    rows_lat{end+1} = {latNames{ri}, predNames{pi}, rLat(ri,pi), pLat(ri,pi), ...
        pFDR_lat_wc(ri,pi), pFDR_lat_all(ri,pi)}; %#ok
end, end
writetable(cell2table(vertcat(rows_lat{:}), ...
    'VariableNames',{'ROI','FC_Component','rho','p_unc','pFDR_within','pFDR_across'}), ...
    fullfile(outputDir,'Table_Structural_Lateralised.csv'));
fprintf('  Saved: Table_Structural_Lateralised.csv\n');

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
%      Note the asymmetry vs §1-§12 cross-sectional, which use the uniform
%      eTIV-as-covariate method — for single-timepoint volumes the two
%      methods are statistically equivalent (§7b confirms |Δρ| ≲ 0.02);
%      for slope outcomes they diverge because the covariate method leaves
%      within-subject eTIV-drift noise in the slope.
%      Spearman partial r | age + sex + FD + EDUC + BL_vol  (M1).
%      Primary: bilateral hippocampus x 4 FC (4-test within-family BH FDR).
%      Sensitivity: L and R separately (each own 4-test within-side FDR).
%      Laterality motivated by preliminary results (L >> R).
% =========================================================================
fprintf('\n[13] Longitudinal hippocampus slope x FC (hypothesis-driven)...\n');

%% 13.0 Defensive extraction: baselineEduc (added 2026-04-20)
%  If Figure4.m was not re-run after the covCandsMap update, baselineEduc
%  will not exist in the workspace. Pull it from Tlong inline here.
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
%   PRIMARY: HC-normative-corrected at both timepoints. The slope
%   constructed from HC-corrected volumes implicitly removes within-subject
%   eTIV drift between BL and 24M scans (FreeSurfer eTIV is not perfectly
%   stable across rescans). The eTIV-as-covariate approach used uniformly
%   in §1-§12 is statistically equivalent for SINGLE-timepoint volumes
%   (§7b: |Δρ| ≲ 0.02) but NOT for slope outcomes — leaving within-subject
%   eTIV drift as residual noise in the slope. So §13/§14/§15 keep the
%   HC-normative residual as primary while §1-§12 use the uniform covariate
%   method (cross-sectional / longitudinal asymmetry is principled, not
%   inconsistent).
%   SENSITIVITY: raw volumes (vBL_*_raw, v24_*_raw) retained for parallel.
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
%       M1 (PRIMARY):     age + sex + FD + EDUC + BL_vol
%       M2 (SENSITIVITY): M1 + duration + LEDD + UPDRS-3
%       Slopes from HC-normative-corrected volumes (within-subject eTIV
%       drift implicitly removed). NO explicit eTIV-as-covariate — that
%       would over-control and add noise (different from §1-§12 cross-
%       sectional, which DO use eTIV-as-covariate; see §13.3 docstring).
%       BL_vol kept as regression-to-mean control.
hippLabels = {'Bilateral','Left','Right'};
slopeSet   = {slope_bi, slope_L, slope_R};
blSet      = {vBL_bi,   vBL_L,   vBL_R};
nHipp = numel(hippLabels);

modelLabels = {'M1_Primary','M2_Sensitivity'};
modelDesc   = {'age+sex+FD+EDUC+BLvol', ...
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
            Zuse = [fcMean.age, sexFC, fcMean.mean_FD, baselineEduc, bl];
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

% Per-subject longitudinal data for Figure4.m panel (g)
TlongHipp = table(fcPATNO, slope_bi, slope_L, slope_R, vBL_bi, vBL_L, vBL_R, fc_eTIV, ...
    'VariableNames', {'PATNO','slope_bilateral','slope_left','slope_right', ...
                      'vBL_bilateral','vBL_left','vBL_right','eTIV_BL'});
writetable(TlongHipp, fullfile(outputDir,'Table_Hippocampus_Longitudinal_Data.csv'));
fprintf('  Saved: Table_Hippocampus_Longitudinal_Data.csv (per-subject for panel g, includes BL eTIV)\n');

%% =========================================================================
%  14. TARGETED LINEAR MODEL — Δ_region ~ BL_FC + (DAT-parity covariates)
%      Hypothesis-driven parametric test on 2 FC components (preliminary
%      signal): Ap(-)1 (RH SM-DA) and Diff(-)1 (Ap-dom SM-DA), each x 3
%      hippocampus regions (bilateral, L, R) = 6 models. No FDR.
%      Covariates (match DAT slope methodology; DAT-specific dynCovs
%      omitted: apathy_bl, QUIP_bl, DAT_bl, dep_bl):
%        age + sex + FD + EDUCYRS + duration + LEDD + UPDRS-3 + BL_vol.
%      Primary uses raw (NoWinsor) slopes from §13. Winsorised copies
%      (slope_*_winsor) are available for sensitivity analyses only.
%      Note: Δ is scalar per subject, so (1|PATNO) would be singleton and
%      unidentified — fitlm (OLS) is the correct collapsed form.
% =========================================================================
fprintf('\n[14] Targeted linear models: Delta-hippocampus ~ BL_FC (6 fits x 2 models)...\n');
fprintf('     M1 Primary: age+sex+FD+EDUC+BLvol   |   M2 Sensitivity: +dur+LEDD+UPDRS3\n');

fcTestIdx = [2, 3];            % 2 = Ap(-)1 RH SM-DA; 3 = Diff(-)1 Ap-dom SM-DA
fcTestLbl = predNames(fcTestIdx);

%% 14.1 Main 6 fits x 2 covariate models with diagnostics
lmRows = {};
for mi = 1:2
    fprintf('\n  ============================================================\n');
    fprintf('  [%s: %s]\n', modelLabels{mi}, modelDesc{mi});
    fprintf('  ============================================================\n');
    for ti = 1:numel(fcTestIdx)
        pi = fcTestIdx(ti);
        for hi = 1:nHipp
            if mi == 1
                Tlm = table(slopeSet{hi}, predData(:,pi), blSet{hi}, ...
                            fcMean.age, sexFC, fcMean.mean_FD, baselineEduc, ...
                    'VariableNames', {'Delta','BL_FC','BL_vol','age','sex','FD','EDUC'});
                mdlFmla = 'Delta ~ BL_FC + BL_vol + age + sex + FD + EDUC';
            else
                Tlm = table(slopeSet{hi}, predData(:,pi), blSet{hi}, ...
                            fcMean.age, sexFC, fcMean.mean_FD, baselineEduc, ...
                            baselineDur, baselineLEDD, baselineUPDRS3, ...
                    'VariableNames', {'Delta','BL_FC','BL_vol','age','sex','FD', ...
                                      'EDUC','dur','LEDD','UPDRS'});
                mdlFmla = ['Delta ~ BL_FC + BL_vol + age + sex + FD + EDUC ' ...
                           '+ dur + LEDD + UPDRS'];
            end
            ok = ~any(isnan(Tlm{:,:}), 2);
            mdl = fitlm(Tlm(ok,:), mdlFmla);

            coef  = mdl.Coefficients;
            b1    = coef{'BL_FC','Estimate'};
            se    = coef{'BL_FC','SE'};
            tval  = coef{'BL_FC','tStat'};
            pval  = coef{'BL_FC','pValue'};
            df    = mdl.DFE;
            n_fit = mdl.NumObservations;
            ci_lo = b1 - tinv(0.975, df) * se;
            ci_hi = b1 + tinv(0.975, df) * se;

            if mi == 1
                Z_part = [blSet{hi}(ok), fcMean.age(ok), sexFC(ok), ...
                          fcMean.mean_FD(ok), baselineEduc(ok)];
                Z_part_full = [blSet{hi}, fcMean.age, sexFC, ...
                               fcMean.mean_FD, baselineEduc];
            else
                Z_part = [blSet{hi}(ok), fcMean.age(ok), sexFC(ok), ...
                          fcMean.mean_FD(ok), baselineEduc(ok), ...
                          baselineDur(ok), baselineLEDD(ok), baselineUPDRS3(ok)];
                Z_part_full = [blSet{hi}, fcMean.age, sexFC, ...
                               fcMean.mean_FD, baselineEduc, ...
                               baselineDur, baselineLEDD, baselineUPDRS3];
            end
            % Pearson partial unchanged (raw-space); Spearman uses partialcorr.
            [rP, pP] = partialcorr(predData(ok,pi), slopeSet{hi}(ok), Z_part, 'Type','Pearson');
            [rS, pS] = partialcorr(predData(:,pi), slopeSet{hi}, Z_part_full, 'Type','Spearman','Rows','complete');

            cooks    = mdl.Diagnostics.CooksDistance;
            cookThr  = 4 / n_fit;
            influIdx = find(cooks > cookThr);
            okSubjIdx = find(ok);
            influSubj = okSubjIdx(influIdx);
            if ~isempty(influSubj)
                influPATNO = fcPATNO(influSubj);
                [~, sortOrd] = sort(cooks(influIdx), 'descend');
                influPATNO = influPATNO(sortOrd);
                cooksSorted = cooks(influIdx(sortOrd));
            else
                influPATNO = []; cooksSorted = [];
            end

            fprintf('\n  [%s Hippocampus] ~ [%s]   (n = %d, df = %d):\n', ...
                hippLabels{hi}, fcTestLbl{ti}, n_fit, df);
            fprintf('    OLS       beta = %+9.3f   SE = %7.3f   95%% CI [%+.3f, %+.3f]\n', ...
                b1, se, ci_lo, ci_hi);
            fprintf('              t(%d) = %+.3f   p = %.4f   R^2_adj = %.3f\n', ...
                df, tval, pval, mdl.Rsquared.Adjusted);
            fprintf('    Pearson   r_partial   = %+.3f   p = %.4f\n', rP, pP);
            fprintf('    Spearman  rho_partial = %+.3f   p = %.4f\n', rS, pS);
            fprintf('    Cook''s D  max = %.3f   thr = 4/n = %.3f   flagged = %d', ...
                max(cooks), cookThr, numel(influIdx));
            if ~isempty(influPATNO)
                topK = min(5, numel(influPATNO));
                pairs = arrayfun(@(k) sprintf('%d (D=%.2f)', influPATNO(k), cooksSorted(k)), ...
                    1:topK, 'UniformOutput', false);
                fprintf('   top: %s\n', strjoin(pairs, ', '));
            else
                fprintf('\n');
            end

            lmRows{end+1} = {modelLabels{mi}, hippLabels{hi}, fcTestLbl{ti}, n_fit, df, ...
                b1, se, ci_lo, ci_hi, tval, pval, mdl.Rsquared.Adjusted, ...
                rP, pP, rS, pS, max(cooks), numel(influIdx)}; %#ok
        end
    end
end

%% 14.2 Save extended CSV (12 fits: 6 per model x 2 models)
writetable(cell2table(vertcat(lmRows{:}), ...
    'VariableNames',{'Model','Hippocampus','FC','n','df','beta_BL_FC','SE', ...
    'CI_lo','CI_hi','t','p','R2_adj', ...
    'Pearson_r_partial','Pearson_p','Spearman_rho_partial','Spearman_p', ...
    'Cook_max','n_influential_4_over_n'}), ...
    fullfile(outputDir,'Table_Hippocampus_LM_Targeted.csv'));
fprintf('\n  Saved: Table_Hippocampus_LM_Targeted.csv (both models, 12 fits)\n');

%% =========================================================================
%  15. UNIFIED PREDICTION — baseline RH SM-DA FC → three outcomes
%      Tests whether baseline FC (Ap(-)1, idxMediation) independently
%      predicts three prospective outcomes in a single common framework:
%        (i)   Global cognitive slope     (slopeAll(:,1))
%        (ii)  DAT slope                  (dynSlopes(:,3))
%        (iii) Left hippocampal slope     (slope_L)
%      Same covariate set across all three (M1: age+sex+FD+EDUC;
%      M2: + duration+LEDD+UPDRS-3). Spearman partial + OLS.
%      Purpose: demonstrate FC as common upstream predictor across
%      dopaminergic, structural, and cognitive trajectories.
% =========================================================================
fprintf('\n[15] Unified 3-outcome prediction: RH SM-DA FC x {cog, DAT, hipp} slopes...\n');

fcV_unif = fcData(:, idxMediation);

outcomeLbl = {'Cog slope (global)','DAT slope','L Hipp slope'};
outcomeY   = {slopeAll(:,1), dynSlopes(:,3), slope_L};

unifRows = {};
for oi = 1:numel(outcomeLbl)
    yU = outcomeY{oi};
    for mi = 1:2
        if mi == 1
            Zuse = [fcMean.age, sexFC, fcMean.mean_FD, baselineEduc];
            covLabels = {'age','sex','FD','EDUC'};
            mdlFmla   = 'Outcome ~ BL_FC + age + sex + FD + EDUC';
        else
            Zuse = [fcMean.age, sexFC, fcMean.mean_FD, baselineEduc, ...
                    baselineDur, baselineLEDD, baselineUPDRS3];
            covLabels = {'age','sex','FD','EDUC','dur','LEDD','UPDRS'};
            mdlFmla   = ['Outcome ~ BL_FC + age + sex + FD + EDUC ' ...
                         '+ dur + LEDD + UPDRS'];
        end
        ok = ~isnan(fcV_unif) & ~isnan(yU) & all(~isnan(Zuse),2);
        n_unif = sum(ok);
        if n_unif < 20
            fprintf('  [%s | %s] n = %d — too small, skipping\n', outcomeLbl{oi}, modelLabels{mi}, n_unif);
            continue
        end

        % Spearman partial
        [rhoU, pU_sp] = partialcorr(fcV_unif(ok), yU(ok), Zuse(ok,:), 'Type','Spearman');

        % OLS with named variables
        if mi == 1
            Tunif = table(yU, fcV_unif, fcMean.age, sexFC, fcMean.mean_FD, baselineEduc, ...
                'VariableNames', {'Outcome','BL_FC','age','sex','FD','EDUC'});
        else
            Tunif = table(yU, fcV_unif, fcMean.age, sexFC, fcMean.mean_FD, baselineEduc, ...
                          baselineDur, baselineLEDD, baselineUPDRS3, ...
                'VariableNames', {'Outcome','BL_FC','age','sex','FD','EDUC', ...
                                  'dur','LEDD','UPDRS'});
        end
        ok_tbl = ~any(isnan(Tunif{:,:}), 2);
        mdlU = fitlm(Tunif(ok_tbl,:), mdlFmla);
        coefU = mdlU.Coefficients;
        b1U   = coefU{'BL_FC','Estimate'};
        seU   = coefU{'BL_FC','SE'};
        pU_ols= coefU{'BL_FC','pValue'};
        dfU   = mdlU.DFE;
        ciLo  = b1U - tinv(0.975, dfU)*seU;
        ciHi  = b1U + tinv(0.975, dfU)*seU;

        fprintf('\n  [%s | %s]  n = %d  df = %d\n', outcomeLbl{oi}, modelLabels{mi}, n_unif, dfU);
        fprintf('    Spearman  rho = %+.3f   p = %.4f\n', rhoU, pU_sp);
        fprintf('    OLS       beta = %+9.4f   SE = %7.4f   95%% CI [%+.4f, %+.4f]   p = %.4f\n', ...
            b1U, seU, ciLo, ciHi, pU_ols);

        unifRows{end+1} = {modelLabels{mi}, outcomeLbl{oi}, n_unif, dfU, ...
            rhoU, pU_sp, b1U, seU, ciLo, ciHi, pU_ols}; %#ok
    end
end

writetable(cell2table(vertcat(unifRows{:}), ...
    'VariableNames',{'Model','Outcome','n','df','rho_Spearman','p_Spearman', ...
    'beta_OLS','SE','CI_lo','CI_hi','p_OLS'}), ...
    fullfile(outputDir,'Table_Unified_FC_Predictions.csv'));
fprintf('\n  Saved: Table_Unified_FC_Predictions.csv\n');

%% =========================================================================
%  16. PARALLEL MEDIATION — FC → {DAT slope, L hipp slope} → Cog slope
%      Two mediators tested simultaneously. Preacher & Hayes percentile
%      bootstrap, 5,000 resamples, 95% CI on indirect effects.
%      Covariates (M1 primary): age + sex + FD + EDUC.
%      Paths:
%        a1: FC → DAT slope            (OLS, partial)
%        a2: FC → L hipp slope         (OLS, partial)
%        b1: DAT slope → cog slope     (partial on FC, hipp, covariates)
%        b2: hipp slope → cog slope    (partial on FC, DAT, covariates)
%        c:  total FC → cog slope
%        c': direct FC → cog slope after both mediators
%      Specific indirects: ab1 = a1*b1, ab2 = a2*b2
%      Total indirect = ab1 + ab2.  Significance via 95% bootstrap CI.
% =========================================================================
fprintf('\n[16] Parallel mediation: FC -> {DAT slope, L hipp slope} -> Cog slope...\n');

nBoot = 5000;
rng(42, 'twister');   % reproducibility

fcV_m   = fcData(:, idxMediation);
cogV_m  = slopeAll(:,1);
datV_m  = dynSlopes(:,3);
hippV_m = slope_L;

C_med = [fcMean.age, sexFC, fcMean.mean_FD, baselineEduc];

ok_m = ~isnan(fcV_m) & ~isnan(cogV_m) & ~isnan(datV_m) & ~isnan(hippV_m) & all(~isnan(C_med),2);
n_med = sum(ok_m);
fprintf('  n (complete on FC + cog + DAT + hipp + M1 covs): %d\n', n_med);

if n_med < 30
    fprintf('  [warn] n too small for stable mediation bootstrap — reporting point estimates only\n');
    nBoot = 0;
end

FC_m   = fcV_m(ok_m);
COG_m  = cogV_m(ok_m);
DAT_m  = datV_m(ok_m);
HIPP_m = hippV_m(ok_m);
Cm     = C_med(ok_m,:);
addIntercept = @(X) [ones(size(X,1),1), X];

% Point estimates
% a paths: mediator ~ FC + C
X_a  = addIntercept([FC_m, Cm]);
ba1  = X_a \ DAT_m;   a1 = ba1(2);
ba2  = X_a \ HIPP_m;  a2 = ba2(2);
% b paths and c': cog ~ FC + DAT + HIPP + C
X_b  = addIntercept([FC_m, DAT_m, HIPP_m, Cm]);
bb   = X_b \ COG_m;
cPrime = bb(2);
b1_m   = bb(3);
b2_m   = bb(4);
% Total effect c: cog ~ FC + C
bc   = X_a \ COG_m;   c_tot = bc(2);

ab1 = a1 * b1_m;
ab2 = a2 * b2_m;
abTot = ab1 + ab2;

% Bootstrap
if nBoot > 0
    ab1_b = nan(nBoot,1);
    ab2_b = nan(nBoot,1);
    abT_b = nan(nBoot,1);
    cP_b  = nan(nBoot,1);
    a1_b  = nan(nBoot,1);
    a2_b  = nan(nBoot,1);
    b1_b  = nan(nBoot,1);
    b2_b  = nan(nBoot,1);
    for bi = 1:nBoot
        idx = randi(n_med, n_med, 1);
        FCb = FC_m(idx); COGb = COG_m(idx);
        DATb = DAT_m(idx); HIPPb = HIPP_m(idx);
        Cb   = Cm(idx,:);
        Xab = addIntercept([FCb, Cb]);
        a1b = Xab \ DATb;   a1_b(bi) = a1b(2);
        a2b = Xab \ HIPPb;  a2_b(bi) = a2b(2);
        Xbb = addIntercept([FCb, DATb, HIPPb, Cb]);
        bbb = Xbb \ COGb;   cP_b(bi) = bbb(2); b1_b(bi) = bbb(3); b2_b(bi) = bbb(4);
        ab1_b(bi) = a1_b(bi) * b1_b(bi);
        ab2_b(bi) = a2_b(bi) * b2_b(bi);
        abT_b(bi) = ab1_b(bi) + ab2_b(bi);
    end
    CI_ab1 = prctile(ab1_b, [2.5 97.5]);
    CI_ab2 = prctile(ab2_b, [2.5 97.5]);
    CI_abT = prctile(abT_b, [2.5 97.5]);
    CI_cP  = prctile(cP_b,  [2.5 97.5]);
    sig_ab1 = (CI_ab1(1)>0) || (CI_ab1(2)<0);
    sig_ab2 = (CI_ab2(1)>0) || (CI_ab2(2)<0);
    sig_abT = (CI_abT(1)>0) || (CI_abT(2)<0);
else
    CI_ab1 = [NaN NaN]; CI_ab2 = [NaN NaN]; CI_abT = [NaN NaN]; CI_cP = [NaN NaN];
    sig_ab1 = false; sig_ab2 = false; sig_abT = false;
end

propMed = NaN;
if abs(c_tot) > eps, propMed = 100 * abTot / c_tot; end

fprintf('\n  Point estimates:\n');
fprintf('    a1  (FC -> DAT slope)        = %+.4f\n', a1);
fprintf('    a2  (FC -> L hipp slope)     = %+.4f\n', a2);
fprintf('    b1  (DAT -> cog | FC, hipp)  = %+.4f\n', b1_m);
fprintf('    b2  (hipp -> cog | FC, DAT)  = %+.4f\n', b2_m);
fprintf('    c   (total FC -> cog)        = %+.4f\n', c_tot);
fprintf('    c'' (direct FC -> cog)       = %+.4f\n', cPrime);

fprintf('\n  Indirect effects (bootstrap 95%% CI, %d reps):\n', nBoot);
if sig_ab1, s1 = 'SIG (CI excludes 0)'; else, s1 = 'n.s.'; end
if sig_ab2, s2 = 'SIG (CI excludes 0)'; else, s2 = 'n.s.'; end
if sig_abT, sT = 'SIG (CI excludes 0)'; else, sT = 'n.s.'; end
fprintf('    ab1 (via DAT slope)   = %+.4f   [%+.4f, %+.4f]   %s\n', ab1, CI_ab1(1), CI_ab1(2), s1);
fprintf('    ab2 (via L hipp)      = %+.4f   [%+.4f, %+.4f]   %s\n', ab2, CI_ab2(1), CI_ab2(2), s2);
fprintf('    ab total              = %+.4f   [%+.4f, %+.4f]   %s\n', abTot, CI_abT(1), CI_abT(2), sT);
fprintf('    Proportion mediated (ab_total / c) = %+.1f%%\n', propMed);

% Save mediation table
medRows = { ...
    {'a1_FC_to_DAT',      a1,     NaN, NaN};
    {'a2_FC_to_Hipp',     a2,     NaN, NaN};
    {'b1_DAT_to_Cog',     b1_m,   NaN, NaN};
    {'b2_Hipp_to_Cog',    b2_m,   NaN, NaN};
    {'c_total',           c_tot,  NaN, NaN};
    {'c_prime_direct',    cPrime, CI_cP(1), CI_cP(2)};
    {'ab1_indirect_DAT',  ab1,    CI_ab1(1), CI_ab1(2)};
    {'ab2_indirect_Hipp', ab2,    CI_ab2(1), CI_ab2(2)};
    {'ab_total_indirect', abTot,  CI_abT(1), CI_abT(2)};
    {'proportion_mediated_pct', propMed, NaN, NaN}};
writetable(cell2table(vertcat(medRows{:}), ...
    'VariableNames', {'Path','Estimate','CI_lo','CI_hi'}), ...
    fullfile(outputDir, 'Table_ParallelMediation_FC_to_Cog.csv'));
fprintf('\n  Saved: Table_ParallelMediation_FC_to_Cog.csv\n');

%% =========================================================================
%  17. DIRECT INDEPENDENCE TESTS — FC → cog, controlling for mediators
%      Does baseline FC retain prognostic value beyond the downstream
%      markers? Four nested models, same n_med sample (intersection of
%      FC + cog + DAT + hipp + M1 covariates):
%        (A) Base:          cog ~ FC + C
%        (B) + DAT slope:   cog ~ FC + DAT + C
%        (C) + hipp slope:  cog ~ FC + HIPP + C
%        (D) + both:        cog ~ FC + DAT + HIPP + C     [= c' from §16]
%      Compare beta_FC across the four fits — if FC effect survives (D),
%      FC has residual prognostic value beyond both mediators.
% =========================================================================
fprintf('\n[17] Direct independence: FC → cog, controlling for mediator(s)...\n');

indepSpecs = { ...
    'A_Base',          []; ...
    'B_plusDAT',       'DAT'; ...
    'C_plusHipp',      'HIPP'; ...
    'D_plusBoth',      'BOTH'};

indepRows = {};
for mi = 1:size(indepSpecs, 1)
    model_name = indepSpecs{mi, 1};
    med_type   = indepSpecs{mi, 2};
    if isempty(med_type)
        Tlm = table(COG_m, FC_m, Cm(:,1), Cm(:,2), Cm(:,3), Cm(:,4), ...
            'VariableNames', {'COG','FC','age','sex','FD','EDUC'});
        mdlFmla = 'COG ~ FC + age + sex + FD + EDUC';
    elseif strcmp(med_type,'DAT')
        Tlm = table(COG_m, FC_m, DAT_m, Cm(:,1), Cm(:,2), Cm(:,3), Cm(:,4), ...
            'VariableNames', {'COG','FC','DAT','age','sex','FD','EDUC'});
        mdlFmla = 'COG ~ FC + DAT + age + sex + FD + EDUC';
    elseif strcmp(med_type,'HIPP')
        Tlm = table(COG_m, FC_m, HIPP_m, Cm(:,1), Cm(:,2), Cm(:,3), Cm(:,4), ...
            'VariableNames', {'COG','FC','HIPP','age','sex','FD','EDUC'});
        mdlFmla = 'COG ~ FC + HIPP + age + sex + FD + EDUC';
    else  % BOTH
        Tlm = table(COG_m, FC_m, DAT_m, HIPP_m, Cm(:,1), Cm(:,2), Cm(:,3), Cm(:,4), ...
            'VariableNames', {'COG','FC','DAT','HIPP','age','sex','FD','EDUC'});
        mdlFmla = 'COG ~ FC + DAT + HIPP + age + sex + FD + EDUC';
    end
    ok_tbl = ~any(isnan(Tlm{:,:}), 2);
    mdl = fitlm(Tlm(ok_tbl,:), mdlFmla);

    coef  = mdl.Coefficients;
    b_FC  = coef{'FC','Estimate'};
    se_FC = coef{'FC','SE'};
    t_FC  = coef{'FC','tStat'};
    p_FC  = coef{'FC','pValue'};
    df    = mdl.DFE;
    n_fit = mdl.NumObservations;
    ci_lo = b_FC - tinv(0.975, df)*se_FC;
    ci_hi = b_FC + tinv(0.975, df)*se_FC;

    fprintf('\n  [%s]   n = %d, df = %d\n', model_name, n_fit, df);
    fprintf('    beta_FC = %+9.4f   SE = %.4f   95%% CI [%+.4f, %+.4f]\n', b_FC, se_FC, ci_lo, ci_hi);
    fprintf('    t(%d) = %+.3f   p = %.4f   R^2_adj = %.3f\n', df, t_FC, p_FC, mdl.Rsquared.Adjusted);

    indepRows{end+1} = {model_name, n_fit, df, b_FC, se_FC, ci_lo, ci_hi, t_FC, p_FC, mdl.Rsquared.Adjusted}; %#ok
end

writetable(cell2table(vertcat(indepRows{:}), ...
    'VariableNames',{'Model','n','df','beta_FC','SE','CI_lo','CI_hi','t','p_FC','R2_adj'}), ...
    fullfile(outputDir,'Table_FC_Cog_Independence.csv'));
fprintf('\n  Saved: Table_FC_Cog_Independence.csv\n');

%% =========================================================================
%  18. SERIAL MEDIATION (EXPLORATORY) — FC → DAT → hipp → cog (+ reverse)
%      Three specific indirects per order:
%        IND1: FC → M1 → cog        (simple)
%        IND2: FC → M2 → cog        (simple, bypass M1)
%        IND3: FC → M1 → M2 → cog   (serial)
%      Total indirect = IND1 + IND2 + IND3.
%      Forward order:  M1 = DAT,  M2 = hipp
%      Reverse order:  M1 = hipp, M2 = DAT   (sensitivity)
%      Covariates: M1 primary = age + sex + FD + EDUC.
%      Preacher & Hayes percentile bootstrap, 5,000 reps.
%      CAVEAT: with 2-timepoint slopes we cannot establish true temporal
%      order — "serial" here is statistical, not causal.
% =========================================================================
fprintf('\n[18] Serial mediation (exploratory): FC → M1 → M2 → cog...\n');

if n_med < 30
    fprintf('  [warn] n = %d too small for serial mediation bootstrap — skipping\n', n_med);
else
    rng(42, 'twister');

    orderLbls = {'forward_FC_DAT_Hipp_Cog','reverse_FC_Hipp_DAT_Cog'};
    serialRows = {};

    for oi = 1:2
        if oi == 1
            M1name = 'DAT'; M2name = 'Hipp';
            M1v = DAT_m; M2v = HIPP_m;
        else
            M1name = 'Hipp'; M2name = 'DAT';
            M1v = HIPP_m; M2v = DAT_m;
        end
        fprintf('\n  === [%s: M1 = %s, M2 = %s] ===\n', orderLbls{oi}, M1name, M2name);

        % Point estimates
        X_FCC  = addIntercept([FC_m, Cm]);
        X_FCM1C = addIntercept([FC_m, M1v, Cm]);
        X_full  = addIntercept([FC_m, M1v, M2v, Cm]);
        a1  = (X_FCC  \ M1v); a1 = a1(2);
        bM2 = (X_FCM1C \ M2v); a2 = bM2(2); d21 = bM2(3);
        bY  = (X_full  \ COG_m); cPrime = bY(2); b1 = bY(3); b2 = bY(4);
        bC  = (X_FCC  \ COG_m); cTot = bC(2);

        IND1 = a1 * b1;            % FC → M1 → cog
        IND2 = a2 * b2;            % FC → M2 → cog (bypass M1)
        IND3 = a1 * d21 * b2;      % FC → M1 → M2 → cog (serial)
        INDT = IND1 + IND2 + IND3;

        % Bootstrap
        IND1_b = nan(nBoot,1); IND2_b = nan(nBoot,1);
        IND3_b = nan(nBoot,1); INDT_b = nan(nBoot,1);
        for bi = 1:nBoot
            idx = randi(n_med, n_med, 1);
            FCb = FC_m(idx); COGb = COG_m(idx);
            M1b = M1v(idx);  M2b = M2v(idx);
            Cb  = Cm(idx,:);
            Xa1 = addIntercept([FCb, Cb]);
            t1  = Xa1 \ M1b;             a1b = t1(2);
            Xa2 = addIntercept([FCb, M1b, Cb]);
            t2  = Xa2 \ M2b;             a2b = t2(2); d21b = t2(3);
            Xb  = addIntercept([FCb, M1b, M2b, Cb]);
            t3  = Xb  \ COGb;            b1b = t3(3); b2b = t3(4);
            IND1_b(bi) = a1b * b1b;
            IND2_b(bi) = a2b * b2b;
            IND3_b(bi) = a1b * d21b * b2b;
            INDT_b(bi) = IND1_b(bi) + IND2_b(bi) + IND3_b(bi);
        end
        CI_IND1 = prctile(IND1_b, [2.5 97.5]);
        CI_IND2 = prctile(IND2_b, [2.5 97.5]);
        CI_IND3 = prctile(IND3_b, [2.5 97.5]);
        CI_INDT = prctile(INDT_b, [2.5 97.5]);

        sigFn = @(ci) (ci(1)>0) || (ci(2)<0);
        if sigFn(CI_IND1), s1 = 'SIG'; else, s1 = 'n.s.'; end
        if sigFn(CI_IND2), s2 = 'SIG'; else, s2 = 'n.s.'; end
        if sigFn(CI_IND3), s3 = 'SIG'; else, s3 = 'n.s.'; end
        if sigFn(CI_INDT), sT = 'SIG'; else, sT = 'n.s.'; end

        fprintf('    a1  (FC -> %s)             = %+.4f\n', M1name, a1);
        fprintf('    a2  (FC -> %s | %s)         = %+.4f\n', M2name, M1name, a2);
        fprintf('    d21 (%s -> %s)              = %+.4f\n', M1name, M2name, d21);
        fprintf('    b1  (%s -> cog)            = %+.4f\n', M1name, b1);
        fprintf('    b2  (%s -> cog)            = %+.4f\n', M2name, b2);
        fprintf('    c   (total FC -> cog)      = %+.4f\n', cTot);
        fprintf('    c'' (direct FC -> cog)     = %+.4f\n', cPrime);
        fprintf('    IND1 (FC -> %s -> cog)         = %+.4f  [%+.4f, %+.4f]   %s\n', M1name, IND1, CI_IND1(1), CI_IND1(2), s1);
        fprintf('    IND2 (FC -> %s -> cog)         = %+.4f  [%+.4f, %+.4f]   %s\n', M2name, IND2, CI_IND2(1), CI_IND2(2), s2);
        fprintf('    IND3 (FC -> %s -> %s -> cog)   = %+.4f  [%+.4f, %+.4f]   %s\n', M1name, M2name, IND3, CI_IND3(1), CI_IND3(2), s3);
        fprintf('    Total indirect                  = %+.4f  [%+.4f, %+.4f]   %s\n', INDT, CI_INDT(1), CI_INDT(2), sT);
        if abs(cTot) > eps
            fprintf('    Proportion mediated (IND_total / c) = %+.1f%%\n', 100*INDT/cTot);
        end

        % Append to save table
        serialRows{end+1} = {orderLbls{oi}, sprintf('a1_FC_to_%s', M1name),               a1,     NaN, NaN}; %#ok
        serialRows{end+1} = {orderLbls{oi}, sprintf('a2_FC_to_%s_given_%s', M2name, M1name), a2,  NaN, NaN}; %#ok
        serialRows{end+1} = {orderLbls{oi}, sprintf('d21_%s_to_%s', M1name, M2name),      d21,    NaN, NaN}; %#ok
        serialRows{end+1} = {orderLbls{oi}, sprintf('b1_%s_to_Cog', M1name),              b1,     NaN, NaN}; %#ok
        serialRows{end+1} = {orderLbls{oi}, sprintf('b2_%s_to_Cog', M2name),              b2,     NaN, NaN}; %#ok
        serialRows{end+1} = {orderLbls{oi}, 'c_total',                                    cTot,   NaN, NaN}; %#ok
        serialRows{end+1} = {orderLbls{oi}, 'c_prime_direct',                             cPrime, NaN, NaN}; %#ok
        serialRows{end+1} = {orderLbls{oi}, sprintf('IND1_FC_%s_Cog', M1name),            IND1,   CI_IND1(1), CI_IND1(2)}; %#ok
        serialRows{end+1} = {orderLbls{oi}, sprintf('IND2_FC_%s_Cog', M2name),            IND2,   CI_IND2(1), CI_IND2(2)}; %#ok
        serialRows{end+1} = {orderLbls{oi}, sprintf('IND3_FC_%s_%s_Cog', M1name, M2name), IND3,   CI_IND3(1), CI_IND3(2)}; %#ok
        serialRows{end+1} = {orderLbls{oi}, 'IND_total',                                  INDT,   CI_INDT(1), CI_INDT(2)}; %#ok
    end

    writetable(cell2table(vertcat(serialRows{:}), ...
        'VariableNames', {'Order','Path','Estimate','CI_lo','CI_hi'}), ...
        fullfile(outputDir, 'Table_SerialMediation_FC.csv'));
    fprintf('\n  Saved: Table_SerialMediation_FC.csv\n');
end

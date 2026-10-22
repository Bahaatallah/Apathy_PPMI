% Role: Connectivity-to-outcome chain, step 1 of 3 for Figure 4. Participant-level global cognitive
%     slopes (bifactor g; reduced LME BLUPs with an OLS fallback, plus a winsorised copy used as a
%     sensitivity), path a (apathy to component FC), path b (component FC to global cognitive slope,
%     Spearman partial correlation adjusted for mean FD and TIV, BH FDR across components) and the
%     apathy to FC to slope mediation with percentile bootstrap CIs and its sensitivity conditions.
%     Writes Table_Fig4_Sensitivity.csv, Table_Fig4_Mediation.csv and
%     Table_Fig4_Mediation_Sensitivity.csv to outputDir.
% Run: on the workspace left by behavioural/cohort_phenotyping.m, after setting outputDir (output
%     folder), EXT_FC_DIR (folder holding nbs_edge_strengths.csv and nbs_suprathreshold_fc.csv),
%     EXT_FEATURE_SET ('P305', primary sample, or 'X343', extended sample) and EXT_MEDIATOR ('diff',
%     the apathy-specific differential component). Step 2 is fc_structural_correlates.m (same
%     workspace), step 3 figure4_connectivity_outcomes.m, which draws the figure.
% Manuscript item: Figure 4a values; Results (connectivity and cognitive decline, mediation);
%     Supplementary Table S10.

fprintf('\n%s\n', repmat('=',1,72));
fprintf('  FIGURE 4: FC MEDIATES APATHY-RELATED COGNITIVE DECLINE\n');
fprintf('%s\n\n', repmat('=',1,72));

%% ========================================================================
%  0.  SETTINGS
% =========================================================================
if ~exist('outputDir','var') || isempty(outputDir), outputDir = pwd; end

fdThreshold = 0.3;   % FD threshold of the feature set

% --- FC source for path b (FC -> cognition) ---
%     true  = suprathreshold FC (mean of NBS-flagged edges only — cleaner signal)
%     false = mean FC (mean of ALL edges in connected component — broader)
%     Path a (apathy -> FC) always uses suprathreshold.
%     When true, only apathy components (which have supra columns) are in
%     the primary FDR family; diff/shared fall back to mean FC for the
%     exploratory heatmap (and remain outside the primary family).
useSupraForPathB = true;

fprintf('  FD threshold: %.1f\n', fdThreshold);
fprintf('  Path b FC source: %s\n\n', ternary(useSupraForPathB,'suprathreshold','mean'));

% feature CSVs from EXT_FC_DIR (set in the workspace); see header.
assert(exist('EXT_FC_DIR','var')==1 && isfolder(EXT_FC_DIR), 'EXT_FC_DIR must be set to the feature CSV folder');
if ~exist('EXT_FEATURE_SET','var') || isempty(EXT_FEATURE_SET), EXT_FEATURE_SET = 'P305'; end
if ~exist('EXT_MEDIATOR','var')    || isempty(EXT_MEDIATOR),    EXT_MEDIATOR    = 'diff'; end
assert(ismember(EXT_FEATURE_SET, {'P305','X343'}) && strcmp(EXT_MEDIATOR, 'diff'));
fprintf('  Feature set %s, mediator %s, dir %s\n', EXT_FEATURE_SET, EXT_MEDIATOR, EXT_FC_DIR);
basePath = EXT_FC_DIR;
fcPath_mean  = fullfile(basePath, 'nbs_edge_strengths.csv');
fcPath_supra = fullfile(basePath, 'nbs_suprathreshold_fc.csv');

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

% Resolve variable names (strict, no fallback). Global cognition uses bifactor_g.
cogVarG = 'bifactor_g';
assert(ismember(cogVarG, Tlong.Properties.VariableNames), ...
    'Tlong is missing column %s — bifactor factor scores must be merged before running Figure 4.', cogVarG);
timeVar=''; for c={'time','YEAR','year_num'}, if ismember(c{1},Tlong.Properties.VariableNames), timeVar=c{1}; break; end, end
assert(~isempty(timeVar), 'Tlong is missing a time column (time/YEAR/year_num).');

% Domain variables (strict)
cogDomains  = {'global'};
cogVarNames = {'bifactor_g'};
cogLabels   = {'Global Cognition'};
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
if fdThreshold == 0.3
    % ---- N = 305 FINAL PRIMARY (protocol-stratified, corrected FD, coverage
    %      rule: 15 dropout parcels blanked, 2 participants dropped).
    %      Components of the primary NBS run:
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
    % ---- N = 343 EXTENSION (protocol-heterogeneous sample):
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
end
nAllFC     = numel(meanFCvars);
nApFC      = numel(apathyIdx);
% --- supraIdx: components with suprathreshold columns (path a and its FDR family). Distinct
%     from apathyIdx, which marks the contrasts that are semantically apathy.
supraIdx = 1:numel(supraFCvars);
nSupraFC   = numel(supraIdx);
assert(numel(supraFCvars) == nSupraFC, ...
    'supraFCvars length must equal numel(supraIdx).');

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
%  2.  EXTRACT SLOPES (reduced LME BLUP with OLS fallback)
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

% --- Sensitivity-only: 1st/99th percentile winsorised copy for section 14 ------
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

% Participant-level global slopes (input of figure4_connectivity_outcomes.m)
writetable(table(fcPATNO, slopeAll(:,1), 'VariableNames', {'subject_id','slope_g'}), ...
    fullfile(outputDir, 'participant_slopes_g.csv'));

% --- Per-subject baseline clinical covariates (for sensitivity analysis) -
%     Pull baseline (min-time row per subject). Column names vary across
%     tables, so try common variants and report which was found.

% Define: display name -> candidate column names in Tlong / Tpd_lme
covCandsMap = { ...
    'duration',  {'duration','disease_duration','PD_duration','pd_duration'}; ...
    'LEDD',      {'ledd','LEDD','ledd_baseline'}; ...
    'UPDRS3',    {'updrs3_score','UPDRS3','updrs3','mds_updrs3'}; ...
    'QUIP',      {'quip','QUIP','quip_score','QUIP_score','quip_total'}; ...
    'depfactor', {'depfactor','depression_factor','depfctr','GDS_depfactor','dep_score'}; ...
    'educ',      {'educ','EDUCYRS','education','education_years','EDUC'}};

tableCandidates = {'Tlong'};
if exist('Tpd_lme','var'), tableCandidates{end+1} = 'Tpd_lme'; end
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

% Primary partial-correlation tests (sections 3 and 4): Spearman partial r with FD + TIV as
% covariates, the covariate set of the section 6 mediation. Section 14 adds the sensitivity
% covariate sets.
sexNum = sexFC;                                           % M=1 / F=0
primaryCovs   = [fcMean.mean_FD, fcMean.TIV];
primaryCovLbl = {'FD','TIV'};
primaryMethod = 'Spearman';   % Capitalised to match MATLAB partialcorr's 'Type' values
nPrimaryComplete = sum(all(~isnan(primaryCovs),2));
fprintf('\n  Primary: %s partial r | %s\n', primaryMethod, strjoin(primaryCovLbl,', '));
fprintf('  Subjects complete on primary covariates: %d / %d\n', nPrimaryComplete, nFC);

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
%     components using suprathreshold FC.
if useSupraForPathB
    primaryIdx = supraIdx;
else
    primaryIdx = 1:nAllFC;
end

% Pre-allocate compR with all fields, length = nAllFC (indexed by fi)
emptyEntry = struct('r',NaN,'p',NaN,'rP',NaN,'pP',NaN,'pFDR',NaN, ...
                    'n',NaN,'nP',NaN,'x',[],'y',[],'mdl',[], ...
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
%  6.  MEDIATION — Apathy → FC → Global Cog Decline
%      Run for every primary component; idxMediation flags the one used
%      for Panel B visualization. Covariates = FD + TIV only.
% =========================================================================
fprintf('\n[6] Mediation: Apathy -> FC -> Global Decline (bootstrap %d reps)...\n', nBoot);

medResults = repmat(struct('fi',NaN,'label','','n',NaN,'a',NaN,'a_p',NaN, ...
    'b',NaN,'b_p',NaN,'cp',NaN,'cp_p',NaN,'c',NaN,'c_p',NaN, ...
    'ab',NaN,'ab_ci',[NaN NaN],'ab_sig',false,'prop',NaN, ...
    'ab_boot',[]), ...
    1, nAllFC);

medIdxList = unique([primaryIdx(:); idxMediation]);  % always include idxMediation

for fi = medIdxList(:)'
    medFC = fcData(:, fi);
    medSlope = slopeAll(:,1);   % global
    medAp = fcMean.apathy_composite;
    % Covariates = FD + TIV only; no per-outcome baseline g, because the BLUP slope already
    % absorbs the level-side covariates.
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
%     Re-run the section 6 mediation for idxMediation across covariate sets and
%     transformations. The 'Baseline' row is intentionally identical to
%     the section 6 primary mediation: same outcome (slopeAll(:,1)), same
%     predictor (apathy_composite), same mediator (fcData(:,idxMediation)),
%     same covariates (FD + TIV), same okM filter (no per-outcome
%     baselineG). Every other row perturbs one or
%     more of {slope transform, covariate set, rank vs linear} relative
%     to that baseline.
fprintf('\n  Mediation sensitivity (idxMediation = fi=%d, %s):\n', ...
    idxMediation, labelMed);
fprintf('  %-22s  %-6s  %10s  %-23s  %s\n', ...
    'Condition', 'N', 'ab', 'Percentile 95% CI', 'prop');
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
    'ExtCovars',        false, 'ext',        false};

medSensRows = {};
rng(42);   % bootstrap stream of the sensitivity conditions, independent of the preceding sections
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

    fprintf('  %-22s  %-6d  %+9.4f  [%+.4f, %+.4f] %-3s  %.1f%%\n', ...
        lbl, n_, ab_, ...
        ci_(1), ci_(2), ternary(sig_,'SIG','n.s.'), ...
        prop_);
    medSensRows{end+1} = {lbl, n_, a_, b_, cp_, c_, ab_, ...
        ci_(1), ci_(2), double(sig_), ...
        prop_}; %#ok
end
if ~isempty(medSensRows)
    writetable(cell2table(vertcat(medSensRows{:}), 'VariableNames', { ...
        'Condition','N','a','b','cp','c','ab', ...
        'Pct_CI_lo','Pct_CI_hi','Pct_Sig', ...
        'PctMediated'}), ...
        fullfile(outputDir,'Table_Fig4_Mediation_Sensitivity.csv'));
    fprintf('  -> Table_Fig4_Mediation_Sensitivity.csv saved.\n');
end

%% ========================================================================
%  LOCAL FUNCTIONS
% =========================================================================

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

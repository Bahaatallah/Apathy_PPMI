% Role: Primary cognition LMEs (g and MoCA), sensitivity models (winsorised, Cook's distance,
%     minimal, baseline-adjusted, pattern mixture), sample-size table and Figure 2.
% Run: on the workspace left by 01_clinical_cognition/Publicaiton_analysis.m, loaded in a fresh MATLAB
%     session (this script standardises columns of lmeT in place), with DESIGN_OVERRIDE absent
%     (primary design); outputs go to outputDir.
% Manuscript item: Figure 2; Results (cognitive decline); Supplementary Tables S3 Panel C, S4, S5,
%     S6.
%% Primary cognition LMEs (g and MoCA), sensitivity models and manuscript Figure 2
%  (Fig2_Primary_Trajectory_Forest). Okabe-Ito palette shared with the other figure scripts.

% --- Palette (Okabe-Ito); overrides workspace COL.
COL.neither  = [0.00 0.60 0.50];   % Okabe-Ito bluish green
COL.apathy   = [0.90 0.60 0.00];   % Okabe-Ito orange
COL.icb      = [0.00 0.45 0.70];   % Okabe-Ito blue
COL.overlap  = [0.45 0.20 0.50];   % deep purple
COL.groups   = [COL.neither; COL.apathy; COL.icb; COL.overlap];
COL.bg       = [1 1 1];
COL.grid     = [0.85 0.85 0.85];
COL.txt      = [0.15 0.15 0.15];
COL.accent   = [0.10 0.10 0.10];

%% ========================================================================
%  SUPPLEMENTARY FIGURE 2  (formerly Figure 2)
%  Combined CFA Cognitive Trajectories + LME Forest Plots
%  Top row  (A–D): Mean ± SEM trajectories for each cognitive domain
%                   across four motivational phenotype groups.
%  Bottom row (E–H): Forest plots of key LME fixed effects (adjusted).
%
%  The PRIMARY Figure 2 (built later in this file) shows: bifactor_g
%  trajectory + paired bifactor_g/MoCA forest with apathy:time vs icb:time
%  bracket + 3-domain heatmap.
%
%  Prerequisites:
%    - CFA factor scores merged into Tpd_lme (Section 23b)
%    - lmeT built with: apathy, icb, time, covariates, PATNO, bifactor_*
%    - okBase logical index
%    - Helpers: safeCol, exportFigure, formatAxis, plotTrajectory,
%               coeffToTable, pToStars, COL struct, groupLabelsShort
% =========================================================================
fprintf('[Supp Fig2] Combined CFA trajectories + LME forest plots...\n');

% -----------------------------wh--------------------------------------------
% 0.  SETTINGS
% -------------------------------------------------------------------------
% Global cognition outcome — flip to 'cfa_g' to use the CFA single-factor
% model instead of the bifactor g. Affects panel A (trajectory), panel B
% (forest, legend / title text), and panel C (covariate the sub-domain
% LMEs partial out).
globalVar      = 'bifactor_g';     % 'bifactor_g' | 'cfa_g'
globalVarLabel = 'factor g';     % display name (used in panel B)
% Override when cfa_g is selected, so plot text matches the model
if strcmp(globalVar, 'cfa_g'), globalVarLabel = 'CFA g'; end

% Panel C sub-domain LMEs: include globalVar as a covariate (true) so the
% apathy / ICB / time effects are conditional on global cognition, OR fit
% without it (false) so they reflect total effects on each sub-domain.
controlForGlobal = false;     % true | false

% cogVars   = {globalVar,'bifactor_att_exec','bifactor_memory','bifactor_visuo'};
cogVars = [{globalVar}, {'obs_att_exec','obs_memory','obs_visuo'}];
% cogVars = [{globalVar}, {'cfa_att_exec','cfa_memory','cfa_visuo'}];

cogLabels = {'Global Cognition (g)','Executive & Attention','Memory','Visuospatial'};

% -------------------------------------------------------------------------
% Re-residualise sub-domain observed scores against the CURRENT globalVar.
% Upstream obs_*_resid columns (built in Publicaiton_analysis.m) are
% residualised against whichever g was active when that script ran. Flip
% globalVar here and those columns become inconsistent with the selected
% global outcome — fix that by overwriting in-script. Both lmeT and
% Tpd_lme are updated so trajectories (Tpd_lme) and LMEs (lmeT) share
% one definition. Simple row-level OLS residualisation: y - (a + b * g).
% -------------------------------------------------------------------------
subdomainRaw = {'obs_att_exec','obs_memory','obs_visuo'};
residualiseTbl = @(T) T;        % placeholder; we update inline below

for i = 1:numel(subdomainRaw)
    raw   = subdomainRaw{i};
    resid = [raw '_resid'];

    % --- lmeT ---
    if all(ismember({raw, globalVar}, lmeT.Properties.VariableNames))
        v = ~isnan(lmeT.(raw)) & ~isnan(lmeT.(globalVar));
        if sum(v) >= 5
            X = [ones(sum(v),1), lmeT.(globalVar)(v)];
            y = lmeT.(raw)(v);
            b = X \ y;
            rcol = nan(height(lmeT), 1);
            rcol(v) = y - X * b;
            lmeT.(resid) = rcol;
        end
    end

    % --- Tpd_lme (used by trajectory plots) ---
    if all(ismember({raw, globalVar}, Tpd_lme.Properties.VariableNames))
        v = ~isnan(Tpd_lme.(raw)) & ~isnan(Tpd_lme.(globalVar));
        if sum(v) >= 5
            X = [ones(sum(v),1), Tpd_lme.(globalVar)(v)];
            y = Tpd_lme.(raw)(v);
            b = X \ y;
            rcol = nan(height(Tpd_lme), 1);
            rcol(v) = y - X * b;
            Tpd_lme.(resid) = rcol;
        end
    end
end
fprintf('  Sub-domain resid columns rebuilt against %s.\n', globalVar);
years     = 0:15;
nGrps     = 4;
minNplot  = 5;   % minimum n per time-point to plot (avoids misleading SEM)

% Top-row visualisation mode
%   'group4' – 4-way phenotype (Neither | Apathy | ICB | Both)
%   'apathy' – 2-way split by clinical apathy flag  (Low vs High)
%   'icb'    – 2-way split by clinical ICB flag     (Low vs High)
% Splits reuse the flags already encoded in group4 so the cut-points
% match the ones used for the group4 analysis.
vizMode = 'group4';

% --- Forest-plot palette (distinct from group colours) ---
fpCol.main   = [0.25 0.30 0.38];   % dark slate  – main effects
fpCol.two    = [0.17 0.50 0.55];   % teal        – 2-way interactions
fpCol.three  = [0.72 0.38 0.20];   % copper      – 3-way interaction
fpCol.ns     = [0.70 0.70 0.70];   % grey        – non-significant

% Explicit LME terms of interest (order = plotting order, top to bottom)
termsOfInterest = {'apathy','icb','time', ...
                   'apathy:icb','apathy:time','icb:time', ...
                   'apathy:icb:time'};
termDisplayNames = {'Apathy','ICB','Time', ...
                    'Apathy \times ICB','Apathy \times Time','ICB \times Time', ...
                    'Apathy \times ICB \times Time'};

% -------------------------------------------------------------------------
% HELPER: robust term matching
% -------------------------------------------------------------------------
% coeffToTable may return names with spaces around ':' or use '*'.
% This normalises both sides before comparing.
normTerm = @(s) strjoin(sort(strsplit(strtrim(regexprep(lower(s), '\s*[:\*]\s*', ':')), ':')), ':');

matchTerm = @(feNames, target) find(cellfun( ...
    @(x) strcmp(normTerm(x), normTerm(target)), feNames), 1);

% -------------------------------------------------------------------------
% 1.  CREATE FIGURE
% -------------------------------------------------------------------------
figSupp2 = figure('Units','centimeters','Position',[1 1 30 12],'Color',COL.bg);
set(figSupp2, 'DefaultAxesFontName','Helvetica','DefaultTextFontName','Helvetica');

nCols = numel(cogVars);
tLayout = tiledlayout(2, nCols, 'TileSpacing','compact', 'Padding','compact');

% -------------------------------------------------------------------------
% 2.  TOP ROW — LONGITUDINAL TRAJECTORIES  (panels A–D)
% -------------------------------------------------------------------------
TpdSub = Tpd_lme;

% Define smoothing factor for plotCPTraj's Gaussian smooth
sigmaVal = 1.0;

% Resolve vizMode into grouping vector / colours / legend labels.
% apathyFlag <=> group4 ∈ {2,4};  icbFlag <=> group4 ∈ {3,4}
switch vizMode
    case 'group4'
        vizGrp    = TpdSub.group4;
        vizLevels = (1:4)';
        vizColors = COL.groups;
        vizLabels = groupLabelsShort;
    case 'apathy'
        vizGrp = nan(height(TpdSub), 1);
        vizGrp(TpdSub.group4 == 1 | TpdSub.group4 == 3) = 1;   % Low apathy
        vizGrp(TpdSub.group4 == 2 | TpdSub.group4 == 4) = 2;   % High apathy
        vizLevels = [1; 2];
        vizColors = [COL.groups(1,:); COL.apathy];
        vizLabels = {'Low apathy','High apathy'};
    case 'icb'
        vizGrp = nan(height(TpdSub), 1);
        vizGrp(TpdSub.group4 == 1 | TpdSub.group4 == 2) = 1;   % Low ICB
        vizGrp(TpdSub.group4 == 3 | TpdSub.group4 == 4) = 2;   % High ICB
        vizLevels = [1; 2];
        vizColors = [COL.groups(1,:); COL.icb];
        vizLabels = {'Low ICB','High ICB'};
    otherwise
        error('Figure2:badVizMode', ...
            'vizMode must be ''group4'', ''apathy'', or ''icb''.');
end
nVizGrps = numel(vizLevels);

for cv = 1:nCols
    ax = nexttile(cv);
    set(ax, 'Color', COL.bg); hold on;

    ph = gobjects(nVizGrps, 1);

    for g = 1:nVizGrps
        rowFilter = (vizGrp == vizLevels(g));
        ph(g) = plotCPTraj(ax, TpdSub, rowFilter, 'YEAR', cogVars{cv}, ...
                           vizColors(g,:), [], sigmaVal, '-');
    end
    
    % Zero reference line
    yline(ax, 0, '--', 'Color', [0.65 0.65 0.65], 'LineWidth', 0.6);
    
    % Panel letter + domain title
    panelLetter = char('A' + cv - 1);
    title(ax, sprintf('\\bf%s\\rm  %s', panelLetter, cogLabels{cv}), ...
        'FontSize', 8, 'Color', COL.txt);
        
    % Axes
    xlim([-0.3  max(years)+0.3]);
    set(ax, 'FontSize', 6.5, 'Box', 'off', 'TickDir', 'out', ...
        'XTick',[0, 5, 10, 15], 'XGrid', 'off', 'YGrid', 'off');
    xlabel('Year', 'FontSize', 7);
    
    if cv == 1
        ylabel('Factor score (z)', 'FontSize', 7);
    end
    
    % Tidy limits
    yl = ylim;
    ylim([yl(1) - 0.02*diff(yl),  yl(2) + 0.06*diff(yl)]);
    
    % Legend (first panel only)
    if cv == 1
        validPh = ph(isgraphics(ph));
        validLg = vizLabels(isgraphics(ph));
        lg = legend(validPh, validLg, ...
            'Location','southwest','FontSize',5.5,'Box','off');
        lg.ItemTokenSize = [12 8];
    end
    
    formatAxis(ax, COL);
    hold off;
end
% -------------------------------------------------------------------------
% 3.  Z-SCORE CONTINUOUS PREDICTORS
% -------------------------------------------------------------------------
% Rationale: In a three-way interaction model (apathy * icb * time),
%   main effects are conditional on the other predictors = 0.  Without
%   centring, "apathy = 0" may be extrapolated or uninterpretable.
%   Grand-mean z-scoring places zero at the sample mean, so main effects
%   represent the effect at an average level of the other predictors.
%   It also puts apathy and ICB on a common SD scale, making betas
%   directly comparable on the forest plot.
%
%   Time is left in its natural metric (years) so that betas read as
%   "change per year", which is the clinical convention.
%
%   Z-scoring uses the FULL lmeT sample (all available observations)
%   so that the reference point is stable regardless of which outcome
%   has missing data.  For time-invariant predictors (apathy, icb, educ,
%   duration) this is equivalent to z-scoring from the baseline distribution.

zScore = @(x) (x - mean(x, 'omitnan')) ./ std(x, 0, 'omitnan');

% Key predictors
lmeT.apathy    = zScore(lmeT.apathy);
lmeT.icb       = zScore(lmeT.icb);

% Continuous covariates
lmeT.educ      = zScore(lmeT.educ);
lmeT.nmiss     = zScore(lmeT.nmiss);
lmeT.duration  = zScore(lmeT.duration);
lmeT.ledd      = zScore(lmeT.ledd);
lmeT.updrs3_score        = zScore(lmeT.updrs3_score);
lmeT.depfactor = zScore(lmeT.depfactor);
lmeT.age = zScore(lmeT.age);
lmeT.ageSq = zScore(lmeT.age.^2);
lmeT.ageCub = zScore(lmeT.age.^3);
lmeT.hy = zScore(lmeT.hy);


% Time stays in years (not z-scored)
fprintf('  Predictors z-scored (apathy, icb, educ, nmiss, duration, ledd, updrs3_score, depfactor, age).\n');
fprintf('  Time retained in years.\n');
% design switch. DESIGN_OVERRIDE = 'D2' (set in the calling workspace) adds the
% baseline apathy and ICB scores as main-effect covariates to every apathy*icb*time formula
% (the runner has already replaced lmeT.apathy / lmeT.icb by the per-visit scores and kept the
% baseline scores as lmeT.apathy_bl / lmeT.icb_bl). Default '' reproduces the published models.
if exist('DESIGN_OVERRIDE', 'var') && strcmpi(DESIGN_OVERRIDE, 'D2'), BLT = ' + apathy_bl + icb_bl'; else, BLT = ''; end
fprintf('  Design: %s\n', ternaryStr(isempty(BLT), 'published (baseline apathy and ICB as the exposures)', ['D2: per-visit exposures' BLT]));

% -------------------------------------------------------------------------
% 4.  RUN ADJUSTED LMEs  (store everything for forest plot + Supp Table 2)
% -------------------------------------------------------------------------
allFE        = cell(1, nCols);    % fixed-effect tables (one per domain)
allModelInfo = cell(1, nCols);    % struct per domain for Supp Table 2
tableNums    = {'Table8','Table9','Table10','Table11'};

for ci = 1:nCols
    outVar = cogVars{ci};

    % --- Filter to complete cases for adjusted model ---
    % NOTE: every fixed-effect predictor (incl. the apathy/icb/time terms
    % that drive the 3-way interaction) and every covariate must be NaN-
    % screened here so that height(lmeD) == mdl.NumObservations downstream.
    % If we miss one, fitlme silently drops those rows but
    % designMatrix(mdl,'Fixed') still returns them as NaN — which then
    % poisons the QR step in §S1 Cook's D (every leverage → NaN → 0 → no
    % observation ever flagged).
    okAdj = okBase & ~isnan(lmeT.(outVar)) & ...
            ~isnan(lmeT.apathy) & ~isnan(lmeT.icb) & ~isnan(lmeT.time) & ...
            ~isnan(lmeT.age) & ~isnan(lmeT.educ) & ~isnan(lmeT.nmiss) & ...
            ~isnan(lmeT.duration) & ~isnan(lmeT.ledd) & ~isnan(lmeT.updrs3_score) & ...
            ~isundefined(lmeT.sex) & ~isnan(lmeT.depfactor);
    lmeD = lmeT(okAdj, :);

    nObs = height(lmeD);
    nPat = numel(unique(lmeD.PATNO));

    try
        % === REML fit for final coefficient estimates ===
        fmla = sprintf(['%s ~ apathy*icb*time' BLT ' + age + sex + educ + nmiss + ' ...
                        'duration + ledd + updrs3_score + depfactor + (1 + time | PATNO) + (1| time)'], outVar);
        mdl  = fitlme(lmeD, fmla, 'FitMethod', 'REML');

        if mdl.NumObservations ~= nObs
            fprintf('  %s — %d obs (filter), %d used by LME, %d patients  [WARN: %d silent NaN drops]\n', ...
                cogLabels{ci}, nObs, mdl.NumObservations, nPat, nObs - mdl.NumObservations);
        else
            fprintf('  %s — %d obs, %d patients\n', cogLabels{ci}, nObs, nPat);
        end

        % Extract coefficients with model-based CIs
        fe        = coeffToTable(mdl);
        modelCI   = coefCI(mdl);
        fe.CI_lo  = modelCI(:, 1);
        fe.CI_hi  = modelCI(:, 2);
        fe.Domain = repmat(cogLabels(ci), height(fe), 1);

        % --- Diagnostic: check which terms of interest are present ---
        for ti = 1:numel(termsOfInterest)
            r = matchTerm(fe.Name, termsOfInterest{ti});
            if isempty(r)
                fprintf('    WARNING: term "%s" not found in %s model.\n', ...
                    termsOfInterest{ti}, cogLabels{ci});
                fprintf('    Available names: %s\n', strjoin(fe.Name, ', '));
            end
        end

        allFE{ci} = fe;

        % Save per-domain CSV
        writetable(fe, fullfile(outputDir, ...
            sprintf('%s_LME_%s_Adjusted.csv', tableNums{ci}, outVar)));

        % === LRT for 3-way interaction (both models must use ML) ===
        fmla_2w = sprintf(['%s ~ apathy*icb + apathy*time + icb*time + ' ...
            'age + sex + educ + nmiss + duration + ledd + updrs3_score + depfactor + ' ...
            '(1 + time | PATNO)'], outVar);
        mdl_2w = fitlme(lmeD, fmla_2w, 'FitMethod', 'ML');
        mdl_3w = fitlme(lmeD, fmla,    'FitMethod', 'ML');
        comp   = compare(mdl_2w, mdl_3w);
        lrt_chi2 = comp.LRStat(2);
        lrt_df   = comp.deltaDF(2);
        lrt_p    = comp.pValue(2);
        fprintf('    LRT 3-way: Chi2=%.2f  df=%d  p=%.4f %s\n', ...
            lrt_chi2, lrt_df, lrt_p, pToStars(lrt_p));

        % === Random effects variance components ===
        [~, ~, reStats] = covarianceParameters(mdl);

        % === Store model-level info for Supp Table 2 ===
        info.domain     = cogLabels{ci};
        info.outcome    = outVar;
        info.formula    = fmla;
        info.nObs       = mdl.NumObservations;   % what fitlme actually used
        info.nObs_filter = nObs;                 % what okAdj would have allowed
        info.nPatients  = nPat;
        info.AIC        = mdl.ModelCriterion.AIC;
        info.BIC        = mdl.ModelCriterion.BIC;
        info.logLik     = mdl.LogLikelihood;
        info.lrt_chi2   = lrt_chi2;
        info.lrt_df     = lrt_df;
        info.lrt_p      = lrt_p;
        info.reStats    = reStats;
        info.mdl        = mdl;
        allModelInfo{ci} = info;

    catch ME
        fprintf('    ERROR (%s): %s\n', cogLabels{ci}, ME.message);
        allFE{ci} = table();
        allModelInfo{ci} = struct('domain', cogLabels{ci}, 'error', ME.message);
    end
end

% -------------------------------------------------------------------------
% 4b.  ADDITIONAL MODELS — sensitivity to covariate choice & dropout
% -------------------------------------------------------------------------
%   (i)   Minimal model    : drops UPDRS3 and depfactor (likely mediators
%                            of apathy -> cognition). Yields the "total
%                            effect" of apathy, uncontaminated by
%                            over-adjustment.
%
%   (ii)  BL-adjusted model: includes the patient's baseline cognitive
%                            score as a covariate. Fit on time > 0 only;
%                            time is re-centred to "years since first
%                            follow-up" per patient so the intercept is
%                            not an extrapolation. Addresses reverse
%                            causation (baseline cognition -> apathy).
%
%   (iii) Pattern-mixture  : adds each patient's last-observed follow-up
%                            year (z-scored) and its interaction with
%                            time. Tests whether trajectories differ by
%                            dropout timing; if apathy terms are stable,
%                            findings are robust to informative dropout.
%
%  Three side-by-side CSVs + a comparison table are written to outputDir.
%  The main forest plot (Section 6) still uses allFE from the full model.
% -------------------------------------------------------------------------
fprintf('[Fig2] Additional sensitivity models (minimal / BL-adj / pattern-mixture)...\n');

altFE = struct();
altFE.minimal = cell(1, nCols);
altFE.blAdj   = cell(1, nCols);
altFE.patmix  = cell(1, nCols);

attrRows = {};

for ci = 1:nCols
    outVar = cogVars{ci};

    % ---------------------------------------------------------------------
    % (i) MINIMAL MODEL  (no UPDRS3, no depfactor)
    % ---------------------------------------------------------------------
    try
        okMin = okBase & ~isnan(lmeT.(outVar)) & ~isnan(lmeT.educ) & ...
                ~isnan(lmeT.nmiss) & ~isnan(lmeT.duration) & ~isnan(lmeT.ledd) & ...
                ~isundefined(lmeT.sex);
        lmeMin = lmeT(okMin, :);

        fmla_min = sprintf(['%s ~ apathy*icb*time' BLT ' + age + sex + educ + nmiss + ' ...
                            'duration + ledd + (1 + time | PATNO)'], outVar);
        mdl_min  = fitlme(lmeMin, fmla_min, 'FitMethod', 'REML');
        fe_min   = coeffToTable(mdl_min);
        ci_m     = coefCI(mdl_min);
        fe_min.CI_lo  = ci_m(:, 1);
        fe_min.CI_hi  = ci_m(:, 2);
        fe_min.Domain = repmat(cogLabels(ci), height(fe_min), 1);
        altFE.minimal{ci} = fe_min;

        writetable(fe_min, fullfile(outputDir, ...
            sprintf('%s_LME_%s_Minimal.csv', tableNums{ci}, outVar)));
        fprintf('  %-24s minimal     : %5d obs, %4d pts\n', ...
            cogLabels{ci}, height(lmeMin), numel(unique(lmeMin.PATNO)));
    catch ME
        fprintf('  ERROR minimal (%s): %s\n', cogLabels{ci}, ME.message);
        altFE.minimal{ci} = table();
    end

    % ---------------------------------------------------------------------
    % (ii) BASELINE-COGNITION-ADJUSTED MODEL
    % ---------------------------------------------------------------------
    try
        okBL0 = okBase & ~isnan(lmeT.(outVar)) & ~isnan(lmeT.educ) & ...
                ~isnan(lmeT.nmiss) & ~isnan(lmeT.duration) & ~isnan(lmeT.ledd) & ...
                ~isnan(lmeT.updrs3_score) & ~isundefined(lmeT.sex) & ...
                ~isnan(lmeT.depfactor);
        lmeBL = lmeT(okBL0, :);

        % Per-patient baseline score (time == 0 rows within the eligible set)
        blRows = lmeBL.time == 0;
        blTab  = lmeBL(blRows, {'PATNO', outVar});
        blTab.Properties.VariableNames{2} = 'bl_outcome';
        % De-duplicate (one row per patient, first baseline)
        [~, uIdx] = unique(blTab.PATNO, 'stable');
        blTab = blTab(uIdx, :);
        % Z-score baseline predictor
        blTab.bl_outcome = (blTab.bl_outcome - mean(blTab.bl_outcome, 'omitnan')) ...
                           ./ std(blTab.bl_outcome, 0, 'omitnan');

        % Drop baseline rows from the fit, attach bl_outcome to follow-ups
        lmeBL = lmeBL(lmeBL.time > 0, :);
        lmeBL = innerjoin(lmeBL, blTab, 'Keys', 'PATNO');

        % Re-centre time to "years since first observed follow-up" per patient.
        % This keeps the per-year slope identical to the main model but makes
        % the intercept interpretable (first follow-up) rather than an
        % extrapolation to baseline.
        [~, ~, gIdx] = unique(lmeBL.PATNO);
        firstFU      = accumarray(gIdx, lmeBL.time, [], @min);
        lmeBL.time_fu = lmeBL.time - firstFU(gIdx);

        fmla_bl = sprintf(['%s ~ apathy*icb*time_fu + age + sex + educ + nmiss + ' ...
                           'duration + ledd + updrs3_score + depfactor + bl_outcome + ' ...
                           '(1 + time_fu | PATNO)'], outVar);
        mdl_bl  = fitlme(lmeBL, fmla_bl, 'FitMethod', 'REML');
        fe_bl   = coeffToTable(mdl_bl);
        % Rename time_fu -> time in coefficient names so the comparison
        % table lines up with the main / minimal / pattern-mixture models.
        fe_bl.Name = regexprep(fe_bl.Name, 'time_fu', 'time');
        ci_b = coefCI(mdl_bl);
        fe_bl.CI_lo  = ci_b(:, 1);
        fe_bl.CI_hi  = ci_b(:, 2);
        fe_bl.Domain = repmat(cogLabels(ci), height(fe_bl), 1);
        altFE.blAdj{ci} = fe_bl;

        writetable(fe_bl, fullfile(outputDir, ...
            sprintf('%s_LME_%s_BLAdjusted.csv', tableNums{ci}, outVar)));
        fprintf('  %-24s BL-adjusted : %5d obs, %4d pts (time>0)\n', ...
            cogLabels{ci}, height(lmeBL), numel(unique(lmeBL.PATNO)));
    catch ME
        fprintf('  ERROR BL-adj  (%s): %s\n', cogLabels{ci}, ME.message);
        altFE.blAdj{ci} = table();
    end

    % ---------------------------------------------------------------------
    % (iii) PATTERN-MIXTURE MODEL (informative dropout sensitivity)
    % ---------------------------------------------------------------------
    try
        okPM = okBase & ~isnan(lmeT.(outVar)) & ~isnan(lmeT.educ) & ...
               ~isnan(lmeT.nmiss) & ~isnan(lmeT.duration) & ~isnan(lmeT.ledd) & ...
               ~isnan(lmeT.updrs3_score) & ~isundefined(lmeT.sex) & ...
               ~isnan(lmeT.depfactor);
        lmePM = lmeT(okPM, :);

        [gid, ~, gIdx] = unique(lmePM.PATNO);
        lastFU  = accumarray(gIdx, lmePM.time, [], @max);
        lmePM.lastFU   = lastFU(gIdx);
        lmePM.lastFU_z = (lmePM.lastFU - mean(lmePM.lastFU, 'omitnan')) ...
                         ./ std(lmePM.lastFU, 0, 'omitnan');

        fmla_pm = sprintf(['%s ~ apathy*icb*time' BLT ' + age + sex + educ + nmiss + ' ...
                           'duration + ledd + updrs3_score + depfactor + ' ...
                           'lastFU_z + lastFU_z:time + (1 + time | PATNO)'], outVar);
        mdl_pm  = fitlme(lmePM, fmla_pm, 'FitMethod', 'REML');
        fe_pm   = coeffToTable(mdl_pm);
        ci_p    = coefCI(mdl_pm);
        fe_pm.CI_lo  = ci_p(:, 1);
        fe_pm.CI_hi  = ci_p(:, 2);
        fe_pm.Domain = repmat(cogLabels(ci), height(fe_pm), 1);
        altFE.patmix{ci} = fe_pm;

        writetable(fe_pm, fullfile(outputDir, ...
            sprintf('%s_LME_%s_PatternMixture.csv', tableNums{ci}, outVar)));

        % Attrition summary for this outcome's sample
        attr.Domain        = cogLabels{ci};
        attr.N_patients    = numel(gid);
        attr.Median_lastFU = median(lastFU);
        attr.IQR_lastFU    = prctile(lastFU, 75) - prctile(lastFU, 25);
        attr.Pct_ge_5y     = round(100 * mean(lastFU >= 5), 1);
        attr.Pct_ge_10y    = round(100 * mean(lastFU >= 10), 1);
        attrRows{end+1}    = attr; %#ok

        fprintf('  %-24s pattern-mix : %4d pts, median lastFU=%.1f y (%.0f%% >=5y, %.0f%% >=10y)\n', ...
            cogLabels{ci}, attr.N_patients, attr.Median_lastFU, ...
            attr.Pct_ge_5y, attr.Pct_ge_10y);
    catch ME
        fprintf('  ERROR pat-mix (%s): %s\n', cogLabels{ci}, ME.message);
        altFE.patmix{ci} = table();
    end
end

% --- Attrition by group4 (Tpd_lme carries the group4 labels; time column is YEAR) ---
attrD = Tpd_lme(~isnan(Tpd_lme.group4), {'PATNO', 'group4', 'YEAR'});
[~, ~, gIdx]  = unique(attrD.PATNO);
lastPat       = accumarray(gIdx, attrD.YEAR, [], @max);
grpPat        = accumarray(gIdx, attrD.group4, [], @(x) x(1));
fprintf('\n  Follow-up by group4:\n');
for g = 1:4
    mask = grpPat == g;
    if sum(mask) == 0, continue; end
    fprintf('    %-10s n=%4d  median lastFU=%4.1f y   >=5y: %4.0f%%   >=10y: %4.0f%%\n', ...
        groupLabelsShort{g}, sum(mask), median(lastPat(mask)), ...
        100 * mean(lastPat(mask) >= 5), 100 * mean(lastPat(mask) >= 10));
end

if ~isempty(attrRows)
    attrTab = struct2table([attrRows{:}]);
    writetable(attrTab, fullfile(outputDir, 'Diagnostics_Attrition_Summary.csv'));
    fprintf('  -> Diagnostics_Attrition_Summary.csv\n');
end

% --- Side-by-side β/p comparison across all four model specifications ---
keyTerms = {'apathy', 'icb', 'time', 'apathy:icb', ...
            'apathy:time', 'icb:time', 'apathy:icb:time'};
cmpRows = {};
for ci = 1:nCols
    for ti = 1:numel(keyTerms)
        row = struct();
        row.Domain = cogLabels{ci};
        row.Term   = keyTerms{ti};
        row.Main_b = altGet(allFE{ci},         keyTerms{ti}, 'Estimate');
        row.Main_p = altGet(allFE{ci},         keyTerms{ti}, 'pValue');
        row.Min_b  = altGet(altFE.minimal{ci}, keyTerms{ti}, 'Estimate');
        row.Min_p  = altGet(altFE.minimal{ci}, keyTerms{ti}, 'pValue');
        row.BL_b   = altGet(altFE.blAdj{ci},   keyTerms{ti}, 'Estimate');
        row.BL_p   = altGet(altFE.blAdj{ci},   keyTerms{ti}, 'pValue');
        row.PM_b   = altGet(altFE.patmix{ci},  keyTerms{ti}, 'Estimate');
        row.PM_p   = altGet(altFE.patmix{ci},  keyTerms{ti}, 'pValue');
        cmpRows{end+1} = row; %#ok
    end
end
if ~isempty(cmpRows)
    cmpTab = struct2table([cmpRows{:}]);
    writetable(cmpTab, fullfile(outputDir, 'Diagnostics_AltModels_Comparison.csv'));
    fprintf('  -> Diagnostics_AltModels_Comparison.csv\n\n');
end

% -------------------------------------------------------------------------
% 5.  FDR CORRECTION ACROSS DOMAINS (same term, four tests)
% -------------------------------------------------------------------------
for ci = 1:nCols
    if ~isempty(allFE{ci})
        allFE{ci}.pFDR = allFE{ci}.pValue;
    end
end

for ti = 1:numel(termsOfInterest)
    rawP   = nan(1, nCols);
    rowIdx = nan(1, nCols);
    for ci = 1:nCols
        fe = allFE{ci};
        if isempty(fe), continue; end
        r = matchTerm(fe.Name, termsOfInterest{ti});
        if ~isempty(r)
            rawP(ci)   = fe.pValue(r);
            rowIdx(ci) = r;
        end
    end
    % Benjamini–Hochberg
    validMask = ~isnan(rawP);
    if sum(validMask) > 1
        pSub  = rawP(validMask);
        m     = numel(pSub);
        [~, sortI] = sort(pSub);
        ranks = 1:m;
        ranks(sortI) = ranks;
        pAdj  = min(1, pSub .* m ./ ranks);
        % Enforce monotonicity
        [~, descI] = sort(ranks, 'descend');
        pAdj_sorted = pAdj(descI);
        pAdj_sorted = cummin(pAdj_sorted);
        pAdj(descI) = pAdj_sorted;

        idxV = find(validMask);
        for vi = 1:numel(idxV)
            ci_idx = idxV(vi);
            allFE{ci_idx}.pFDR(rowIdx(ci_idx)) = pAdj(vi);
        end
    end
end

% -------------------------------------------------------------------------
% 6.  BOTTOM ROW — FOREST PLOTS  (panels E–H)
% -------------------------------------------------------------------------
for ci = 1:nCols
    ax = nexttile(nCols + ci);
    set(ax, 'Color', COL.bg); hold on;

    fe = allFE{ci};
    if isempty(fe), continue; end

    nTerms = numel(termsOfInterest);
    plotted = false(nTerms, 1);

    for ti = 1:nTerms
        r = matchTerm(fe.Name, termsOfInterest{ti});
        if isempty(r), continue; end

        plotted(ti) = true;

        est   = fe.Estimate(r);
        ci_lo = fe.CI_lo(r);
        ci_hi = fe.CI_hi(r);
        pCorr = fe.pFDR(r);
        yp    = nTerms - ti + 1;

        % Colour by interaction order
        nColons = sum(termsOfInterest{ti} == ':');
        switch nColons
            case 0, baseCol = fpCol.main;
            case 1, baseCol = fpCol.two;
            case 2, baseCol = fpCol.three;
        end

        isSig = (pCorr < 0.05);

        % CI whisker
        whiskerCol = baseCol;
        if ~isSig, whiskerCol = fpCol.ns; end
        plot(ax, [ci_lo ci_hi], [yp yp], '-', ...
            'Color', whiskerCol, 'LineWidth', 1.6);

        % Point estimate (filled if significant, open otherwise)
        if isSig
            plot(ax, est, yp, 'o', 'Color', baseCol, ...
                'MarkerFaceColor', baseCol, 'MarkerSize', 5, 'LineWidth', 0.6);
        else
            plot(ax, est, yp, 'o', 'Color', fpCol.ns, ...
                'MarkerFaceColor', COL.bg, 'MarkerSize', 5, 'LineWidth', 1);
        end
    end

    % Null-effect reference line
    xline(ax, 0, '-', 'Color', [0.78 0.78 0.78], 'LineWidth', 0.5);

    % --- Significance stars placed AFTER x-limits are finalised ---
    xRng = xlim(ax);
    starOffset = 0.03 * diff(xRng);

    for ti = 1:nTerms
        if ~plotted(ti), continue; end
        r     = matchTerm(fe.Name, termsOfInterest{ti});
        pCorr = fe.pFDR(r);
        if pCorr < 0.05
            yp = nTerms - ti + 1;
            nColons = sum(termsOfInterest{ti} == ':');
            switch nColons
                case 0, sCol = fpCol.main;
                case 1, sCol = fpCol.two;
                case 2, sCol = fpCol.three;
            end
            text(fe.CI_hi(r) + starOffset, yp, pToStars(pCorr), ...
                'FontSize', 7, 'Color', sCol, 'FontWeight', 'bold', ...
                'VerticalAlignment', 'middle');
        end
    end

    % Warn if any terms were missing
    if any(~plotted)
        missing = termsOfInterest(~plotted);
        fprintf('    Panel %s: terms not found: %s\n', ...
            char('E'+ci-1), strjoin(missing, ', '));
    end

    % Y-axis labels
    set(ax, 'YTick', 1:nTerms, 'YTickLabel', flip(termDisplayNames), ...
        'YLim', [0.3, nTerms + 0.7], ...
        'FontSize', 6.5, 'TickDir', 'out', 'Box', 'off', ...
        'XGrid', 'off', 'YGrid', 'off');

    xlabel('\beta  (95% CI)', 'FontSize', 7);
    if ci == 1
        ylabel('Fixed effect', 'FontSize', 7);
    end

    panelLetter = char('E' + ci - 1);
    title(ax, sprintf('\\bf%s\\rm  %s', panelLetter, cogLabels{ci}), ...
        'FontSize', 8, 'Color', COL.txt);

    formatAxis(ax, COL);
    hold off;
end

% -------------------------------------------------------------------------
% 7.  EXPORT FIGURE
% -------------------------------------------------------------------------
exportFigure(figSupp2, fullfile(outputDir, 'Fig2_Supp_4Domain_Trajectories_ForestPlots'));
fprintf('  -> Fig2_Supp_4Domain_Trajectories_ForestPlots saved.\n\n');


%% ========================================================================
%  PRIMARY FIGURE 2  (new layout)
%  A. bifactor_g trajectory by group4 (or vizMode split)
%  B. Forest plot with paired markers for bifactor_g and MoCA;
%     bracket compares Apathy x Time vs ICB x Time on bifactor_g (Wald)
%  C. Heatmap of three sub-domain LMEs (exec / memory / visuo); cells
%     coloured by beta; FDR-significant cells (BH within term, k = 3)
%     annotated with stars.
%
%  NOTE: bifactor_g and MoCA are BOTH treated as primary global outcomes
%  in panel B and use raw p-values (not the 4-domain FDR computed in
%  Section 5). The 3-domain FDR for the heatmap is computed below.
% =========================================================================
fprintf('[Fig2 Primary] Trajectory + paired forest + heatmap...\n');

% -------------------------------------------------------------------------
% Step 1: fit MoCA LMM (same covariate set as the main bifactor_g model)
% -------------------------------------------------------------------------
feMoca   = table();
mdl_moca = [];
try
    okM = okBase & ~isnan(lmeT.moca) & ...
          ~isnan(lmeT.apathy) & ~isnan(lmeT.icb) & ~isnan(lmeT.time) & ...
          ~isnan(lmeT.age) & ~isnan(lmeT.educ) & ~isnan(lmeT.nmiss) & ...
          ~isnan(lmeT.duration) & ~isnan(lmeT.ledd) & ~isnan(lmeT.updrs3_score) & ...
          ~isundefined(lmeT.sex) & ~isnan(lmeT.depfactor);
    lmeM   = lmeT(okM, :);
    % Match the bifactor-g random-effects specification exactly so that
    % the manuscript primary MoCA table and the supplementary sensitivity
    % table report the same β's. The (1|time) term adds a year-level random
    % intercept that absorbs cohort-wave variance; without it the linear
    % time slope is biased (MoCA shrinks from -0.24 to -0.15 in our data).
    fmla_m = ['moca ~ apathy*icb*time' BLT ' + age + sex + educ + nmiss + duration + ' ...
              'ledd + updrs3_score + depfactor + (1 + time | PATNO) + (1| time)'];
    mdl_moca = fitlme(lmeM, fmla_m, 'FitMethod', 'REML');
    feMoca   = coeffToTable(mdl_moca);
    ciM      = coefCI(mdl_moca);
    feMoca.CI_lo = ciM(:, 1);
    feMoca.CI_hi = ciM(:, 2);
    fprintf('  MoCA LME fitted: %d obs (filter), %d used by LME, %d patients\n', ...
        height(lmeM), mdl_moca.NumObservations, numel(unique(lmeM.PATNO)));
    writetable(feMoca, fullfile(outputDir, 'Table12_LME_moca_Adjusted.csv'));
catch ME
    fprintf('  ERROR MoCA LME: %s\n', ME.message);
end

% -------------------------------------------------------------------------
% Step 2: locate the global cognition outcome (used by panels A and C)
% -------------------------------------------------------------------------
gIdx = find(strcmp(cogVars, globalVar), 1);
if isempty(gIdx), gIdx = 1; end

% -------------------------------------------------------------------------
% Step 3: 3-domain FDR (sub-domains only) for the heatmap
% -------------------------------------------------------------------------
domIdx = setdiff(1:nCols, gIdx);     % indices of 3 sub-domains
nD     = numel(domIdx);
nT     = numel(termsOfInterest);

betaMat3 = nan(nD, nT);
pMat3    = nan(nD, nT);
for di = 1:nD
    fe = allFE{domIdx(di)};
    if isempty(fe), continue; end
    for ti = 1:nT
        r = matchTerm(fe.Name, termsOfInterest{ti});
        if isempty(r), continue; end
        betaMat3(di, ti) = fe.Estimate(r);
        pMat3(di, ti)    = fe.pValue(r);
    end
end

% BH adjustment within each term (k = 3). Keep everything as column
% vectors so the broadcast doesn't accidentally produce a matrix.
pAdjMat3 = pMat3;
for ti = 1:nT
    rawP  = pMat3(:, ti);              % column
    valid = ~isnan(rawP);
    if sum(valid) > 1
        pSub = rawP(valid);            % column
        m    = numel(pSub);
        [~, sortI]   = sort(pSub);     % column
        ranks        = (1:m)';         % column (was row -> broadcast bug)
        ranks(sortI) = ranks;
        pAdj         = min(1, pSub .* m ./ ranks);
        % monotonicity
        [~, descI]   = sort(ranks, 'descend');
        sortedAdj    = pAdj(descI);
        sortedAdj    = cummin(sortedAdj);
        pAdj(descI)  = sortedAdj;
        pAdjMat3(valid, ti) = pAdj;
    end
end

% -------------------------------------------------------------------------
% Step 3b: fit OBSERVED sub-domain LMEs (parallel set, no residualisation)
% and compute their 3-domain FDR for the second heatmap in panel C.
% -------------------------------------------------------------------------
obsVars   = {'obs_att_exec', 'obs_memory', 'obs_visuo'};
obsLabels = {'Executive & Attention','Memory','Visuospatial'};
nObs      = numel(obsVars);

allFEObs    = cell(1, nObs);
betaMat3Obs = nan(nObs, nT);
pMat3Obs    = nan(nObs, nT);

for di = 1:nObs
    outVar = obsVars{di};
    if ~ismember(outVar, lmeT.Properties.VariableNames)
        fprintf('  WARNING: %s not in lmeT, skipping observed-score model.\n', outVar);
        continue;
    end
    try
        okO = okBase & ~isnan(lmeT.(outVar)) & ...
              ~isnan(lmeT.educ) & ~isnan(lmeT.nmiss) & ...
              ~isnan(lmeT.duration) & ~isnan(lmeT.ledd) & ~isnan(lmeT.updrs3_score) & ...
              ~isundefined(lmeT.sex) & ~isnan(lmeT.depfactor);
        if controlForGlobal
            okO = okO & ~isnan(lmeT.(globalVar));
        end
        lmeO   = lmeT(okO, :);
        % Optionally control for global cognition (globalVar) inside the LMM:
        %   - true  -> apathy / ICB / time effects are conditional on g
        %   - false -> total effects on each raw sub-domain (no g covariate)
        if controlForGlobal
            fmla_o = sprintf(['%s ~ apathy*icb*time' BLT ' + %s + age + sex + educ + nmiss + ' ...
                              'duration + ledd + updrs3_score + depfactor + (1 + time | PATNO)'], ...
                              outVar, globalVar);
        else
            fmla_o = sprintf(['%s ~ apathy*icb*time' BLT ' + age + sex + educ + nmiss + ' ...
                              'duration + ledd + updrs3_score + depfactor + (1 + time | PATNO)'], ...
                              outVar);
        end
        mdl_o  = fitlme(lmeO, fmla_o, 'FitMethod', 'REML');
        feO    = coeffToTable(mdl_o);
        ciO    = coefCI(mdl_o);
        feO.CI_lo  = ciO(:, 1);
        feO.CI_hi  = ciO(:, 2);
        feO.Domain = repmat(obsLabels(di), height(feO), 1);
        allFEObs{di} = feO;
        for ti = 1:nT
            r = matchTerm(feO.Name, termsOfInterest{ti});
            if isempty(r), continue; end
            betaMat3Obs(di, ti) = feO.Estimate(r);
            pMat3Obs(di, ti)    = feO.pValue(r);
        end
        suffix = ternary(controlForGlobal, ...
            sprintf('Adjusted_ctrl_%s', globalVar), 'Adjusted_noCtrl');
        writetable(feO, fullfile(outputDir, ...
            sprintf('Table_LME_%s_%s.csv', outVar, suffix)));
        fprintf('  %-22s LME fitted: %d obs, %d patients\n', ...
            outVar, height(lmeO), numel(unique(lmeO.PATNO)));
    catch ME
        fprintf('  ERROR LME %s: %s\n', outVar, ME.message);
        allFEObs{di} = table();
    end
end

% BH adjustment within each term for OBSERVED (k = 3)
pAdjMat3Obs = pMat3Obs;
for ti = 1:nT
    rawP  = pMat3Obs(:, ti);
    valid = ~isnan(rawP);
    if sum(valid) > 1
        pSub = rawP(valid);
        m    = numel(pSub);
        [~, sortI]   = sort(pSub);
        ranks        = (1:m)';
        ranks(sortI) = ranks;
        pAdj         = min(1, pSub .* m ./ ranks);
        [~, descI]   = sort(ranks, 'descend');
        sortedAdj    = pAdj(descI);
        sortedAdj    = cummin(sortedAdj);
        pAdj(descI)  = sortedAdj;
        pAdjMat3Obs(valid, ti) = pAdj;
    end
end

% -------------------------------------------------------------------------
% Step 4: build the figure
%   Outer 1x2 layout. Only panels A and B in the primary figure.
%     Panel A: 2x3 nested grid -> main trajectory spans 2 rows x 2 cols
%              (left); two smaller sub-panels stack on the right
%              (apathy split top, ICB split bottom).
%     Panel B: standalone forest plot (globalVar + MoCA, paired markers).
%   Panel C (sub-domain heatmap) moved to its own supplementary figure
%   built right after fig2 is exported below.
% -------------------------------------------------------------------------
fig2 = figure('Units','centimeters','Position',[1 1 26 12.5],'Color',COL.bg);
% Helvetica everywhere; bold panel titles, but NORMAL-weight axis tick
% numbers and text annotations (n=N / endpoint labels). Tick numbers are
% deliberately not bold.
set(fig2, ...
    'DefaultAxesFontName','Helvetica','DefaultTextFontName','Helvetica', ...
    'DefaultAxesFontWeight','normal','DefaultTextFontWeight','normal', ...
    'DefaultAxesTitleFontWeight','bold');
% 1x5 outer: panel A spans 3 tiles, panel B spans 2 — gives A more width
% than B while keeping a clean side-by-side layout (each tile = width/5).
outerTL = tiledlayout(fig2, 1, 5, 'TileSpacing', 'loose', 'Padding', 'compact');

% Diverging colormap (used by the supplementary heatmap below; defined
% here so the recipe lives in one place).
nCM = 256; cmap_div = zeros(nCM, 3);
for ci = 1:nCM
    t = (ci - 1) / (nCM - 1);
    if t < 0.5
        s = t / 0.5;
        cmap_div(ci, :) = [0.12 + 0.88*s, 0.30 + 0.70*s, 0.85 + 0.15*s];
    else
        s = (t - 0.5) / 0.5;
        cmap_div(ci, :) = [1.0, 1.0 - 0.70*s, 1.0 - 0.85*s];
    end
end

% =========================================================================
% PANEL A : main trajectory (left, 2x2 of inner) + sub-panels stacked on
%           the right column (apathy split top, ICB split bottom).
% =========================================================================
gA = tiledlayout(outerTL, 2, 3, 'TileSpacing', 'compact', 'Padding', 'compact');
gA.Layout.Tile = 1;
gA.Layout.TileSpan = [1 3];

% --- Main trajectory (left column, spans 2 rows x 2 cols of inner grid) --
axMain = nexttile(gA, [2 2]);
set(axMain, 'Color', COL.bg); hold(axMain, 'on');
ph = gobjects(nVizGrps, 1);
for g = 1:nVizGrps
    rowFilter = (vizGrp == vizLevels(g));
    ph(g) = plotCPTraj(axMain, TpdSub, rowFilter, 'YEAR', cogVars{gIdx}, ...
                       vizColors(g, :), [], sigmaVal, '-');
end
yline(axMain, 0, '--', 'Color', [0.65 0.65 0.65], 'LineWidth', 0.6);
% Extra right margin (max+1.0 rather than +0.3) so the year-15 n=N label in
% the footer band clears the axes/tile edge instead of being clipped.
xlim(axMain, [-0.5 max(years)+1.0]);
title(axMain, sprintf('\\bfa\\rm  %s', cogLabels{gIdx}), ...
    'FontSize', 11, 'Color', COL.txt);
xlabel(axMain, 'Year', 'FontSize', 12, 'FontWeight', 'bold');
ylabel(axMain, 'Factor score (z)', 'FontSize', 12, 'FontWeight', 'bold');
set(axMain, 'FontSize', 10, 'Box', 'off', 'TickDir', 'out', ...
    'XTick', [0 5 10 15], 'XGrid', 'off', 'YGrid', 'off');
% Reserve a footer band INSIDE the axes for the n=N labels by extending
% ylim downward. Drawing inside the axes guarantees the labels render
% (Clipping=off outside a tiledlayout tile is unreliable in some setups).
yl = ylim(axMain);
range0 = diff(yl);
ylim(axMain, [yl(1) - 0.18*range0, yl(2) + 0.06*range0]);

years5 = [0 5 10 15];
nVec   = zeros(size(years5));
for yi = 1:numel(years5)
    % Cohort-based count (matches Table 1): patients in the four
    % phenotype groups (group4 in 1..4) with a visit at this year.
    % No outcome-missingness filter, so n reflects who is *in* the
    % cohort at that visit, not who contributes to the curve estimate.
    sel = abs(TpdSub.YEAR - years5(yi)) < 0.5 ...
        & ~isnan(TpdSub.group4);
    nVec(yi) = numel(unique(TpdSub.PATNO(sel)));
end

ylNew  = ylim(axMain);
ySep   = yl(1) - 0.04 * range0;     % thin separator (just below original data range)
yLabel = yl(1) - 0.12 * range0;     % n=N centred in the footer band

% Subtle separator above the label row
plot(axMain, [-0.5, max(years)+1.0], [ySep, ySep], '-', ...
    'Color', [0.82 0.82 0.82], 'LineWidth', 0.4);

% n=N labels (italic mid-grey, sitting inside the axes footer)
for yi = 1:numel(years5)
    if nVec(yi) == 0, continue; end
    text(axMain, years5(yi), yLabel, sprintf('n=%d', nVec(yi)), ...
        'HorizontalAlignment','center', 'VerticalAlignment','middle', ...
        'FontSize', 8.5, 'FontAngle', 'italic', ...
        'Color', [0.40 0.40 0.40]);
end

validPh = ph(isgraphics(ph));
validLg = vizLabels(isgraphics(ph));
lg = legend(validPh, validLg, 'Location', 'southwest', 'FontSize', 9, 'Box', 'off');
lg.ItemTokenSize = [14 10];
% Anchor the legend horizontally at southwest but raise it vertically so
% its bottom sits at data y ~= -2 (clear of the n=N footer band).
drawnow;                             % ensure layout is finalised
ylNow   = ylim(axMain);
axPos   = axMain.Position;           % normalized figure coords
yTarget = -2;
yFrac   = (yTarget - ylNow(1)) / diff(ylNow);
lg.Position(2) = axPos(2) + yFrac * axPos(4);
formatAxis(axMain, COL, 10); hold(axMain, 'off');

% Helper: plot a 2-group trajectory on a small sub-axes
plotSubTraj = @(ax, gVec, cols, hiLbl, loLbl) deal();  % (placeholder, computed below)

% --- Sub-panel a1: apathy split (bottom-left) ----------------------------
axA1 = nexttile(gA);
set(axA1, 'Color', COL.bg); hold(axA1, 'on');
gApathy = nan(height(TpdSub), 1);
gApathy(TpdSub.group4 == 1 | TpdSub.group4 == 3) = 1;   % Low apathy
gApathy(TpdSub.group4 == 2 | TpdSub.group4 == 4) = 2;   % High apathy
colApathy = [COL.groups(1, :); COL.apathy];
hPhA1 = gobjects(1, 2);
for g = 1:2
    rowFilter = (gApathy == g);
    hPhA1(g) = plotCPTraj(axA1, TpdSub, rowFilter, 'YEAR', cogVars{gIdx}, ...
                          colApathy(g, :), [], sigmaVal, '-');
end
yline(axA1, 0, '--', 'Color', [0.65 0.65 0.65], 'LineWidth', 0.5);
title(axA1, 'Apathy split', 'FontSize', 10, 'FontWeight', 'bold', 'Color', COL.txt);
xlabel(axA1, 'Year', 'FontSize', 9, 'FontWeight', 'bold');
set(axA1, 'FontSize', 9, 'Box', 'off', 'TickDir', 'out', 'XTick', [0 5 10 15]);
formatAxis(axA1, COL, 9); hold(axA1, 'off');

% --- Sub-panel a2: ICB split (bottom-right) ------------------------------
axA2 = nexttile(gA);
set(axA2, 'Color', COL.bg); hold(axA2, 'on');
gICB = nan(height(TpdSub), 1);
gICB(TpdSub.group4 == 1 | TpdSub.group4 == 2) = 1;      % Low ICB
gICB(TpdSub.group4 == 3 | TpdSub.group4 == 4) = 2;      % High ICB
colICB = [COL.groups(1, :); COL.icb];
hPhA2 = gobjects(1, 2);
for g = 1:2
    rowFilter = (gICB == g);
    hPhA2(g) = plotCPTraj(axA2, TpdSub, rowFilter, 'YEAR', cogVars{gIdx}, ...
                          colICB(g, :), [], sigmaVal, '-');
end
yline(axA2, 0, '--', 'Color', [0.65 0.65 0.65], 'LineWidth', 0.5);
title(axA2, 'ICB split', 'FontSize', 10, 'FontWeight', 'bold', 'Color', COL.txt);
xlabel(axA2, 'Year', 'FontSize', 9, 'FontWeight', 'bold');
set(axA2, 'FontSize', 9, 'Box', 'off', 'TickDir', 'out', 'XTick', [0 5 10 15]);
formatAxis(axA2, COL, 9); hold(axA2, 'off');

% -----------------------------------------------------------------------
% Sync sub-panel y-limits (consistent y-ticks across A1/A2) and place the
% legend labels at the y-position of each line's last plotted point.
% Pin the lower bound to -1.7; share the upper bound across both panels.
% -----------------------------------------------------------------------
yl1 = ylim(axA1); yl2 = ylim(axA2);
yMax = max(yl1(2), yl2(2));
ylSync = [-1.7, yMax];
ylim(axA1, ylSync); ylim(axA2, ylSync);

% Sparse, explicit y-ticks (instead of the dense auto every-0.2 set). Set
% AFTER ylim so they aren't reset. Values chosen to span the panel range.
subYTicks = [-1.5 -1 -0.5 0];
set(axA1, 'YTick', subYTicks); set(axA2, 'YTick', subYTicks);

% Extend xlim slightly so endpoint labels don't get clipped
xExt = max(years) + 1.6;
xlim(axA1, [-0.3 xExt]); xlim(axA2, [-0.3 xExt]);

% Endpoint-anchored labels: use each line's actual final y-coordinate
labA1   = {'Low','High'};
labA2   = {'Low','High'};
colA1   = [COL.groups(1,:)*0.6; COL.apathy];
colA2   = [COL.groups(1,:)*0.6; COL.icb];
lblOff  = 0.25;            % small horizontal offset past the last data point

for g = 1:2
    if isgraphics(hPhA1(g))
        xd = get(hPhA1(g), 'XData'); yd = get(hPhA1(g), 'YData');
        if ~isempty(xd)
            text(axA1, xd(end) + lblOff, yd(end), labA1{g}, ...
                'HorizontalAlignment','left', 'VerticalAlignment','middle', ...
                'FontSize', 8.5, 'FontWeight', 'bold', 'Color', colA1(g,:));
        end
    end
    if isgraphics(hPhA2(g))
        xd = get(hPhA2(g), 'XData'); yd = get(hPhA2(g), 'YData');
        if ~isempty(xd)
            text(axA2, xd(end) + lblOff, yd(end), labA2{g}, ...
                'HorizontalAlignment','left', 'VerticalAlignment','middle', ...
                'FontSize', 8.5, 'FontWeight', 'bold', 'Color', colA2(g,:));
        end
    end
end

% =========================================================================
% PANEL B : paired forest plot (bifactor_g + MoCA)
% =========================================================================
gB = tiledlayout(outerTL, 1, 1, 'TileSpacing', 'compact', 'Padding', 'compact');
gB.Layout.Tile = 4;
gB.Layout.TileSpan = [1 2];
ax = nexttile(gB);
set(ax, 'Color', COL.bg); hold(ax, 'on');

fe_bg = allFE{gIdx};                      % bifactor_g fixed-effects table
fe_mc = feMoca;                           % MoCA fixed-effects table
bgCol = [0.25 0.30 0.38];                 % dark slate  – bifactor_g (primary)
mcCol = [0.50 0.20 0.55];                 % deep purple – MoCA
yOff  = 0.18;

for ti = 1:nT
    yp = nT - ti + 1;

    if ~isempty(fe_bg)
        rB = matchTerm(fe_bg.Name, termsOfInterest{ti});
        if ~isempty(rB)
            eB  = fe_bg.Estimate(rB); loB = fe_bg.CI_lo(rB); hiB = fe_bg.CI_hi(rB);
            sigB = fe_bg.pValue(rB) < 0.05;
            plot(ax, [loB hiB], [yp+yOff yp+yOff], '-', 'Color', bgCol, 'LineWidth', 1.4);
            if sigB
                plot(ax, eB, yp+yOff, 'o', 'Color', bgCol, ...
                    'MarkerFaceColor', bgCol, 'MarkerSize', 5);
            else
                plot(ax, eB, yp+yOff, 'o', 'Color', bgCol, ...
                    'MarkerFaceColor', COL.bg, 'MarkerSize', 5, 'LineWidth', 1);
            end
        end
    end

    if ~isempty(fe_mc)
        rM = matchTerm(fe_mc.Name, termsOfInterest{ti});
        if ~isempty(rM)
            eM  = fe_mc.Estimate(rM); loM = fe_mc.CI_lo(rM); hiM = fe_mc.CI_hi(rM);
            sigM = fe_mc.pValue(rM) < 0.05;
            plot(ax, [loM hiM], [yp-yOff yp-yOff], '-', 'Color', mcCol, 'LineWidth', 1.4);
            if sigM
                plot(ax, eM, yp-yOff, 'd', 'Color', mcCol, ...
                    'MarkerFaceColor', mcCol, 'MarkerSize', 5);
            else
                plot(ax, eM, yp-yOff, 'd', 'Color', mcCol, ...
                    'MarkerFaceColor', COL.bg, 'MarkerSize', 5, 'LineWidth', 1);
            end
        end
    end
end

xline(ax, 0, '-', 'Color', [0.78 0.78 0.78], 'LineWidth', 0.5);
set(ax, 'YTick', 1:nT, 'YTickLabel', flip(termDisplayNames), ...
    'YLim', [0.3, nT+0.7], 'FontSize', 9.5, 'TickDir', 'out', 'Box', 'off');
xlabel(ax, '\beta_z (95% CI)', 'FontSize', 12, 'FontWeight', 'bold');
title(ax, sprintf('\\bfb\\rm  Global cognition: %s vs MoCA', globalVarLabel), ...
    'FontSize', 11, 'Color', COL.txt);

% Significance stars (raw p; both outcomes are primary)
xRng = xlim(ax);
starOffset = 0.02 * diff(xRng);
for ti = 1:nT
    yp = nT - ti + 1;
    if ~isempty(fe_bg)
        rB = matchTerm(fe_bg.Name, termsOfInterest{ti});
        if ~isempty(rB)
            stB = pToStars(fe_bg.pValue(rB));
            if ~isempty(stB)
                text(ax, fe_bg.CI_hi(rB) + starOffset, yp + yOff, stB, ...
                    'FontSize', 10, 'Color', bgCol, 'FontWeight', 'bold', ...
                    'VerticalAlignment', 'middle');
            end
        end
    end
    if ~isempty(fe_mc)
        rM = matchTerm(fe_mc.Name, termsOfInterest{ti});
        if ~isempty(rM)
            stM = pToStars(fe_mc.pValue(rM));
            if ~isempty(stM)
                text(ax, fe_mc.CI_hi(rM) + starOffset, yp - yOff, stM, ...
                    'FontSize', 10, 'Color', mcCol, 'FontWeight', 'bold', ...
                    'VerticalAlignment', 'middle');
            end
        end
    end
end
xlim(ax, [xRng(1), xRng(2) + 0.06 * diff(xRng)]);

hBg = plot(nan, nan, 'o', 'Color', bgCol, 'MarkerFaceColor', bgCol, 'MarkerSize', 5);
hMc = plot(nan, nan, 'd', 'Color', mcCol, 'MarkerFaceColor', mcCol, 'MarkerSize', 5);
lgB = legend([hBg hMc], {globalVarLabel, 'MoCA'}, ...
    'Location', 'southwest', 'FontSize', 9, 'Box', 'off');
lgB.ItemTokenSize = [12 9];
% Panel B height: pbaspect ratio [w h 1] constrains the plot box. Larger
% h relative to w => taller panel. Comment this line out to fill the tile.
% [6 5 1] => 9.6 x 8 cm   (the previous "too short" version)
% [6 6 1] => 9.6 x 9.6 cm (square, mild shrink)
% Currently no constraint -> panel B fills its tile.
% pbaspect(ax, [6 6 1]);
formatAxis(ax, COL, 9.5); hold(ax, 'off');

exportFigure(fig2, fullfile(outputDir, 'Fig2_Primary_Trajectory_Forest'));
fprintf('  -> Fig2_Primary_Trajectory_Forest saved.\n\n');

%%
% =========================================================================
% SAMPLE SIZE TABLE: total + per-group counts at each visit year.
% Cohort-based counts (matches Table 1 logic): unique patients with a
% visit at each integer year, regardless of whether the global outcome
% is observed. Group columns follow groupLabelsShort
% (Neither / Apathy / ICB / Both); Ungroupable captures patients with
% group4 == NaN (missing baseline apathy_composite or QUIP), so on every
% row Total = Neither + Apathy + ICB + Both + Ungroupable.
% =========================================================================
yrAll = unique(TpdSub.YEAR(~isnan(TpdSub.YEAR)));
yrAll = yrAll(yrAll >= 0);
nszRows = cell(numel(yrAll), 1);
for yi = 1:numel(yrAll)
    yr   = yrAll(yi);
    selY = abs(TpdSub.YEAR - yr) < 0.5;
    row  = struct();
    row.Year  = yr;
    row.Total = numel(unique(TpdSub.PATNO(selY)));
    for g = 1:4
        fld       = matlab.lang.makeValidName(groupLabelsShort{g});
        selG      = selY & TpdSub.group4 == g;
        row.(fld) = numel(unique(TpdSub.PATNO(selG)));
    end
    selU            = selY & isnan(TpdSub.group4);
    row.Ungroupable = numel(unique(TpdSub.PATNO(selU)));
    nszRows{yi} = row;
end
nszTab = struct2table([nszRows{:}]);

% Sanity check: per-row balance (Total == Neither+Apathy+ICB+Both+Ungroupable)
grpCols  = cellfun(@matlab.lang.makeValidName, groupLabelsShort, 'UniformOutput', false);
groupSum = sum(nszTab{:, grpCols}, 2) + nszTab.Ungroupable;
mismatch = nszTab.Total - groupSum;
if any(mismatch ~= 0)
    fprintf('  WARN: %d rows where Total != groups + Ungroupable. Max diff = %d.\n', ...
        sum(mismatch ~= 0), max(abs(mismatch)));
end

writetable(nszTab, fullfile(outputDir, 'Diagnostics_SampleSize_ByYear_ByGroup.csv'));
fprintf('  Sample size by year x group (Total + group4 + Ungroupable):\n');
disp(nszTab);
fprintf('  -> Diagnostics_SampleSize_ByYear_ByGroup.csv\n\n');


%% ========================================================================
%  SUPPLEMENTARY FIGURE (panel-C heatmap): sub-domain LMEs
%  Standalone figure so it can be referenced/printed independently.
%  Same data and aesthetics as the previous panel C: 3 sub-domains × 7
%  terms, β_z heatmap with BH-FDR within-term overlays. Whether the
%  sub-domain models partial out global cognition follows controlForGlobal.
% =========================================================================
fprintf('[Fig2 Supp - Heatmap] Sub-domain LME heatmap...\n');

figSupp2c = figure('Units','centimeters','Position',[1 1 14 9],'Color',COL.bg);
set(figSupp2c, 'DefaultAxesFontName','Helvetica','DefaultTextFontName','Helvetica');

axH = axes('Parent', figSupp2c, 'Color', COL.bg);
bp  = betaMat3Obs; bp(isnan(bp)) = 0;
mab = max(abs(bp(:))); if mab == 0, mab = 1; end
clims = [-1 1] * mab * 1.05;
imagesc(axH, bp, clims);
colormap(axH, cmap_div);
cb = colorbar(axH); cb.Label.String = '\beta_z'; cb.FontSize = 6;
hold(axH, 'on');

% Cell value + significance overlay
for di = 1:nObs
    for ti = 1:nT
        if isnan(betaMat3Obs(di, ti)), continue; end
        pF = pAdjMat3Obs(di, ti);
        if      isnan(pF),     st = '';
        elseif  pF < 0.001,    st = '***';
        elseif  pF < 0.01,     st = '**';
        elseif  pF < 0.05,     st = '*';
        elseif  pF < 0.10,     st = char(8224);
        else,                  st = '';
        end
        nB = (bp(di, ti) - clims(1)) / (clims(2) - clims(1));
        tc = ternary(nB > 0.3 && nB < 0.7, [0.15 0.15 0.15], [1 1 1]);
        text(axH, ti, di, sprintf('%.3f', betaMat3Obs(di, ti)), ...
            'HorizontalAlignment','center','VerticalAlignment','middle', ...
            'FontSize', 6, 'Color', tc);
        if ~isempty(st)
            text(axH, ti, di + 0.3, st, ...
                'HorizontalAlignment','center','VerticalAlignment','middle', ...
                'FontSize', 6.5, 'Color', tc, 'FontWeight', 'bold');
        end
        if ~isnan(pF) && pF < 0.05
            rectangle(axH, 'Position', [ti-0.5, di-0.5, 1, 1], ...
                'EdgeColor', [0.1 0.1 0.1], 'LineWidth', 2);
        elseif ~isnan(pF) && pF < 0.10
            rectangle(axH, 'Position', [ti-0.5, di-0.5, 1, 1], ...
                'EdgeColor', [0.4 0.4 0.4], 'LineWidth', 1, 'LineStyle', '--');
        end
    end
end

set(axH, 'XTick', 1:nT, 'XTickLabel', termDisplayNames, ...
    'YTick', 1:nObs, 'YTickLabel', obsLabels, ...
    'FontSize', 6.5, 'Box', 'off', 'TickDir', 'out');
xtickangle(axH, 35);
title(axH, '\bfc\rm  Sub-domain LMEs (\beta_z, FDR within term, k=3)', ...
    'FontSize', 9, 'Color', COL.txt);
if controlForGlobal
    subC = sprintf('observed scores, controlling for %s', globalVarLabel);
else
    subC = 'observed scores, no global covariate';
end
subtitle(axH, subC, 'FontSize', 7, 'FontAngle', 'italic', 'Color', COL.txt);
xlim(axH, [0.5, nT+0.5]); ylim(axH, [0.5, nObs+0.5]);
hold(axH, 'off');

exportFigure(figSupp2c, fullfile(outputDir, 'Fig2_Supp_SubdomainHeatmap'));
fprintf('  -> Fig2_Supp_SubdomainHeatmap saved.\n\n');


%% ========================================================================
%  FIGURE S1: MODEL DIAGNOSTICS & VARIABLE DISTRIBUTIONS  (OPTIMISED)
%  Run AFTER the LME models in Fig2_CFA_Combined.m
%  Scope: bifactor_g (primary) + MoCA only — sub-domain models are skipped.
%  Requires: lmeT, allModelInfo, allFE, okBase, cogVars, cogLabels, COL,
%            mdl_moca, feMoca (from §7 of Figure 2).
% =========================================================================
fprintf('[Diagnostics] Variable distributions and LME residual checks...\n');
tStart = tic;

% Clean up any workspace variables that shadow built-in functions
clear fitted residuals dd diagData;

maxPts = 1500;  % subsample cap for scatter-type plots (fast & visually identical)
rng(42);        % reproducible subsampling

% -------------------------------------------------------------------------
% 0.  SHARED SETUP
% -------------------------------------------------------------------------
contVars   = {'apathy','icb','depfactor','educ','nmiss','duration','ledd','updrs3_score'};
contLabels = {'Apathy (z)','ICB (z)','Depression (z)','Education (z)', ...
              'Neuro. miss. (z)','Disease dur. (z)','LEDD (z)','UPDRS III (z)'};

% Diagnostics focus: bifactor_g (primary) + MoCA (secondary global outcome).
% Sub-domain outcomes (executive, memory, visuospatial) are intentionally
% excluded — Figure S1 is meant to validate the global-cognition story.
% MoCA is fit upstream (§7 of Figure 2) as mdl_moca / feMoca; we wrap it
% into the same {info, fe-table} shape used by allModelInfo / allFE so the
% loops below can treat the two outcomes uniformly.
dxVars   = {globalVar, 'moca'};
dxLabels = {'Global Cognition (g)', 'MoCA'};
dxColors = [COL.groups(1,:); 0.50 0.20 0.55];   % slate + purple (matches Fig 2b)

dxInfo = cell(1, numel(dxVars));
dxFE   = cell(1, numel(dxVars));

% (1) bifactor_g comes straight from the §4 pipeline outputs.
gIdxDx = find(strcmp(cogVars, globalVar), 1);
if ~isempty(gIdxDx)
    dxInfo{1} = allModelInfo{gIdxDx};
    dxFE{1}   = allFE{gIdxDx};
end

% (2) MoCA: rebuild a matching info struct from mdl_moca / feMoca.
if exist('mdl_moca','var') && ~isempty(mdl_moca)
    mocaInfo.domain  = 'MoCA';
    mocaInfo.outcome = 'moca';
    mocaInfo.mdl     = mdl_moca;
    dxInfo{2} = mocaInfo;
    if exist('feMoca','var') && ~isempty(feMoca)
        fe_m = feMoca;
        fe_m.Domain = repmat({'MoCA'}, height(fe_m), 1);
        dxFE{2} = fe_m;
    end
else
    fprintf('  WARNING: mdl_moca not found; MoCA diagnostics will be skipped.\n');
end

blIdx = okBase & (lmeT.time == 0);
lmeB  = lmeT(blIdx, :);
fprintf('  Baseline rows: %d\n', height(lmeB));

% Pre-extract all baseline vectors once (avoids repeated table indexing)
blData = struct();
for vi = 1:numel(contVars)
    blData.(contVars{vi}) = lmeB.(contVars{vi});
end
for ci = 1:numel(dxVars)
    blData.(dxVars{ci}) = lmeB.(dxVars{ci});
end

% Fast subsampling helper
subsamp = @(x, n) x(randperm(numel(x), min(n, numel(x))));

% -------------------------------------------------------------------------
% 1.  PREDICTOR & OUTCOME DISTRIBUTIONS
% -------------------------------------------------------------------------
figDist = figure('Units','centimeters','Position',[1 1 30 20],'Color',COL.bg);
set(figDist,'DefaultAxesFontName','Helvetica','DefaultTextFontName','Helvetica');
tl = tiledlayout(3, 4, 'TileSpacing','compact', 'Padding','compact');
title(tl, 'Predictor & Outcome Distributions (baseline)', ...
    'FontSize', 10, 'FontWeight', 'bold', 'Color', COL.accent);

% Predictors
for vi = 1:numel(contVars)
    ax = nexttile; hold on;
    v = blData.(contVars{vi});
    v = v(~isnan(v));

    histogram(v, 30, 'Normalization','pdf', ...
        'FaceColor', COL.groups(1,:), 'EdgeColor','none', 'FaceAlpha', 0.5);
    if numel(v) > 10
        [f, xi] = ksdensity(v);
        plot(xi, f, '-', 'Color', COL.accent, 'LineWidth', 1.5);
    end

    sk = skewness(v); ku = kurtosis(v);
    text(0.97, 0.95, sprintf('n=%d\nsk=%.2f\nku=%.2f', numel(v), sk, ku), ...
        'Units','normalized', 'FontSize', 6, 'HorizontalAlignment','right', ...
        'VerticalAlignment','top', 'Color', COL.txt);

    title(contLabels{vi}, 'FontSize', 7.5, 'Color', COL.txt);
    set(ax, 'FontSize', 6, 'Box', 'off', 'TickDir', 'out', ...
        'XGrid', 'off', 'YGrid', 'off');
    ylabel('Density', 'FontSize', 6);
    formatAxis(ax, COL); hold off;
end

% Outcomes (all time-points) — bifactor_g + MoCA only
for ci = 1:numel(dxVars)
    ax = nexttile; hold on;
    v = lmeT.(dxVars{ci});
    v = v(okBase & ~isnan(v));

    histogram(v, 40, 'Normalization','pdf', ...
        'FaceColor', dxColors(ci,:), 'EdgeColor','none', 'FaceAlpha', 0.5);
    if numel(v) > 10
        [f, xi] = ksdensity(v);
        plot(xi, f, '-', 'Color', COL.accent, 'LineWidth', 1.5);
    end

    sk = skewness(v); ku = kurtosis(v);
    text(0.97, 0.95, sprintf('n=%d\nsk=%.2f\nku=%.2f', numel(v), sk, ku), ...
        'Units','normalized', 'FontSize', 6, 'HorizontalAlignment','right', ...
        'VerticalAlignment','top', 'Color', COL.txt);

    title(dxLabels{ci}, 'FontSize', 8.5, 'FontWeight', 'bold', 'Color', COL.txt);
    set(ax, 'FontSize', 6, 'Box', 'off', 'TickDir', 'out', ...
        'XGrid', 'off', 'YGrid', 'off');
    ylabel('Density', 'FontSize', 6);
    formatAxis(ax, COL); hold off;
end

exportFigure(figDist, fullfile(outputDir, 'FigS1a_Variable_Distributions'));
fprintf('  -> FigS1a_Variable_Distributions (%.1fs)\n', toc(tStart));

% -------------------------------------------------------------------------
% 2.  LME RESIDUAL DIAGNOSTICS
% -------------------------------------------------------------------------
% Pre-extract ALL residuals and fitted values before opening the figure.
% This separates the slow model-query step from the fast plotting step.
fprintf('  Extracting residuals from %d models...\n', numel(dxVars));

diagData = struct();
for ci = 1:numel(dxVars)
    info = dxInfo{ci};
    if isempty(info) || ~isfield(info, 'mdl') || isempty(info.mdl), continue; end
    mdl = info.mdl;
    dd = struct();
    dd.condResid = residuals(mdl, 'ResidualType', 'Pearson');
    dd.fittedVals = predict(mdl);
    dd.timeVec   = mdl.Variables.time;
    diagData(ci).d = dd;
end
fprintf('  Residuals extracted (%.1fs)\n', toc(tStart));

figResid = figure('Units','centimeters','Position',[1 1 30 14],'Color',COL.bg);
set(figResid,'DefaultAxesFontName','Helvetica','DefaultTextFontName','Helvetica');
tl2 = tiledlayout(numel(dxVars), 4, 'TileSpacing','compact', 'Padding','compact');
title(tl2, 'LME Residual Diagnostics', ...
    'FontSize', 10, 'FontWeight', 'bold', 'Color', COL.accent);

for ci = 1:numel(dxVars)
    if ci > numel(diagData) || ~isfield(diagData(ci), 'd'), continue; end
    dd = diagData(ci).d;
    condResid  = dd.condResid;
    fittedVals = dd.fittedVals;
    timeVec    = dd.timeVec;
    nAll       = numel(condResid);
    domCol     = dxColors(ci, :);

    % Subsample indices (shared across panels for this model)
    if nAll > maxPts
        sIdx = sort(randperm(nAll, maxPts));
    else
        sIdx = 1:nAll;
    end

    % --- 2a. Residuals vs Fitted (subsampled, opaque dots) ---
    ax = nexttile((ci-1)*4 + 1); hold on;
    plot(fittedVals(sIdx), condResid(sIdx), '.', ...
        'Color', [domCol 0.25], 'MarkerSize', 3);
    yline(0, '-', 'Color', [0.6 0.6 0.6], 'LineWidth', 0.5);

    % Smoother on FULL data (cheap: movmedian on sorted vector)
    [fSort, sI] = sort(fittedVals);
    winSz = max(50, round(nAll/30));
    plot(fSort, movmedian(condResid(sI), winSz), '-', ...
        'Color', COL.accent, 'LineWidth', 1.2);

    title(sprintf('%s — Resid vs Fitted', dxLabels{ci}), 'FontSize', 7, 'Color', COL.txt);
    xlabel('Fitted', 'FontSize', 6); ylabel('Pearson residual', 'FontSize', 6);
    set(ax, 'FontSize', 5.5, 'Box', 'off', 'TickDir', 'out');
    formatAxis(ax, COL); hold off;

    % --- 2b. Q–Q plot (subsampled for display, line from full data) ---
    ax = nexttile((ci-1)*4 + 2); hold on;

    % Subsample for dots
    rSub   = sort(condResid(sIdx));
    nSub   = numel(rSub);
    thSub  = norminv(((1:nSub)' - 0.5) / nSub);
    plot(thSub, rSub, '.', 'Color', [domCol 0.25], 'MarkerSize', 3);

    % Reference line from full-data quartiles
    q1 = prctile(condResid, 25); q3 = prctile(condResid, 75);
    t1 = norminv(0.25); t3 = norminv(0.75);
    sl_qq = (q3 - q1) / (t3 - t1);
    in_qq = q1 - sl_qq * t1;
    xl = xlim;
    plot(xl, in_qq + sl_qq*xl, '-', 'Color', COL.accent, 'LineWidth', 1);

    title('Q–Q (Pearson)', 'FontSize', 7, 'Color', COL.txt);
    xlabel('Theoretical', 'FontSize', 6); ylabel('Sample', 'FontSize', 6);
    set(ax, 'FontSize', 5.5, 'Box', 'off', 'TickDir', 'out');
    formatAxis(ax, COL); hold off;

    % --- 2c. Residuals vs Time (IQR summary, no scatter) ---
    ax = nexttile((ci-1)*4 + 3); hold on;
    uYears = unique(timeVec);
    for yi = 1:numel(uYears)
        yr = uYears(yi);
        rY = condResid(timeVec == yr);
        if numel(rY) < 5, continue; end
        q25 = prctile(rY, 25); q75 = prctile(rY, 75);
        q05 = prctile(rY, 5);  q95 = prctile(rY, 95);
        med = median(rY);
        % 5th–95th whisker
        plot([yr yr], [q05 q95], '-', 'Color', [domCol 0.3], 'LineWidth', 0.8);
        % IQR box
        rectangle('Position', [yr-0.15, q25, 0.3, q75-q25], ...
            'FaceColor', [domCol 0.4], 'EdgeColor', domCol, 'LineWidth', 0.6);
        % Median
        plot(yr, med, 's', 'Color', COL.accent, 'MarkerFaceColor', COL.accent, ...
            'MarkerSize', 3.5);
    end
    yline(0, '-', 'Color', [0.6 0.6 0.6], 'LineWidth', 0.5);

    title('Resid vs Time', 'FontSize', 7, 'Color', COL.txt);
    xlabel('Year', 'FontSize', 6); ylabel('Pearson residual', 'FontSize', 6);
    set(ax, 'FontSize', 5.5, 'Box', 'off', 'TickDir', 'out');
    formatAxis(ax, COL); hold off;

    % --- 2d. Scale–location ---
    ax = nexttile((ci-1)*4 + 4); hold on;
    sqrtAbsR = sqrt(abs(condResid));
    plot(fittedVals(sIdx), sqrtAbsR(sIdx), '.', ...
        'Color', [domCol 0.25], 'MarkerSize', 3);
    plot(fSort, movmedian(sqrtAbsR(sI), winSz), '-', ...
        'Color', COL.accent, 'LineWidth', 1.2);

    title('Scale–Location', 'FontSize', 7, 'Color', COL.txt);
    xlabel('Fitted', 'FontSize', 6); ylabel('\surd|Pearson resid|', 'FontSize', 6);
    set(ax, 'FontSize', 5.5, 'Box', 'off', 'TickDir', 'out');
    formatAxis(ax, COL); hold off;
end

exportFigure(figResid, fullfile(outputDir, 'FigS1b_LME_Residual_Diagnostics'));
fprintf('  -> FigS1b_LME_Residual_Diagnostics (%.1fs)\n', toc(tStart));

% -------------------------------------------------------------------------
% 3.  RANDOM EFFECTS Q–Q (fast: ~1400 points per panel, no transparency)
% -------------------------------------------------------------------------
nDx = numel(dxVars);
figRE = figure('Units','centimeters','Position',[1 1 18 14],'Color',COL.bg);
set(figRE,'DefaultAxesFontName','Helvetica','DefaultTextFontName','Helvetica');
tl3 = tiledlayout(2, nDx, 'TileSpacing','compact', 'Padding','compact');
title(tl3, 'Random Effects: Intercept & Slope Distributions', ...
    'FontSize', 10, 'FontWeight', 'bold', 'Color', COL.accent);

for ci = 1:nDx
    info = dxInfo{ci};
    if isempty(info) || ~isfield(info, 'mdl') || isempty(info.mdl), continue; end
    mdl = info.mdl;

    reTab   = randomEffects(mdl);
    reInt   = reTab(1:2:end);
    reSlope = reTab(2:2:end);
    domCol  = dxColors(ci,:);

    % --- Intercept Q–Q ---
    ax = nexttile(ci); hold on;
    n = numel(reInt);
    sortedI = sort(reInt);
    theorQ  = norminv(((1:n)' - 0.5) / n);
    plot(theorQ, sortedI, '.', 'Color', domCol, 'MarkerSize', 4);

    q1 = prctile(reInt, 25); q3 = prctile(reInt, 75);
    t1 = norminv(0.25); t3 = norminv(0.75);
    sl = (q3 - q1)/(t3 - t1); intc = q1 - sl*t1;
    xl = xlim; plot(xl, intc + sl*xl, '-', 'Color', COL.accent, 'LineWidth', 1);

    title(sprintf('%s — Intercept', dxLabels{ci}), 'FontSize', 7, 'Color', COL.txt);
    xlabel('Theoretical', 'FontSize', 6); ylabel('Random int.', 'FontSize', 6);
    text(0.03, 0.95, sprintf('sk=%.2f\nku=%.2f', skewness(reInt), kurtosis(reInt)), ...
        'Units','normalized', 'FontSize', 5.5, 'VerticalAlignment','top', 'Color', COL.txt);
    set(ax, 'FontSize', 5.5, 'Box', 'off', 'TickDir', 'out');
    formatAxis(ax, COL); hold off;

    % --- Slope Q–Q ---
    ax = nexttile(nDx + ci); hold on;
    n = numel(reSlope);
    sortedS = sort(reSlope);
    theorQ  = norminv(((1:n)' - 0.5) / n);
    plot(theorQ, sortedS, '.', 'Color', domCol, 'MarkerSize', 4);

    q1 = prctile(reSlope, 25); q3 = prctile(reSlope, 75);
    sl = (q3 - q1)/(t3 - t1); intc = q1 - sl*t1;
    xl = xlim; plot(xl, intc + sl*xl, '-', 'Color', COL.accent, 'LineWidth', 1);

    title(sprintf('%s — Slope', dxLabels{ci}), 'FontSize', 7, 'Color', COL.txt);
    xlabel('Theoretical', 'FontSize', 6); ylabel('Random slope', 'FontSize', 6);
    text(0.03, 0.95, sprintf('sk=%.2f\nku=%.2f', skewness(reSlope), kurtosis(reSlope)), ...
        'Units','normalized', 'FontSize', 5.5, 'VerticalAlignment','top', 'Color', COL.txt);
    set(ax, 'FontSize', 5.5, 'Box', 'off', 'TickDir', 'out');
    formatAxis(ax, COL); hold off;
end

exportFigure(figRE, fullfile(outputDir, 'FigS1c_RandomEffects_QQ'));
fprintf('  -> FigS1c_RandomEffects_QQ (%.1fs)\n', toc(tStart));

% -------------------------------------------------------------------------
% 4.  COLLINEARITY (VIF) — instantaneous, no plotting
% -------------------------------------------------------------------------
fprintf('\n  Variance Inflation Factors (baseline, complete cases):\n');

vifVars   = {'apathy','icb','depfactor','educ','nmiss','duration','ledd','updrs3_score'};
blComplete = blIdx;
for vi = 1:numel(vifVars)
    blComplete = blComplete & ~isnan(lmeT.(vifVars{vi}));
end
blComplete = blComplete & ~isundefined(lmeT.sex);

Xmat = zeros(sum(blComplete), numel(vifVars));
for vi = 1:numel(vifVars)
    Xmat(:, vi) = lmeT.(vifVars{vi})(blComplete);
end

R = corrcoef(Xmat);
VIF = diag(inv(R));

vifTab = table(vifVars(:), round(VIF(:),2), 'VariableNames', {'Predictor','VIF'});
for vi = 1:height(vifTab)
    flag = '';
    if vifTab.VIF(vi) > 10, flag = ' *** SEVERE';
    elseif vifTab.VIF(vi) > 5, flag = ' ** MODERATE';
    elseif vifTab.VIF(vi) > 2.5, flag = ' * MILD';
    end
    fprintf('    %-20s  VIF = %5.2f%s\n', vifTab.Predictor{vi}, vifTab.VIF(vi), flag);
end
writetable(vifTab, fullfile(outputDir, 'Diagnostics_VIF.csv'));

% -------------------------------------------------------------------------
% 5.  BIVARIATE RELATIONSHIPS — HEATMAP + KEY SCATTER PAIRS ONLY
% -------------------------------------------------------------------------
% The old 4×8 scatter matrix had 32 panels of ~1400 transparent dots each
% — extremely slow.  Replace with:
%   (a) A Spearman correlation heatmap (all predictors × all outcomes)
%   (b) Scatter plots for apathy & ICB vs each outcome only (8 panels)

allVars   = [contVars, dxVars];
allLabels = [contLabels, dxLabels];
nV = numel(allVars);

% Build correlation matrix (baseline, pairwise complete)
corrMat = nan(nV);
pMat    = nan(nV);
for i = 1:nV
    vi = blData.(allVars{i});
    for j = i:nV
        vj = blData.(allVars{j});
        ok = ~isnan(vi) & ~isnan(vj);
        if sum(ok) > 10
            [corrMat(i,j), pMat(i,j)] = corr(vi(ok), vj(ok), 'Type','Spearman');
            corrMat(j,i) = corrMat(i,j);
            pMat(j,i) = pMat(i,j);
        end
    end
end

% --- 5a. Heatmap ---
figHeat = figure('Units','centimeters','Position',[1 1 22 18],'Color',COL.bg);
set(figHeat,'DefaultAxesFontName','Helvetica','DefaultTextFontName','Helvetica');
ax = axes('Color', COL.bg); hold on;

imagesc(corrMat, [-0.5 0.5]);
colormap(ax, [linspace(0.17,1,128)' linspace(0.50,1,128)' linspace(0.55,1,128)'; ...
              linspace(1,0.72,128)' linspace(1,0.38,128)' linspace(1,0.20,128)']);
cb = colorbar; cb.Label.String = 'Spearman \rho';
cb.Label.FontSize = 8;

% Annotate significant cells
for i = 1:nV
    for j = 1:nV
        if i == j, continue; end
        if ~isnan(corrMat(i,j))
            col = COL.txt;
            if abs(corrMat(i,j)) > 0.35, col = [1 1 1]; end
            str = sprintf('%.2f', corrMat(i,j));
            if pMat(i,j) < 0.001, str = [str '***'];
            elseif pMat(i,j) < 0.01, str = [str '**'];
            elseif pMat(i,j) < 0.05, str = [str '*'];
            end
            text(j, i, str, 'HorizontalAlignment','center', ...
                'FontSize', 5, 'Color', col);
        end
    end
end

set(ax, 'XTick', 1:nV, 'XTickLabel', allLabels, 'XTickLabelRotation', 45, ...
    'YTick', 1:nV, 'YTickLabel', allLabels, ...
    'FontSize', 6.5, 'TickDir', 'out', 'Box', 'off');
title('Spearman Correlation Matrix (baseline)', ...
    'FontSize', 9, 'FontWeight', 'bold', 'Color', COL.accent);
formatAxis(ax, COL); hold off;

exportFigure(figHeat, fullfile(outputDir, 'FigS1d_Correlation_Heatmap'));
fprintf('  -> FigS1d_Correlation_Heatmap (%.1fs)\n', toc(tStart));

% --- 5b. Key scatter pairs: apathy & ICB vs each outcome (8 panels) ---
keyPreds = {'apathy','icb'};
keyLabels = {'Apathy (z)','ICB (z)'};

figKey = figure('Units','centimeters','Position',[1 1 18 12],'Color',COL.bg);
set(figKey,'DefaultAxesFontName','Helvetica','DefaultTextFontName','Helvetica');
tl5 = tiledlayout(2, numel(dxVars), 'TileSpacing','compact', 'Padding','compact');
title(tl5, 'Key Predictor–Outcome Relationships (baseline)', ...
    'FontSize', 9, 'FontWeight', 'bold', 'Color', COL.accent);

for pi = 1:numel(keyPreds)
    xAll = blData.(keyPreds{pi});
    for ci = 1:numel(dxVars)
        ax = nexttile; hold on;
        yAll = blData.(dxVars{ci});
        ok   = ~isnan(xAll) & ~isnan(yAll);
        x = xAll(ok); y = yAll(ok);

        % Subsample for dots
        if numel(x) > maxPts
            si = randperm(numel(x), maxPts);
        else
            si = 1:numel(x);
        end
        plot(x(si), y(si), '.', 'Color', [dxColors(ci,:) 0.3], 'MarkerSize', 4);

        % Smoother on full data
        [xS, sI] = sort(x);
        winSz = max(30, round(numel(xS)/20));
        plot(xS, movmedian(y(sI), winSz), '-', 'Color', COL.accent, 'LineWidth', 1.3);

        % Annotation
        [rho, pval] = corr(x, y, 'Type', 'Spearman');
        text(0.03, 0.95, sprintf('r_s=%.2f, p=%.2g', rho, pval), ...
            'Units','normalized', 'FontSize', 5.5, ...
            'VerticalAlignment','top', 'Color', COL.txt);

        if pi == 1
            title(dxLabels{ci}, 'FontSize', 7.5, 'Color', COL.txt);
        end
        if ci == 1
            ylabel(keyLabels{pi}, 'FontSize', 7, 'FontWeight', 'bold', 'Color', COL.accent);
        end

        set(ax, 'FontSize', 5.5, 'Box', 'off', 'TickDir', 'out');
        formatAxis(ax, COL); hold off;
    end
end

exportFigure(figKey, fullfile(outputDir, 'FigS1e_Key_Scatter'));
fprintf('  -> FigS1e_Key_Scatter (%.1fs)\n', toc(tStart));

% -------------------------------------------------------------------------
% 6.  DESCRIPTIVE STATISTICS TABLE
% -------------------------------------------------------------------------
fprintf('\n  Compiling descriptive statistics table...\n');

descRows = cell(numel(allVars), 1);
for vi = 1:numel(allVars)
    v = blData.(allVars{vi});
    v = v(~isnan(v));

    row.Variable = allLabels{vi};
    row.N        = numel(v);
    row.Mean     = round(mean(v), 4);
    row.SD       = round(std(v), 4);
    row.Median   = round(median(v), 4);
    row.IQR_lo   = round(prctile(v, 25), 4);
    row.IQR_hi   = round(prctile(v, 75), 4);
    row.Min      = round(min(v), 4);
    row.Max      = round(max(v), 4);
    row.Skewness = round(skewness(v), 3);
    row.Kurtosis = round(kurtosis(v), 3);

    descRows{vi} = row;
end

descTab = struct2table([descRows{:}]);
writetable(descTab, fullfile(outputDir, 'Diagnostics_DescriptiveStats.csv'));
fprintf('  -> Diagnostics_DescriptiveStats.csv\n');

% -------------------------------------------------------------------------
% 7.  CONSOLE SUMMARY
% -------------------------------------------------------------------------
fprintf('\n  === TRANSFORMATION RECOMMENDATIONS ===\n');
for vi = 1:numel(allVars)
    v = blData.(allVars{vi});
    v = v(~isnan(v));
    sk = skewness(v); ku = kurtosis(v);

    flag = '  OK';
    if abs(sk) > 2
        flag = '  >> SEVERE SKEW — consider log or sqrt transform';
    elseif abs(sk) > 1
        flag = '  > MODERATE SKEW — inspect residuals; transform if residuals poor';
    end
    if ku > 7, flag = [flag, ' | HEAVY TAILS']; end %#ok

    fprintf('    %-30s  sk=%+6.2f  ku=%5.2f  %s\n', allLabels{vi}, sk, ku, flag);
end

fprintf('\n  NOTE: With n > 1000 and LME, moderate predictor skew is tolerable\n');
fprintf('  if RESIDUALS are well-behaved. Prioritise FigS1b over FigS1a.\n');

% -------------------------------------------------------------------------
% 8.  SENSITIVITY ANALYSIS: WINSORISED RE-FIT
% -------------------------------------------------------------------------
% Cap all continuous predictors at their 1st and 99th percentiles, re-fit
% each LME, and compare coefficients to the original models.  If the key
% terms (apathy, time:apathy, time:apathy:icb) shift materially, the
% original results are sensitive to extreme observations.

fprintf('\n  === SENSITIVITY: WINSORISED PREDICTORS (1st–99th pctile) ===\n');

winsorise = @(x) max(min(x, prctile(x, 99)), prctile(x, 1));

% Build winsorised copy of lmeT
lmeW = lmeT;
wVars = {'apathy','icb','depfactor','educ','nmiss','duration','ledd','updrs3_score'};
for vi = 1:numel(wVars)
    v = lmeW.(wVars{vi});
    okV = ~isnan(v);
    v(okV) = winsorise(v(okV));
    lmeW.(wVars{vi}) = v;
    nClipped = sum(lmeT.(wVars{vi})(okV) ~= v(okV));
    fprintf('    %-15s  %d values clipped\n', wVars{vi}, nClipped);
end

% Re-fit each domain model with winsorised predictors
normTerm = @(s) strjoin(sort(strsplit(strtrim(regexprep(lower(s), '\s*[:\*]\s*', ':')), ':')), ':');
matchTerm = @(feNames, target) find(cellfun( ...
    @(x) strcmp(normTerm(x), normTerm(target)), feNames), 1);

termsOfInterest = {'apathy','icb','time', ...
                   'apathy:icb','apathy:time','icb:time', ...
                   'apathy:icb:time'};
termDisplayNames = {'Apathy','ICB','Time', ...
                    'Apathy \times ICB','Apathy \times Time','ICB \times Time', ...
                    'Apathy \times ICB \times Time'};

sensRows = {};
for ci = 1:numel(dxVars)
    outVar = dxVars{ci};

    okAdj = okBase & ~isnan(lmeW.(outVar)) & ...
            ~isnan(lmeW.apathy) & ~isnan(lmeW.icb) & ~isnan(lmeW.time) & ...
            ~isnan(lmeW.age) & ~isnan(lmeW.educ) & ~isnan(lmeW.nmiss) & ...
            ~isnan(lmeW.duration) & ~isnan(lmeW.ledd) & ~isnan(lmeW.updrs3_score) & ...
            ~isundefined(lmeW.sex) & ~isnan(lmeW.depfactor);
    lmeD_w = lmeW(okAdj, :);

    fprintf('\n    %s (winsorised) — %d obs, %d patients\n', ...
        dxLabels{ci}, height(lmeD_w), numel(unique(lmeD_w.PATNO)));

    try
        % Use the primary model's full formula (incl. (1|time) random
        % intercept) so β_primary vs β_winsor isolates the winsorisation
        % effect, not the RE-spec change. char(mdl.Formula) returns the
        % expanded formula string that fitlme will re-accept.
        info = dxInfo{ci};
        if isempty(info) || ~isfield(info,'mdl') || isempty(info.mdl)
            fprintf('    SKIP (%s): no primary model to inherit formula from.\n', dxLabels{ci});
            continue;
        end
        fmla  = char(info.mdl.Formula);
        mdl_w = fitlme(lmeD_w, fmla, 'FitMethod', 'REML');
        fe_w  = coeffToTable(mdl_w);
        ci_w  = coefCI(mdl_w);
        fe_w.CI_lo = ci_w(:,1);
        fe_w.CI_hi = ci_w(:,2);

        % Compare with original model
        fe_orig = dxFE{ci};
        if isempty(fe_orig)
            fprintf('    SKIP (%s): no original fixed-effects table to compare against.\n', dxLabels{ci});
            continue;
        end

        for ti = 1:numel(termsOfInterest)
            r_o = matchTerm(fe_orig.Name, termsOfInterest{ti});
            r_w = matchTerm(fe_w.Name, termsOfInterest{ti});
            if isempty(r_o) || isempty(r_w), continue; end

            row = struct();
            row.Domain        = dxLabels{ci};
            row.Term          = termDisplayNames{ti};
            row.Beta_orig     = fe_orig.Estimate(r_o);
            row.SE_orig       = fe_orig.SE(r_o);
            row.p_orig        = fe_orig.pValue(r_o);
            row.Beta_winsor   = fe_w.Estimate(r_w);
            row.SE_winsor     = fe_w.SE(r_w);
            row.p_winsor      = fe_w.pValue(r_w);
            row.Beta_change_pct = 100 * (fe_w.Estimate(r_w) - fe_orig.Estimate(r_o)) / ...
                                  max(abs(fe_orig.Estimate(r_o)), 1e-6);
            % Flag if direction flipped or significance changed
            row.Direction_stable = sign(fe_orig.Estimate(r_o)) == sign(fe_w.Estimate(r_w));
            row.Sig_orig   = fe_orig.pValue(r_o) < 0.05;
            row.Sig_winsor = fe_w.pValue(r_w) < 0.05;
            row.Sig_stable = (row.Sig_orig == row.Sig_winsor);

            sensRows{end+1} = row; %#ok

            % Console output
            stabFlag = '';
            if ~row.Direction_stable, stabFlag = ' ** DIRECTION FLIP'; end
            if ~row.Sig_stable, stabFlag = [stabFlag ' * SIG CHANGE']; end %#ok
            fprintf('      %-28s  orig=%+.4f  winsor=%+.4f  (%+.1f%%)%s\n', ...
                termDisplayNames{ti}, row.Beta_orig, row.Beta_winsor, ...
                row.Beta_change_pct, stabFlag);
        end

    catch ME
        fprintf('    ERROR (%s): %s\n', dxLabels{ci}, ME.message);
    end
end

sensTab = struct2table([sensRows{:}]);
writetable(sensTab, fullfile(outputDir, 'Diagnostics_Sensitivity_Winsorised.csv'));
fprintf('\n  -> Diagnostics_Sensitivity_Winsorised.csv\n');

% -------------------------------------------------------------------------
% 9.  SENSITIVITY FIGURE: ORIGINAL vs WINSORISED COEFFICIENTS
% -------------------------------------------------------------------------
figSens = figure('Units','centimeters','Position',[1 1 18 14],'Color',COL.bg);
set(figSens,'DefaultAxesFontName','Helvetica','DefaultTextFontName','Helvetica');
tlS = tiledlayout(1, numel(dxVars), 'TileSpacing','compact', 'Padding','compact');
title(tlS, 'Sensitivity: Original vs Winsorised (1st–99th) Coefficients', ...
    'FontSize', 9, 'FontWeight', 'bold', 'Color', COL.accent);

for ci = 1:numel(dxVars)
    ax = nexttile; hold on;

    domRows = sensTab(strcmp(sensTab.Domain, dxLabels{ci}), :);
    nT = height(domRows);

    for ti = 1:nT
        yp = nT - ti + 1;
        b_o = domRows.Beta_orig(ti);
        b_w = domRows.Beta_winsor(ti);

        % Original (filled circle)
        origCol = [0.25 0.30 0.38];
        plot(b_o, yp + 0.1, 'o', 'Color', origCol, ...
            'MarkerFaceColor', origCol, 'MarkerSize', 5, 'LineWidth', 0.5);

        % Winsorised (open diamond)
        winCol = [0.72 0.38 0.20];
        plot(b_w, yp - 0.1, 'd', 'Color', winCol, ...
            'MarkerFaceColor', 'none', 'MarkerSize', 5.5, 'LineWidth', 1);

        % Connecting line to show shift
        plot([b_o b_w], [yp+0.1 yp-0.1], '-', 'Color', [0.7 0.7 0.7], ...
            'LineWidth', 0.6);

        % Flag instabilities
        if ~domRows.Direction_stable(ti)
            text(max(b_o, b_w) + 0.005, yp, 'FLIP', ...
                'FontSize', 5, 'Color', [0.8 0.2 0.2], 'FontWeight', 'bold');
        end
    end

    xline(0, '-', 'Color', [0.78 0.78 0.78], 'LineWidth', 0.5);

    set(ax, 'YTick', 1:nT, 'YTickLabel', flip(domRows.Term), ...
        'YLim', [0.3 nT+0.7], ...
        'FontSize', 6.5, 'TickDir', 'out', 'Box', 'off', ...
        'XGrid', 'off', 'YGrid', 'off');
    xlabel('\beta', 'FontSize', 7);
    if ci == 1
        ylabel('Fixed effect', 'FontSize', 7);
    end

    panelLetter = char('A' + ci - 1);
    title(ax, sprintf('\\bf%s\\rm  %s', panelLetter, dxLabels{ci}), ...
        'FontSize', 8, 'Color', COL.txt);

    % Legend on first panel
    if ci == 1
        hO = plot(nan, nan, 'o', 'Color', origCol, 'MarkerFaceColor', origCol, 'MarkerSize', 5);
        hW = plot(nan, nan, 'd', 'Color', winCol, 'MarkerFaceColor', 'none', 'MarkerSize', 5.5, 'LineWidth', 1);
        lg = legend([hO hW], {'Original','Winsorised'}, ...
            'Location','southeast', 'FontSize', 5.5, 'Box', 'off');
        lg.ItemTokenSize = [10 8];
    end

    formatAxis(ax, COL); hold off;
end

exportFigure(figSens, fullfile(outputDir, 'FigS2a_Sensitivity_Winsorised'));
fprintf('  -> FigS2a_Sensitivity_Winsorised\n');

% -------------------------------------------------------------------------
% 10.  INFLUENTIAL OBSERVATIONS: CONDITIONAL COOK'S DISTANCE
% -------------------------------------------------------------------------
% For each model, compute observation-level Cook's distance from the
% conditional (subject-specific) residuals and the hat matrix diagonal.
% Flag observations exceeding the conventional 4/n threshold.

fprintf('\n  === INFLUENTIAL OBSERVATIONS (Cook''s D) ===\n');

figCook = figure('Units','centimeters','Position',[1 1 18 10],'Color',COL.bg);
set(figCook,'DefaultAxesFontName','Helvetica','DefaultTextFontName','Helvetica');
tlC = tiledlayout(1, numel(dxVars), 'TileSpacing','compact', 'Padding','compact');
title(tlC, 'Cook''s Distance (conditional)', ...
    'FontSize', 9, 'FontWeight', 'bold', 'Color', COL.accent);

cookSummary = {};
for ci = 1:numel(dxVars)
    info = dxInfo{ci};
    if isempty(info) || ~isfield(info, 'mdl') || isempty(info.mdl), continue; end
    mdl = info.mdl;

    % Observation-level leverage and residuals.
    % NOTE: fixed-effects-only approximation to Cook's D — the hat matrix
    % ignores the random-effects design, and Pearson residuals are already
    % conditional on the random effects. Screening heuristic only.
    %
    % Leverages via thin QR: X = Q*R (orthonormal Q) ⇒ H = Q*Q', so
    %   hii = sum(Q.^2, 2). Avoids forming the n-by-n hat matrix and is
    %   numerically far more stable than inverting X'*X.
    %
    % Robustness: on some bifactor_g fits, designMatrix(mdl,'Fixed') returns
    % a row whose interaction column evaluates to NaN even though the LME
    % itself excluded that observation. NaN then propagates through qr →
    % every diagR entry is NaN → rank check returns 0 → every cookD is 0
    % → nothing is ever flagged influential. We therefore restrict the QR
    % (and Cook's D) to the rows that are finite end-to-end, mapping the
    % result back to the full-row index space at the end.
    X  = designMatrix(mdl, 'Fixed');
    if issparse(X), X = full(X); end
    eP = residuals(mdl, 'ResidualType', 'Pearson');

    goodRows = all(isfinite(X), 2) & isfinite(eP);
    nBad = sum(~goodRows);
    if nBad > 0
        fprintf('    NOTE (%s): %d/%d design rows non-finite — excluded from Cook''s D.\n', ...
            dxLabels{ci}, nBad, numel(goodRows));
    end

    Xc = X(goodRows, :);
    eC = eP(goodRows);

    [Q, R] = qr(Xc, 0);
    diagR  = abs(diag(R));
    maxDR  = max([diagR(:); eps]);
    rankX  = sum(diagR > max(size(Xc)) * eps(maxDR));
    if rankX < size(Xc, 2)
        fprintf('    WARNING (%s): fixed-effects design rank-deficient (%d/%d cols).\n', ...
            dxLabels{ci}, rankX, size(Xc, 2));
    end

    hii_c = sum(Q.^2, 2);
    p     = max(rankX, 1);
    n     = numel(eC);

    % Guard against pathological 1-hii values (leverage ~1 means deletion is undefined).
    denom_c = p .* (1 - hii_c).^2;
    cookD_c = (eC.^2 .* hii_c) ./ denom_c;
    cookD_c(~isfinite(cookD_c)) = 0;

    % Re-project to the full row-index space (length numel(goodRows)) so
    % the plot and §11 keepObs mask line up with the original observation
    % indices. Bad rows get hii=NaN / cookD=0 so they aren't flagged.
    nFull = numel(goodRows);
    hii   = nan(nFull, 1);  hii(goodRows)   = hii_c;
    cookD = zeros(nFull, 1); cookD(goodRows) = cookD_c;

    thresh = 4 / n;
    nInfluential = sum(cookD_c > thresh);
    pctInfluential = 100 * nInfluential / n;

    fprintf('    %s:  threshold=%.5f  n_influential=%d (%.1f%%)\n', ...
        dxLabels{ci}, thresh, nInfluential, pctInfluential);

    cs.Domain = dxLabels{ci};
    cs.N = n;                     % observations actually used for Cook's D
    cs.N_excluded_nan = nBad;     % design rows skipped due to NaN
    cs.Threshold = thresh;
    cs.N_influential = nInfluential;
    cs.Pct_influential = round(pctInfluential, 2);
    cs.Max_CookD = max(cookD_c);
    cs.Mean_CookD = mean(cookD_c);
    cs.Median_CookD = median(cookD_c);
    cookSummary{end+1} = cs; %#ok

    % Plot
    ax = nexttile; hold on;
    domCol = dxColors(ci,:);

    % Subsample non-influential for speed; plot ALL influential
    isInfl = cookD > thresh;
    nonInfl = find(~isInfl);
    if numel(nonInfl) > maxPts
        plotNon = nonInfl(randperm(numel(nonInfl), maxPts));
    else
        plotNon = nonInfl;
    end

    plot(plotNon, cookD(plotNon), '.', 'Color', [domCol 0.2], 'MarkerSize', 3);
    if any(isInfl)
        inflIdx = find(isInfl);
        plot(inflIdx, cookD(inflIdx), '.', 'Color', [0.8 0.2 0.2], 'MarkerSize', 5);
    end
    yline(thresh, '--', 'Color', COL.accent, 'LineWidth', 0.8);

    text(0.97, 0.95, sprintf('n_{infl}=%d\n(%.1f%%)', nInfluential, pctInfluential), ...
        'Units','normalized', 'FontSize', 6, 'HorizontalAlignment','right', ...
        'VerticalAlignment','top', 'Color', COL.txt);

    panelLetter = char('A' + ci - 1);
    title(ax, sprintf('\\bf%s\\rm  %s', panelLetter, dxLabels{ci}), ...
        'FontSize', 8, 'Color', COL.txt);
    xlabel('Observation', 'FontSize', 6); 
    if ci == 1, ylabel('Cook''s D', 'FontSize', 7); end
    set(ax, 'FontSize', 5.5, 'Box', 'off', 'TickDir', 'out');
    formatAxis(ax, COL); hold off;
end

exportFigure(figCook, fullfile(outputDir, 'FigS2b_Cooks_Distance'));
fprintf('  -> FigS2b_Cooks_Distance\n');

cookTab = struct2table([cookSummary{:}]);
writetable(cookTab, fullfile(outputDir, 'Diagnostics_CooksD_Summary.csv'));
fprintf('  -> Diagnostics_CooksD_Summary.csv\n');

% -------------------------------------------------------------------------
% 11.  SENSITIVITY RE-FIT: EXCLUDING HIGH-INFLUENCE OBSERVATIONS
% -------------------------------------------------------------------------
% Re-fit each model after dropping observations exceeding 4/n Cook's D,
% and compare to original coefficients.

fprintf('\n  === SENSITIVITY: EXCLUDING HIGH-COOK''s D OBSERVATIONS ===\n');

sensExcl = {};
for ci = 1:numel(dxVars)
    info = dxInfo{ci};
    if isempty(info) || ~isfield(info, 'mdl') || isempty(info.mdl), continue; end
    mdl = info.mdl;
    outVar = dxVars{ci};

    % Recompute Cook's D with NaN-robust QR (see §10 for rationale).
    X  = designMatrix(mdl, 'Fixed');
    if issparse(X), X = full(X); end
    eP = residuals(mdl, 'ResidualType', 'Pearson');

    goodRows = all(isfinite(X), 2) & isfinite(eP);
    nBad = sum(~goodRows);
    if nBad > 0
        fprintf('    NOTE (%s): %d/%d design rows non-finite — kept in refit, fitlme will re-exclude.\n', ...
            dxLabels{ci}, nBad, numel(goodRows));
    end

    Xc = X(goodRows, :);
    eC = eP(goodRows);

    [Q, R] = qr(Xc, 0);
    diagR  = abs(diag(R));
    maxDR  = max([diagR(:); eps]);
    rankX  = sum(diagR > max(size(Xc)) * eps(maxDR));
    hii_c  = sum(Q.^2, 2);
    p      = max(rankX, 1);
    n      = numel(eC);
    denom_c = p .* (1 - hii_c).^2;
    cookD_c = (eC.^2 .* hii_c) ./ denom_c;
    cookD_c(~isfinite(cookD_c)) = 0;
    thresh  = 4 / n;

    % Identify the rows in the ORIGINAL data table that were used.
    okAdj = okBase & ~isnan(lmeT.(outVar)) & ...
            ~isnan(lmeT.apathy) & ~isnan(lmeT.icb) & ~isnan(lmeT.time) & ...
            ~isnan(lmeT.age) & ~isnan(lmeT.educ) & ~isnan(lmeT.nmiss) & ...
            ~isnan(lmeT.duration) & ~isnan(lmeT.ledd) & ~isnan(lmeT.updrs3_score) & ...
            ~isundefined(lmeT.sex) & ~isnan(lmeT.depfactor);
    adjIdx = find(okAdj);

    % Drop only the high-Cook clean rows. Bad rows stay in lmeD_excl;
    % fitlme will exclude them again because they still carry NaN.
    if numel(adjIdx) ~= numel(goodRows)
        fprintf('    WARNING (%s): adjIdx (%d) ≠ design rows (%d); skipping Cook''s D refit.\n', ...
            dxLabels{ci}, numel(adjIdx), numel(goodRows));
        continue;
    end
    keepObs = true(numel(adjIdx), 1);
    keepObs(goodRows) = cookD_c <= thresh;
    lmeD_excl = lmeT(adjIdx(keepObs), :);

    nDropped = sum(~keepObs);
    fprintf('    %s:  dropped %d obs (%.1f%%), refitting on %d obs...\n', ...
        dxLabels{ci}, nDropped, 100*nDropped/n, height(lmeD_excl));

    try
        % Inherit the primary model's formula so the refit isolates the
        % effect of dropping Cook's-D-influential observations rather than
        % conflating it with a change in random-effects specification.
        fmla = char(mdl.Formula);
        mdl_excl = fitlme(lmeD_excl, fmla, 'FitMethod', 'REML');
        fe_excl  = coeffToTable(mdl_excl);

        fe_orig = dxFE{ci};
        if isempty(fe_orig)
            fprintf('    SKIP (%s): no original fixed-effects table to compare against.\n', dxLabels{ci});
            continue;
        end

        for ti = 1:numel(termsOfInterest)
            r_o = matchTerm(fe_orig.Name, termsOfInterest{ti});
            r_e = matchTerm(fe_excl.Name, termsOfInterest{ti});
            if isempty(r_o) || isempty(r_e), continue; end

            row = struct();
            row.Domain          = dxLabels{ci};
            row.Term            = termDisplayNames{ti};
            row.Beta_orig       = fe_orig.Estimate(r_o);
            row.p_orig          = fe_orig.pValue(r_o);
            row.Beta_excl       = fe_excl.Estimate(r_e);
            row.p_excl          = fe_excl.pValue(r_e);
            row.Beta_change_pct = 100 * (fe_excl.Estimate(r_e) - fe_orig.Estimate(r_o)) / ...
                                  max(abs(fe_orig.Estimate(r_o)), 1e-6);
            row.N_dropped       = nDropped;
            row.Direction_stable = sign(fe_orig.Estimate(r_o)) == sign(fe_excl.Estimate(r_e));
            row.Sig_orig   = fe_orig.pValue(r_o) < 0.05;
            row.Sig_excl   = fe_excl.pValue(r_e) < 0.05;
            row.Sig_stable = (row.Sig_orig == row.Sig_excl);

            sensExcl{end+1} = row; %#ok

            stabFlag = '';
            if ~row.Direction_stable, stabFlag = ' ** DIRECTION FLIP'; end
            if ~row.Sig_stable, stabFlag = [stabFlag ' * SIG CHANGE']; end %#ok
            fprintf('      %-28s  orig=%+.4f  excl=%+.4f  (%+.1f%%)%s\n', ...
                termDisplayNames{ti}, row.Beta_orig, row.Beta_excl, ...
                row.Beta_change_pct, stabFlag);
        end

    catch ME
        fprintf('    ERROR (%s): %s\n', dxLabels{ci}, ME.message);
    end
end

if isempty(sensExcl)
    fprintf('  WARNING: no Cook''s-D exclusion refits produced rows; CSV not written.\n');
    sensExclTab = table('Size', [0 12], ...
        'VariableTypes', {'string','string','double','double','double','double', ...
                          'double','double','logical','logical','logical','logical'}, ...
        'VariableNames', {'Domain','Term','Beta_orig','p_orig','Beta_excl','p_excl', ...
                          'Beta_change_pct','N_dropped','Direction_stable', ...
                          'Sig_orig','Sig_excl','Sig_stable'});
else
    sensExclTab = struct2table([sensExcl{:}]);
    writetable(sensExclTab, fullfile(outputDir, 'Diagnostics_Sensitivity_ExclCooksD.csv'));
    fprintf('  -> Diagnostics_Sensitivity_ExclCooksD.csv\n');
end

% -------------------------------------------------------------------------
% 12.  COMBINED SENSITIVITY SUMMARY  (console)
% -------------------------------------------------------------------------
fprintf('\n  ================================================================\n');
fprintf('  SENSITIVITY SUMMARY\n');
fprintf('  ================================================================\n');
fprintf('  Key question: Does the Apathy × Time decline survive on both\n');
fprintf('  global cognition outcomes?\n');
fprintf('    1. Apathy × Time on bifactor-g (primary)\n');
fprintf('    2. Apathy × Time on MoCA (secondary global outcome)\n\n');

keyChecks = {
    'Global Cognition (g)', 'Apathy \times Time';
    'MoCA',                  'Apathy \times Time'
};

for k = 1:size(keyChecks, 1)
    dom  = keyChecks{k, 1};
    term = keyChecks{k, 2};

    % Winsorised
    wRow = sensTab(strcmp(sensTab.Domain, dom) & strcmp(sensTab.Term, term), :);
    % Exclusion
    eRow = sensExclTab(strcmp(sensExclTab.Domain, dom) & strcmp(sensExclTab.Term, term), :);

    fprintf('  %s — %s\n', dom, term);
    if ~isempty(wRow)
        fprintf('    Original:    β=%+.4f  p=%.4f\n', wRow.Beta_orig, wRow.p_orig);
        fprintf('    Winsorised:  β=%+.4f  p=%.4f  (Δ=%+.1f%%  dir=%s  sig=%s)\n', ...
            wRow.Beta_winsor, wRow.p_winsor, wRow.Beta_change_pct, ...
            ternary(wRow.Direction_stable, 'STABLE', 'FLIPPED'), ...
            ternary(wRow.Sig_stable, 'STABLE', 'CHANGED'));
    end
    if ~isempty(eRow)
        fprintf('    Excl Cook D: β=%+.4f  p=%.4f  (Δ=%+.1f%%  dir=%s  sig=%s  n_drop=%d)\n', ...
            eRow.Beta_excl, eRow.p_excl, eRow.Beta_change_pct, ...
            ternary(eRow.Direction_stable, 'STABLE', 'FLIPPED'), ...
            ternary(eRow.Sig_stable, 'STABLE', 'CHANGED'), eRow.N_dropped);
    end
    fprintf('\n');
end

tElapsed = toc(tStart);
fprintf('\n[Diagnostics] Complete in %.1f seconds.\n', tElapsed);
%%
% =========================================================================
% LOCAL FUNCTIONS (must appear after all script-level executable code)
% =========================================================================
function out = ternary(cond, a, b)
    if cond, out = a; else, out = b; end
end

function v = altGet(fe, term, col)
    % Fetch fe.(col) for the row matching `term` (order-invariant).
    v = NaN;
    if isempty(fe) || ~ismember(col, fe.Properties.VariableNames), return; end
    nrm = @(s) strjoin(sort(strsplit(strtrim( ...
            regexprep(lower(s), '\s*[:\*]\s*', ':')), ':')), ':');
    r = find(cellfun(@(x) strcmp(nrm(x), nrm(term)), fe.Name), 1);
    if isempty(r), return; end
    v = fe.(col)(r);
end
%% ----------------------------
% -------- Helper Functions -------
%% ----------------------%%%%%%%%%%


function vals = safeNumeric(v)
    if isnumeric(v), vals = double(v);
    elseif iscell(v) || isstring(v), vals = str2double(string(v));
    elseif iscategorical(v), vals = double(v);
    elseif islogical(v), vals = double(v);
    else, vals = nan(numel(v), 1);
    end
end

function vals = safeCol(T, name)
    if ismember(name, T.Properties.VariableNames)
        vals = safeNumeric(T.(name));
    else
        vals = nan(height(T), 1);
    end
end



function s = pToStars(p)
    if isnan(p), s = ''; return; end
    if p < 0.001, s = '***';
    elseif p < 0.01, s = '**';
    elseif p < 0.05, s = '*';
    else, s = '';
    end
end





function ph = plotTrajectory(ax, years, mY, semY, nY, vYrs, col, showN)
    % Plot a single group trajectory with SEM band and sample size annotations.
    % Returns a single line handle for use in legends (fill is hidden).
    % showN: 'all' = annotate every year, 'endpoints' = first+last only, 'none'
    ph = gobjects(1);
    if ~any(vYrs), return; end
    yV = years(vYrs); mV = mY(vYrs); sV = semY(vYrs); nV = nY(vYrs);
    % SEM band (hidden from legend)
    fill(ax, [yV, fliplr(yV)], [mV(:)'+sV(:)', fliplr(mV(:)'-sV(:)')], ...
        col, 'FaceAlpha',0.12, 'EdgeColor','none', 'HandleVisibility','off');
    % Main line + markers (this is the legend entry)
    ph = plot(ax, yV, mV, '-o', 'Color',col, 'LineWidth',1.5, ...
        'MarkerFaceColor',col, 'MarkerSize',3.5, 'MarkerEdgeColor','w','MarkerEdgeColor',col);
    % Sample size annotations
    if nargin >= 8 && ~strcmp(showN, 'none')
        for yi = 1:numel(yV)
            if strcmp(showN, 'all') || ...
                    (strcmp(showN, 'endpoints') && (yi==1 || yi==numel(yV)))
                text(ax, yV(yi), mV(yi)-sV(yi), sprintf('n=%d',nV(yi)), ...
                    'FontSize',3.8, 'Color',col*0.6, 'HorizontalAlignment','center', ...
                    'VerticalAlignment','top');
            end
        end
    end
end





function exportFigure(fig, filepath)
    set(fig,'PaperPositionMode','auto','Renderer','painters');
    print(fig, [filepath '.pdf'], '-dpdf', '-bestfit');
    print(fig, [filepath '.png'], '-dpng', '-r600');
    close(fig);
end

function formatAxis(ax, COL, fs)
    % fs (optional) sets the axes font size; defaults to 7 so existing
    % callers are unchanged. Pass a larger value to enlarge tick labels.
    % formatAxis runs after any per-axes 'FontSize' set() call, so this is
    % the value that ultimately takes effect on the tick labels.
    if nargin < 3 || isempty(fs), fs = 7; end
    set(ax,'FontSize',fs,'Box','off','TickDir','out','XGrid','off','YGrid','off', ...
        'GridColor',COL.grid,'GridAlpha',0.4);
end


function hLine = plotCPTraj(ax, Tg, rowFilter, xVar, yVar, col, nBins, sigma, lineStyle)
%PLOTCPTRAJ Plot a Gaussian-smoothed trajectory with within-subject SE.
%   Implements the same logic as conditionalPlot (StandardBinning +
%   WithinSubjectError) directly, avoiding external dependencies (sq,
%   nancat, removeIdenticals).
%
%   Pipeline:
%     1. Bin each subject's Y by integer Year (standard binning)
%     2. Within-subject error correction (Loftus & Masson 1994):
%        subtract each subject's grand mean across timepoints
%     3. Correction factor sqrt(k/(k-1)) for the lost DF
%     4. Gaussian-smooth M and SE, render on target axes
    
    if nargin < 9, lineStyle = '-'; end
    hLine = gobjects(1);
    
    % --- Extract data for this group ---
    tbl = Tg(rowFilter, :);
    xAll = tbl.(xVar);
    yAll = double(safeCol(tbl, yVar));
    patAll = tbl.PATNO;
    
    % Unique timepoints (bins) — for integer Year this is just 0,1,2,...
    uYears = unique(xAll(~isnan(xAll)));
    uYears = sort(uYears);
    nT = numel(uYears);
    if nT < 2, return; end
    
    subjects = unique(patAll);
    nSubj = numel(subjects);
    if nSubj < 5, return; end
    
    % --- Standard binning: mean Y at each Year for each subject ---
    Yb = nan(nSubj, nT);   % (subject, timepoint)
    for si = 1:nSubj
        sIdx = patAll == subjects(si);
        for ti = 1:nT
            vals = yAll(sIdx & xAll == uYears(ti));
            vals = vals(~isnan(vals));
            if ~isempty(vals)
                Yb(si, ti) = nanmean(vals);
            end
        end
    end
    
    % --- Within-subject error correction (Loftus & Masson 1994) ---
    % For each subject, subtract their grand mean across all timepoints.
    % This removes between-subject variance in overall level, leaving
    % only variance in trajectory *shape* — appropriate for comparing
    % time courses within a repeated-measures design.
    subjGrandMean = nanmean(Yb, 2);                          % (nSubj x 1)
    Yb_demeaned = bsxfun(@minus, Yb, subjGrandMean);
    
    % Morey (2008) correction factor: sqrt(k/(k-1)) where k = n timepoints
    % with data for at least one subject
    k = sum(~all(isnan(Yb), 1));
    corrFactor = sqrt(k / max(1, k - 1));
    
    % Group mean and within-subject SE at each timepoint
    M  = nanmean(Yb, 1);
    nPerBin = sum(~isnan(Yb_demeaned), 1);
    SE = nanstd(Yb_demeaned, [], 1) ./ sqrt(max(1, nPerBin)) * corrFactor;
    Xm = uYears(:)';
    
    % Remove empty bins
    valid = ~isnan(M) & nPerBin >= 5;
    M = M(valid); SE = SE(valid); Xm = Xm(valid);
    if numel(Xm) < 2, return; end
    
    % --- Gaussian smoothing ---
    yrsD = linspace(min(Xm), max(Xm), 200);
    mSmooth  = gaussSmooth(yrsD, interp1(Xm, M,  yrsD, 'pchip'), sigma);
    seSmooth = gaussSmooth(yrsD, interp1(Xm, SE, yrsD, 'pchip'), sigma);
    
    % --- Plot on target axes ---
    axes(ax); hold on; %#ok<LAXES>
    
    % Within-subject SE band
    fill([yrsD fliplr(yrsD)], [(mSmooth+seSmooth) fliplr(mSmooth-seSmooth)], ...
        col, 'FaceAlpha', 0.12, 'EdgeColor', 'none');
    
    % Smoothed trajectory
    hLine = plot(yrsD, mSmooth, lineStyle, 'Color', col, 'LineWidth', 1.5);
    
    % Raw bin means (data points)
    plot(Xm, M, '.', 'Color', col, 'MarkerSize', 6);
end


function yS = gaussSmooth(x, y, sigma)
    yS = nan(size(y));
    for i = 1:numel(x)
        w = exp(-0.5 * ((x - x(i)) / sigma).^2);
        w = w / sum(w);
        yS(i) = sum(w .* y);
    end
end

function s = ternaryStr(c, a, b)
    if c, s = a; else, s = b; end
end

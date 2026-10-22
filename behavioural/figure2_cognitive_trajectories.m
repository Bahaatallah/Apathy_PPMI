% Role: Primary cognition LMEs (g and MoCA), sensitivity models (pattern mixture, winsorised
%     predictors, conditional Cook's distance and refit excluding influential observations),
%     sample-size table and Figure 2.
% Run: on the workspace left by behavioural/cohort_phenotyping.m, loaded in a fresh MATLAB
%     session (this script standardises columns of lmeT in place); outputs go to outputDir.
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

% -------------------------------------------------------------------------
% 0.  SETTINGS
% -------------------------------------------------------------------------
globalVar      = 'bifactor_g';     % global cognition outcome (bifactor g)
globalVarLabel = 'factor g';     % display name (used in panel B)


cogVars = {globalVar};
nCols = numel(cogVars);

cogLabels = {'Global Cognition (g)'};

years     = 0:15;

% Top-row visualisation mode
%   'group4' – 4-way phenotype (Neither | Apathy | ICB | Both)
%   'apathy' – 2-way split by clinical apathy flag  (Low vs High)
%   'icb'    – 2-way split by clinical ICB flag     (Low vs High)
% Splits reuse the flags already encoded in group4 so the cut-points
% match the ones used for the group4 analysis.
vizMode = 'group4';

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

% Trajectory grouping shared by the Figure 2 panels
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

% -------------------------------------------------------------------------
% 4.  RUN ADJUSTED LMEs  (store everything for forest plot + Supplementary Table S5)
% -------------------------------------------------------------------------
allFE        = cell(1, nCols);    % fixed-effect tables (one per domain)
allModelInfo = cell(1, nCols);    % struct per domain for Supplementary Table S5
tableNums    = {'Table8'};

for ci = 1:nCols
    outVar = cogVars{ci};

    % --- Filter to complete cases for adjusted model ---
    % NOTE: every fixed-effect predictor (incl. the apathy/icb/time terms
    % that drive the 3-way interaction) and every covariate must be NaN-
    % screened here so that height(lmeD) == mdl.NumObservations downstream.
    % If we miss one, fitlme silently drops those rows but
    % designMatrix(mdl,'Fixed') still returns them as NaN — which then
    % poisons the QR step of the Cook's distance section (every leverage → NaN → 0 → no
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
        fmla = sprintf(['%s ~ apathy*icb*time + age + sex + educ + nmiss + ' ...
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

        % === Store model-level info for Supplementary Table S5 ===
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
% 4b.  PATTERN-MIXTURE MODEL (sensitivity to informative dropout)
% -------------------------------------------------------------------------
%   Adds each patient's last-observed follow-up year (z-scored) and its
%   interaction with time. Tests whether trajectories differ by dropout
%   timing; if apathy terms are stable, findings are robust to informative
%   dropout. The model CSV and a comparison table with the primary model
%   are written to outputDir.
% -------------------------------------------------------------------------
fprintf('[Fig2] Pattern-mixture sensitivity model...\n');

altFE = struct();
altFE.patmix  = cell(1, nCols);

for ci = 1:nCols
    outVar = cogVars{ci};

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

        fmla_pm = sprintf(['%s ~ apathy*icb*time + age + sex + educ + nmiss + ' ...
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

        fprintf('  %-24s pattern-mix : %4d pts, median lastFU=%.1f y (%.0f%% >=5y, %.0f%% >=10y)\n', ...
            cogLabels{ci}, attr.N_patients, attr.Median_lastFU, ...
            attr.Pct_ge_5y, attr.Pct_ge_10y);
    catch ME
        fprintf('  ERROR pat-mix (%s): %s\n', cogLabels{ci}, ME.message);
        altFE.patmix{ci} = table();
    end
end

% --- Side-by-side β/p comparison of the primary and pattern-mixture models ---
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

%% ========================================================================
%  PRIMARY FIGURE 2  (new layout)
%  A. bifactor_g trajectory by group4 (or vizMode split)
%  B. Forest plot with paired markers for bifactor_g and MoCA;
%
%  NOTE: bifactor_g and MoCA are BOTH treated as primary global outcomes
%  in panel B and use raw p-values.
% =========================================================================
fprintf('[Fig2 Primary] Trajectory + paired forest...\n');

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
    fmla_m = ['moca ~ apathy*icb*time + age + sex + educ + nmiss + duration + ' ...
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
% Step 2: locate the global cognition outcome (used by panels A and B)
% -------------------------------------------------------------------------
gIdx = find(strcmp(cogVars, globalVar), 1);
if isempty(gIdx), gIdx = 1; end

% -------------------------------------------------------------------------
% Step 4: build the figure
%   Outer 1x2 layout. Only panels A and B in the primary figure.
%     Panel A: 2x3 nested grid -> main trajectory spans 2 rows x 2 cols
%              (left); two smaller sub-panels stack on the right
%              (apathy split top, ICB split bottom).
%     Panel B: standalone forest plot (globalVar + MoCA, paired markers).
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
nT     = numel(termsOfInterest);

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
% Shared setup for the sensitivity refits
% =========================================================================
tStart = tic;

% Clean up any workspace variables that shadow built-in functions
clear fitted residuals dd diagData;

% Diagnostics focus: bifactor_g (primary) + MoCA (secondary global outcome).
% MoCA is fitted above as mdl_moca / feMoca; we wrap it
% into the same {info, fe-table} shape used by allModelInfo / allFE so the
% loops below can treat the two outcomes uniformly.
dxVars   = {globalVar, 'moca'};
dxLabels = {'Global Cognition (g)', 'MoCA'};

dxInfo = cell(1, numel(dxVars));
dxFE   = cell(1, numel(dxVars));

% (1) bifactor_g comes straight from the primary model outputs.
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
% 10.  INFLUENTIAL OBSERVATIONS: CONDITIONAL COOK'S DISTANCE
% -------------------------------------------------------------------------
% For each model, compute observation-level Cook's distance from the
% conditional (subject-specific) residuals and the hat matrix diagonal.
% Flag observations exceeding the conventional 4/n threshold.

fprintf('\n  === INFLUENTIAL OBSERVATIONS (Cook''s D) ===\n');

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
end

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

    % Recompute Cook's D with NaN-robust QR (see the Cook's distance section for the rationale).
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


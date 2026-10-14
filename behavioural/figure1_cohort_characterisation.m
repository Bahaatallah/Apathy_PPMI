% Role: Figure 1 renderer (Okabe-Ito palette); writes the Table 1 and Supplementary Table S1 CSVs
%     (Table1_Baseline_Demographics.csv, SuppTable1_PostHoc_Pairwise.csv) and the dimensional-view
%     figure.
% Run: on the workspace left by behavioural/cohort_phenotyping.m, loaded in a fresh MATLAB
%     session; outputs go to outputDir. Set FIG_SUPPRESS_ID_EXPORTS = true before running to skip
%     the two participant-identifier exports.
% Manuscript item: Figure 1; Table 1; Supplementary Figure S2; Supplementary Table S1.
%% Generates manuscript Figure 1 (Fig1_SampleCharacterisation_AllPD), Supplementary Figure S2
%  (Fig1_DimensionalView), Table 1 and Supplementary Table S1. Text is sized for print at the
%  170 mm inclusion width. Optional flags: FIG_SUPPRESS_ID_EXPORTS = true skips the two
%  participant-identifier exports; FIG_FONT_AUDIT = true writes a CSV of printed font sizes.

%% ========================================================================
%  RESULTS SECTION 1 — SAMPLE CHARACTERISATION
%  Table 1  : Baseline demographics (KW / chi-square, with effect sizes)
%  Supp T1  : Post-hoc pairwise comparisons (Dunn / Fisher, Bonferroni)
%  Figure 1 : (a) Prevalence bar,  (b) Clinical radar,  (c) Sankey transitions
% =========================================================================
% --- Publication colour palette (Okabe-Ito, colourblind-safe) -------------
% Okabe-Ito values, shared with figure2_cognitive_trajectories.m, figure5_genetic_moderation.m and fc_cognitive_slope_mediation.m.
COL.neither  = [0.00 0.60 0.50];   % Okabe-Ito bluish green
COL.apathy   = [0.90 0.60 0.00];   % Okabe-Ito orange      
COL.icb      = [0.00 0.45 0.70];   % Okabe-Ito blue        
COL.overlap  = [0.45 0.20 0.50];   % deep purple           
COL.groups   = [COL.neither; COL.apathy; COL.icb; COL.overlap];
COL.bg       = [1 1 1];            % white
COL.grid     = [0.85 0.85 0.85];
COL.txt      = [0.15 0.15 0.15];
COL.accent   = [0.10 0.10 0.10];
%% ========================================================================
%  RESULTS SECTION 1 — SAMPLE CHARACTERISATION
%  Table 1  : Baseline demographics (KW / chi-square, with effect sizes)
%  Supp T1  : Post-hoc pairwise comparisons (Dunn / Fisher, Bonferroni)
%  Figure 1 : (a) Prevalence bar,  (b) Clinical radar,  (c) Sankey transitions
% =========================================================================

%% ===  TABLE 1 + SUPP TABLE 1  ==========================================
%  Run both together so KW stats are computed once and reused for post-hocs.
fprintf('\n[Table 1 / Supp Table 1] Demographics + post-hoc comparisons...\n');

% --- Derived binary: on dopaminergic therapy (LEDD>0) --------------------
%  PPMI is a de novo cohort: ~77% of participants have LEDD == 0 at baseline,
%  so median[IQR] collapses to "0 [0-0]" for the larger groups and hides the
%  real group differences. Show treatment prevalence + LEDD mean+/-SD instead.
%
%  Encoding mirrors SEX (where 0 = female and pctN reports "Sex (% female)"):
%  on_therapy = 0 means ON therapy (LEDD > 0); pctN will then display "(%)"
%  as "% on therapy" by counting zeros.
bl.on_therapy = double(bl.LEDD == 0);
bl.on_therapy(isnan(bl.LEDD)) = NaN;

% --- Variable definitions ------------------------------------------------
t1Defs = {
    'age_at_visit',        'Age (years)',               'cont', 'meansd';
    'SEX',                 'Sex (% female)',             'cat',  'pctN';
    'EDUCYRS',             'Education (years)',          'cont', 'meansd';
    'duration_yrs',        'Disease duration (years)',   'cont', 'mediqr';
    'moca',                'MoCA total (0–30)',          'cont', 'meansd';
    'apathy_composite',    'Apathy composite (0–7)',     'cont', 'mediqr';
    'depression_factor',   'Depression factor (GDS-9)',  'cont', 'mediqr';
    'gds',                 'GDS-15 total',               'cont', 'mediqr';
    'quip',                'QUIP-S ICB domains',         'cont', 'mediqr';
    'stai_trait',          'STAI trait anxiety',          'cont', 'meansd';
    'updrs3_score',        'MDS-UPDRS Part III (motor)', 'cont', 'meansd';
    'updrs2_score',        'MDS-UPDRS Part II (ADL)',    'cont', 'meansd';
    'MSEADLG',             'Schwab–England ADL (%)',     'cont', 'meansd';
    'NHY',                 'Hoehn & Yahr stage',         'cont', 'mediqr';
    'ess',                 'ESS sleepiness (0–24)',      'cont', 'meansd';
    'rem',                 'RBDSQ total',                'cont', 'meansd';
    'scopa',               'SCOPA-AUT total',            'cont', 'meansd';
    'MIA_CAUDATE_BILAT',   'DaT caudate (bilateral)',    'cont', 'meansd';
    'MIA_PUTAMEN_BILAT',   'DaT putamen (bilateral)',    'cont', 'meansd';
    'on_therapy',          'On dopaminergic therapy (%)','cat',  'pctN';
    'LEDD',                'LEDD (mg, mean ± SD)',       'cont', 'meansd';
};

nVars  = size(t1Defs, 1);
nG     = 4;
grpN   = arrayfun(@(g) sum(bl.group4 == g), 1:nG);

% --- Pre-allocate outputs ------------------------------------------------
headerT1 = ['Variable', ...
    arrayfun(@(g) sprintf('%s (n=%d)', groupLabelsShort{g}, grpN(g)), ...
             1:nG, 'UniformOutput', false), ...
    'p', 'Effect size'];
summaryRows = cell(nVars, nG + 3);

pairLabels = {'1v2','1v3','1v4','2v3','2v4','3v4'};
pairIdx    = [1 2; 1 3; 1 4; 2 3; 2 4; 3 4];
nPairs     = size(pairIdx, 1);
headerSup  = ['Variable', pairLabels];
suppRows   = cell(nVars, 1 + nPairs);

% --- Single pass: Table 1 + Supp Table 1 --------------------------------
for v = 1:nVars
    vname  = t1Defs{v,1};
    vlabel = t1Defs{v,2};
    vtype  = t1Defs{v,3};
    vfmt   = t1Defs{v,4};

    allV = safeCol(bl, vname);
    grpV = bl.group4;
    ok   = ~isnan(allV) & ~isnan(grpV);

    summaryRows{v,1} = vlabel;
    suppRows{v,1}    = vlabel;

    % --- Per-group summaries (Table 1) -----------------------------------
    for g = 1:nG
        vals = allV(grpV == g);
        vals = vals(~isnan(vals));
        n    = numel(vals);
        if n == 0
            summaryRows{v, g+1} = '—';
        else
            switch vfmt
                case 'meansd'
                    summaryRows{v, g+1} = sprintf('%.1f +/- %.1f', mean(vals), std(vals));
                case 'mediqr'
                    summaryRows{v, g+1} = sprintf('%.1f [%.1f–%.1f]', ...
                        median(vals), prctile(vals,25), prctile(vals,75));
                case 'pctN'
                    nEvent = sum(vals == 0);
                    summaryRows{v, g+1} = sprintf('%.1f%% (%d/%d)', ...
                        100*nEvent/n, nEvent, n);
            end
        end
    end

    % --- Skip if too few observations ------------------------------------
    if sum(ok) <= 10
        summaryRows{v, nG+2} = '—';
        summaryRows{v, nG+3} = '—';
        for pp = 1:nPairs, suppRows{v, pp+1} = '—'; end
        continue
    end

    % --- Omnibus + post-hoc in one block ---------------------------------
    if strcmp(vtype, 'cat')
        % Omnibus: chi-square
        observed = zeros(2, nG);
        for g = 1:nG
            gv = allV(ok & grpV == g);
            observed(1, g) = sum(gv == 0);
            observed(2, g) = sum(gv ~= 0);
        end
        totalN   = sum(observed(:));
        rowSum   = sum(observed, 2);
        colSum   = sum(observed, 1);
        expected = (rowSum * colSum) / totalN;

        % MATLAB's crosstab returns [tbl, chi2, p, labels] -- name carefully.
        [~, chi2val, pOmni] = crosstab(allV(ok), grpV(ok));
        k = min(size(observed));
        cramV = sqrt(chi2val / (totalN * (k - 1)));
        summaryRows{v, nG+3} = sprintf('V=%.2f', cramV);

        if any(expected(:) < 5)
            summaryRows{v, nG+3} = [summaryRows{v, nG+3}, ' (low n)'];
        end

        % Post-hoc: pairwise chi-square, Bonferroni
        for pp = 1:nPairs
            g1 = pairIdx(pp,1); g2 = pairIdx(pp,2);
            mask = ok & (grpV == g1 | grpV == g2);
            if sum(mask) < 5
                suppRows{v, pp+1} = '—';
                continue
            end
            % crosstab returns [tbl, chi2, p, labels]; we want the third out.
            [~, ~, pRaw] = crosstab(allV(mask), grpV(mask));
            if isnan(pRaw), pRaw = 1; end
            suppRows{v, pp+1} = formatPwith(min(1, pRaw * nPairs));
        end

    else
        % Omnibus: Kruskal-Wallis (run ONCE, reuse stats for post-hoc)
        [pOmni, tbl, stats] = kruskalwallis(allV(ok), grpV(ok), 'off');

        % Epsilon-squared: H statistic is in tbl{2,5}
        H     = tbl{2,5};
        N     = sum(ok);
        eta2H = max(0, (H - nG + 1) / (N - nG));
        summaryRows{v, nG+3} = sprintf('eta2=%.3f', eta2H);

        % Post-hoc: Dunn via multcompare (Bonferroni), reusing stats
        mc = multcompare(stats, 'CType','bonferroni', 'Display','off');
        for pp = 1:nPairs
            g1 = pairIdx(pp,1); g2 = pairIdx(pp,2);
            row = (mc(:,1)==g1 & mc(:,2)==g2) | (mc(:,1)==g2 & mc(:,2)==g1);
            if any(row)
                suppRows{v, pp+1} = formatPwith(mc(row, 6));
            else
                suppRows{v, pp+1} = '—';
            end
        end
    end

    % Format omnibus p
    if pOmni < 0.001
        summaryRows{v, nG+2} = '< .001';
    else
        summaryRows{v, nG+2} = sprintf('%.3f', pOmni);
    end
end

% --- Export ---------------------------------------------------------------
summaryT1 = cell2table(summaryRows, 'VariableNames', headerT1);
disp(summaryT1);
writetable(summaryT1, fullfile(outputDir, 'Table1_Baseline_Demographics.csv'));
fprintf('    -> Table1_Baseline_Demographics.csv\n');

suppT1 = cell2table(suppRows, 'VariableNames', headerSup);
disp(suppT1);
writetable(suppT1, fullfile(outputDir, 'SuppTable1_PostHoc_Pairwise.csv'));
fprintf('    -> SuppTable1_PostHoc_Pairwise.csv\n\n');


%% ===  FIGURE 1: THREE-PANEL HORIZONTAL  =================================
%  (a) Baseline prevalence  (b) Clinical radar  (c) Sankey transitions
%
%  Built TWICE, using one shared helper (`buildFigure1Panel`) — DRY:
%    - Full PD cohort  (overview of phenotype distribution)
%    - Incident cohort (<1 yr disease duration at baseline; cleaner natural-
%      history view of phenotype emergence and transitions)
%
%  Aesthetic conventions follow Figure 4:
%    * tiledlayout, formatAxis, exportFigure, COL palette
%    * Helvetica + bold axes/titles, "(a)/(b)/(c)" baked into title
%    * Stats annotations on white background ([1 1 1 0.75])
%
%  Bug fix vs. previous version: Sankey flow X-coordinate is now LINEAR
%  (cubic smoothstep applied to Y only). The previous all-PD block applied
%  the smoothstep to BOTH X and Y, which warped the time axis.
% =========================================================================
fprintf('[Figure 1] Building three-panel horizontal: prevalence + radar + Sankey...\n');

% ---- ALL-PD figure (cohort overview) ------------------------------------
[fig1all, sk_all] = buildFigure1Panel(bl, T, isPD, COL, ...
    groupLabels, groupLabelsShort, ...
    'cohort','all', 'sankeyYears',[0 2 4], 'minComplete',30);
if exist('FIG_FONT_AUDIT', 'var') && isequal(FIG_FONT_AUDIT, true)   % legibility QA (printed at 170 mm)
    auditFigureText(fig1all, 17/26, fullfile(outputDir, 'FontAudit_Fig1_SampleCharacterisation_AllPD.csv'));
end
exportFigure(fig1all, fullfile(outputDir, 'Fig1_SampleCharacterisation_AllPD'));
fprintf('    -> Fig1_SampleCharacterisation_AllPD\n');

% ---- INCIDENT cohort figure (primary natural-history transitions) -------
[fig1inc, sk_inc] = buildFigure1Panel(bl, T, isPD, COL, ...
    groupLabels, groupLabelsShort, ...
    'cohort','incident', 'incidentDurationYrs',1.0, ...
    'sankeyYears',[0 2 4], 'minComplete',10);
exportFigure(fig1inc, fullfile(outputDir, 'Fig1_SampleCharacterisation_Incident'));
fprintf('    -> Fig1_SampleCharacterisation_Incident\n\n');

% ---- Expose incident-cohort variables for downstream transition stats ---
% (matches the original behaviour where the second build overwrote vars)
patGroupComplete = sk_inc.patGroupComplete;
sankeyPatIDs     = sk_inc.sankeyPatIDs;
sankeyYears      = sk_inc.sankeyYears;
complete         = sk_inc.complete;
nComplete        = sk_inc.nComplete;
nTrans           = numel(sankeyYears) - 1;
nG               = 4;

% ---- DIMENSIONAL companion figure (continuous, no thresholding) --------
%   Honest visualisation of the apathy and ICB scores that the
%   statistical models actually use as continuous predictors. Designed
%   as a Supplementary / "honesty check" companion to the categorical
%   4-group Figure 1 above (which is descriptive only).
fig1d = buildFigure1Continuous(bl, T, isPD, COL);
exportFigure(fig1d, fullfile(outputDir, 'Fig1_DimensionalView'));
fprintf('    -> Fig1_DimensionalView\n\n');

%% ===  TRANSITION STATISTICS  ============================================
% (Transition stats code remains unchanged but will naturally run on the newly filtered subset)
%% ===  TRANSITION STATISTICS  ============================================
%  Transition matrices, stability rates, McNemar-Bowker symmetry test,
%  chi-square tests of independence, and per-group retention/emergence.
fprintf('[Transition Stats] Analysing group transitions...\n');

if nComplete > 30

    transLabels = groupLabelsShort(1:nG);

    for ti = 1:nTrans
        yrFrom = sankeyYears(ti);
        yrTo   = sankeyYears(ti+1);
        fprintf('\n--- Transition: Year %d -> Year %d ---\n', yrFrom, yrTo);

        g_from = patGroupComplete(:, ti);
        g_to   = patGroupComplete(:, ti+1);

        % ----- 1. Transition count matrix (rows = from, cols = to) -------
        transMat = zeros(nG, nG);
        for g1 = 1:nG
            for g2 = 1:nG
                transMat(g1, g2) = sum(g_from == g1 & g_to == g2);
            end
        end

        % ----- 2. Row-normalised transition probability matrix -----------
        rowSums  = sum(transMat, 2);
        rowSums(rowSums == 0) = 1;   % guard against empty rows
        transPct = 100 * transMat ./ rowSums;

        % Display
        fprintf('  Count matrix:\n');
        countT = array2table(transMat, 'RowNames', transLabels, 'VariableNames', transLabels);
        disp(countT);

        fprintf('  Row %% (transition probabilities):\n');
        pctT = array2table(round(transPct, 1), 'RowNames', transLabels, 'VariableNames', transLabels);
        disp(pctT);

        % ----- 3. Per-group stability (% remaining in same group) --------
        stability = zeros(nG, 1);
        for g = 1:nG
            nFrom = sum(g_from == g);
            if nFrom > 0
                stability(g) = 100 * transMat(g,g) / nFrom;
            end
        end
        overallStability = 100 * sum(diag(transMat)) / nComplete;

        fprintf('  Stability rates (%%remaining in same group):\n');
        for g = 1:nG
            fprintf('    %s: %.1f%% (%d/%d)\n', transLabels{g}, stability(g), ...
                transMat(g,g), sum(g_from==g));
        end
        fprintf('    Overall: %.1f%% (%d/%d)\n', overallStability, ...
            sum(diag(transMat)), nComplete);

        % ----- 4. Chi-square test of independence ------------------------
        %  H0: group at time ti+1 is independent of group at time ti
        %  crosstab returns [tbl, chi2, p, labels] -> chi2 is 2nd, p is 3rd
        [~, statIndep, pIndep] = crosstab(g_from, g_to);
        fprintf('  Chi-square test of independence: chi2=%.2f, p=%s\n', ...
            statIndep, formatPwith(pIndep));

        % ----- 5. McNemar-Bowker test of symmetry ------------------------
        %  H0: transition matrix is symmetric (P(i->j) = P(j->i))
        %  Uses the upper triangle of (transMat - transMat') differences
        bowkerStat = 0;
        bowkerDF   = 0;
        for g1 = 1:(nG-1)
            for g2 = (g1+1):nG
                nij = transMat(g1, g2);
                nji = transMat(g2, g1);
                if (nij + nji) > 0
                    bowkerStat = bowkerStat + (nij - nji)^2 / (nij + nji);
                    bowkerDF   = bowkerDF + 1;
                end
            end
        end
        if bowkerDF > 0
            pBowker = 1 - chi2cdf(bowkerStat, bowkerDF);
            fprintf('  McNemar-Bowker symmetry: chi2=%.2f, df=%d, p=%s\n', ...
                bowkerStat, bowkerDF, formatPwith(pBowker));
        else
            pBowker = NaN;
            fprintf('  McNemar-Bowker symmetry: not computable (no off-diagonal flow)\n');
        end

        % ----- 6. Marginal homogeneity (Stuart-Maxwell) -----------------
        %  H0: row marginals = column marginals (group sizes unchanged)
        margFrom = sum(transMat, 2)';   % 1 x nG
        margTo   = sum(transMat, 1);    % 1 x nG
        d = margTo(1:end-1) - margFrom(1:end-1);  % first nG-1 differences
        if any(d ~= 0)
            % Variance-covariance of d under H0
            V = zeros(nG-1);
            for a = 1:(nG-1)
                for bb = 1:(nG-1)
                    if a == bb
                        V(a,bb) = margFrom(a) + margTo(a) - 2*transMat(a,a);
                    else
                        V(a,bb) = -(transMat(a,bb) + transMat(bb,a));
                    end
                end
            end
            if rcond(V) > eps
                smStat = d * (V \ d');
                smDF   = nG - 1;
                pSM    = 1 - chi2cdf(smStat, smDF);
                fprintf('  Stuart-Maxwell marginal homogeneity: chi2=%.2f, df=%d, p=%s\n', ...
                    smStat, smDF, formatPwith(pSM));
            else
                pSM = NaN;
                fprintf('  Stuart-Maxwell: V matrix singular, test not computable\n');
            end
        else
            pSM = NaN;
            fprintf('  Stuart-Maxwell: marginals identical, no test needed\n');
        end

        % ----- 7. Export transition matrices as CSVs ---------------------
        % Count matrix
        countTbl = array2table(transMat, 'RowNames', transLabels, 'VariableNames', transLabels);
        fname = sprintf('Transition_Count_Yr%d_to_Yr%d.csv', yrFrom, yrTo);
        writetable(countTbl, fullfile(outputDir, fname), 'WriteRowNames', true);
        fprintf('    -> %s\n', fname);

        % Probability matrix
        pctTbl = array2table(round(transPct, 1), 'RowNames', transLabels, 'VariableNames', transLabels);
        fname = sprintf('Transition_Pct_Yr%d_to_Yr%d.csv', yrFrom, yrTo);
        writetable(pctTbl, fullfile(outputDir, fname), 'WriteRowNames', true);
        fprintf('    -> %s\n', fname);
    end

    % ----- 8. Summary stats table across all intervals -------------------
    %  One row per transition interval: N, overall stability, chi2, Bowker p, SM p
    fprintf('\n  Compiling transition summary table...\n');
    summTransRows = cell(nTrans, 7);
    for ti = 1:nTrans
        yrFrom = sankeyYears(ti); yrTo = sankeyYears(ti+1);
        g_from = patGroupComplete(:, ti);
        g_to   = patGroupComplete(:, ti+1);

        transMat_s = zeros(nG, nG);
        for g1 = 1:nG
            for g2 = 1:nG
                transMat_s(g1,g2) = sum(g_from==g1 & g_to==g2);
            end
        end

        stab = 100 * sum(diag(transMat_s)) / nComplete;
        % crosstab returns [tbl, chi2, p, labels] -> chi2 is 2nd, p is 3rd
        [~, chiInd, pInd] = crosstab(g_from, g_to);

        % Bowker
        bS = 0; bDF = 0;
        for g1 = 1:(nG-1)
            for g2 = (g1+1):nG
                nij = transMat_s(g1,g2); nji = transMat_s(g2,g1);
                if (nij+nji)>0, bS = bS + (nij-nji)^2/(nij+nji); bDF = bDF+1; end
            end
        end
        if bDF > 0, pB = 1 - chi2cdf(bS, bDF); else, pB = NaN; end

        % Stuart-Maxwell
        mF = sum(transMat_s,2)'; mT = sum(transMat_s,1);
        dd = mT(1:end-1) - mF(1:end-1);
        Vs = zeros(nG-1);
        for a=1:(nG-1), for bb=1:(nG-1)
            if a==bb, Vs(a,bb) = mF(a)+mT(a)-2*transMat_s(a,a);
            else,     Vs(a,bb) = -(transMat_s(a,bb)+transMat_s(bb,a)); end
        end, end
        if rcond(Vs)>eps, smS=dd*(Vs\dd'); pSMs=1-chi2cdf(smS,nG-1);
        else, pSMs=NaN; end

        summTransRows(ti,:) = {
            sprintf('Yr %d–%d', yrFrom, yrTo), ...
            nComplete, ...
            sprintf('%.1f%%', stab), ...
            sprintf('%.1f', chiInd), formatPwith(pInd), ...
            formatPwith(pB), formatPwith(pSMs)
        };
    end

    summTransT = cell2table(summTransRows, 'VariableNames', ...
        {'Interval','N','Stability','Chi2_indep','p_indep','p_Bowker','p_StuartMaxwell'});
    disp(summTransT);
    writetable(summTransT, fullfile(outputDir, 'Transition_Summary_Stats.csv'));
    fprintf('    -> Transition_Summary_Stats.csv\n\n');

else
    fprintf('    Fewer than 30 complete cases — transition stats skipped.\n\n');
end


%% ===  HELPER: format p-values  ==========================================
function s = formatPwith(p)
    if isnan(p)
        s = '—';
    elseif p < 0.001
        s = '< .001***';
    elseif p < 0.01
        s = sprintf('%.3f**', p);
    elseif p < 0.05
        s = sprintf('%.3f*', p);
    else
        s = sprintf('%.3f', p);
    end
end

%%
%% ===  IDENTIFY PERSISTENT SUBJECTS  ======================================
fprintf('\n[Persistent Subjects] Extracting PATNOs for stable phenotypes...\n');

% 1. Extract the exact PATNOs that match the rows in patGroupComplete
% (The 'complete' logical array was created during the Sankey prep)
if exist('complete', 'var') && exist('sankeyPatIDs', 'var')
    completePatIDs = sankeyPatIDs(complete);
else
    error('Required variables (complete, sankeyPatIDs) not found in workspace.');
end

% 2. Define group indices based on your standard 1-4 mapping
% (Adjust these numbers if your internal group coding differs)
GRP_NEITHER = 1;
GRP_APATHY  = 2;
GRP_ICB     = 3;
GRP_OVERLAP = 4;

% 3. Find boolean vectors for patients who were strictly in one group across ALL columns
is_persist_neither = all(patGroupComplete == GRP_NEITHER, 2);
is_persist_apathy  = all(patGroupComplete == GRP_APATHY, 2);
is_persist_icb     = all(patGroupComplete == GRP_ICB, 2);
is_persist_overlap = all(patGroupComplete == GRP_OVERLAP, 2);

% 4. Extract the actual PATNOs
patnos_persist_neither = completePatIDs(is_persist_neither);
patnos_persist_apathy  = completePatIDs(is_persist_apathy);
patnos_persist_icb     = completePatIDs(is_persist_icb);
patnos_persist_overlap = completePatIDs(is_persist_overlap);

% 5. Count the number of persistent subjects per group
n_persist_neither = numel(patnos_persist_neither);
n_persist_apathy  = numel(patnos_persist_apathy);
n_persist_icb     = numel(patnos_persist_icb);
n_persist_overlap = numel(patnos_persist_overlap);

% 6. Print the summary to the command window
fprintf('---------------------------------------------------\n');
fprintf('Persistent Subgroup Counts (across all timepoints):\n');
fprintf('---------------------------------------------------\n');
fprintf('  Persistent Neither: n = %d\n', n_persist_neither);
fprintf('  Persistent Apathy:  n = %d\n', n_persist_apathy);
fprintf('  Persistent ICB:     n = %d\n', n_persist_icb);
fprintf('  Persistent Overlap: n = %d\n', n_persist_overlap);
fprintf('---------------------------------------------------\n');

% 7. Optional: Save to a structured variable or CSV for easy exporting later
PersistentCohort = struct(...
    'Neither', patnos_persist_neither, ...
    'Apathy',  patnos_persist_apathy, ...
    'ICB',     patnos_persist_icb, ...
    'Overlap', patnos_persist_overlap);


%%
%% ===  EXTRACT BASELINE DEMOGRAPHICS BY GROUP  ============================
% this section and the next write participant IDs (PATNO) to
% outputDir and the log; a calling driver can skip both with
% FIG_SUPPRESS_ID_EXPORTS = true. Default (flag absent): unchanged.
if ~(exist('FIG_SUPPRESS_ID_EXPORTS', 'var') && isequal(FIG_SUPPRESS_ID_EXPORTS, true))
fprintf('\n[Baseline Extraction] Identifying subjects and demographics per group...\n');

% 1. Define your exact column names here (update these if necessary)
col_patno = 'PATNO';
col_group = 'group4';
col_age   = 'age_at_visit';     % e.g., 'age', 'age_at_bl', 'ENROLL_AGE'
col_sex   = 'SEX';     % e.g., 'sex', 'GENDER'
col_site  = 'SITE';    % e.g., 'site', 'CEN'

% 2. Ensure we are using the baseline table 
% (Assuming 'bl' from your Figure 1 code contains the baseline PD cohort)
if ~exist('bl', 'var')
    error('Baseline table ''bl'' not found in workspace.');
end

% 3. Verify all requested columns exist in the table
reqCols = {col_patno, col_group, col_age, col_sex, col_site};
missingCols = reqCols(~ismember(reqCols, bl.Properties.VariableNames));
if ~isempty(missingCols)
    error('The following required columns are missing from ''bl'': %s', strjoin(missingCols, ', '));
end

% 4. Extract the relevant columns
T_bl_demog = bl(:, reqCols);

% 5. Filter out anyone without a valid group assignment (1 through 4)
isValidGroup = T_bl_demog.(col_group) >= 1 & T_bl_demog.(col_group) <= 4;
T_bl_demog = T_bl_demog(isValidGroup, :);

% 6. Map the numeric groups to descriptive text labels for readability
groupNamesMap = {'Neither', 'Apathy', 'ICB', 'Overlap'};
mappedLabels = cell(height(T_bl_demog), 1);
for i = 1:height(T_bl_demog)
    mappedLabels{i} = groupNamesMap{T_bl_demog.(col_group)(i)};
end
T_bl_demog.GroupLabel = mappedLabels;

% 7. Sort the table so it's grouped cleanly 1 to 4
T_bl_demog = sortrows(T_bl_demog, col_group);

% 8. Display a quick summary of what was extracted
fprintf('Extracted %d baseline subjects.\n', height(T_bl_demog));
for g = 1:4
    n_grp = sum(T_bl_demog.(col_group) == g);
    fprintf('  %s: n = %d\n', groupNamesMap{g}, n_grp);
end

% 9. Export to CSV for easy sharing or viewing in Excel
outName = 'Baseline_Subjects_Demographics.csv';
writetable(T_bl_demog, fullfile(outputDir, outName));
fprintf('    -> Saved to %s\n', outName);

% 10. Optional: Create separate tables if you want to access them in MATLAB
BL_Neither = T_bl_demog(T_bl_demog.(col_group) == 1, :);
BL_Apathy  = T_bl_demog(T_bl_demog.(col_group) == 2, :);
BL_ICB     = T_bl_demog(T_bl_demog.(col_group) == 3, :);
BL_Overlap = T_bl_demog(T_bl_demog.(col_group) == 4, :);

%%
%% ===  EXPORT COMMA-SEPARATED PATNO LISTS  ============================
fprintf('\n[Comma-Separated Lists] Generating PATNO strings for each baseline group...\n');

% Ensure col_patno and col_group are still defined from the previous step
% (e.g., col_patno = 'PATNO'; col_group = 'group4';)
if ~exist('T_bl_demog', 'var')
    error('Please run the baseline extraction code first so T_bl_demog is available.');
end

groupNamesMap = {'Neither', 'Apathy', 'ICB', 'Overlap'};

% Create a text file to save the full strings
outTextFile = fullfile(outputDir, 'Baseline_PATNO_Lists.txt');
fid = fopen(outTextFile, 'w');

if fid == -1
    error('Could not open file for writing. Check your outputDir path.');
end

for g = 1:4
    % Extract the PATNOs for the current group
    current_patnos = T_bl_demog.(col_patno)(T_bl_demog.(col_group) == g);
    
    % Convert the IDs to strings and join them with a comma and a space
    % (The ' symbol transposes it to ensure it joins correctly whether it's a row or column vector)
    patno_str = strjoin(string(current_patnos'), ', ');
    
    % Print a preview to the MATLAB command window
    fprintf('\n--- %s (n=%d) ---\n', groupNamesMap{g}, numel(current_patnos));
    if strlength(patno_str) > 100
        fprintf('%s ... [preview truncated]\n', extractBefore(patno_str, 100));
    else
        fprintf('%s\n', patno_str);
    end
    
    % Write the full, untruncated string to the text file
    fprintf(fid, '--- %s (n=%d) ---\n', groupNamesMap{g}, numel(current_patnos));
    fprintf(fid, '%s\n\n', patno_str);
end

fclose(fid);
fprintf('\n    -> Full comma-separated lists successfully saved to:\n       %s\n\n', outTextFile);
else
    fprintf('\n[Baseline Extraction / PATNO lists] skipped (FIG_SUPPRESS_ID_EXPORTS = true): no participant IDs written\n');
end   % FIG_SUPPRESS_ID_EXPORTS guard
%% ----------------------------
% -------- Helper Functions -------
%% ----------------------%%%%%%%%%%
function feT = lme_coeff_to_table(mdl)
    % Convert fitlme Coefficients (dataset or table) to a standard table
    % with a 'Name' column containing coefficient names.
    % Works across all MATLAB versions regardless of whether Coefficients
    % returns a dataset or a table.
    raw = mdl.Coefficients;
    if isa(raw, 'dataset')
        feT = dataset2table(raw);
    else
        feT = raw;  % already a table
    end
    % Coefficient names may be in RowNames or ObsNames depending on type
    if ~ismember('Name', feT.Properties.VariableNames)
        rn = feT.Properties.RowNames;
        if isempty(rn)
            try rn = feT.Properties.ObsNames; catch, end
        end
        if isempty(rn)
            rn = strcat('Term', arrayfun(@num2str, (1:height(feT))', 'Uni', false));
        end
        feT.Name = rn;
        feT = [feT(:,end), feT(:,1:end-1)];  % move Name to first column
    end
end

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

function s = nanSumIfComplete(mat, minValid)
    nValid = sum(~isnan(mat), 2);
    s = sum(mat, 2, 'omitnan');
    s(nValid < minValid) = NaN;
end

function alpha = cronbach_alpha(X)
    [~, k] = size(X);
    itemVar = var(X, 0, 1);
    totalVar = var(sum(X, 2));
    alpha = (k / (k - 1)) * (1 - sum(itemVar) / totalVar);
end

function pAdj = fdr_bh(pvals)
    pAdj = nan(size(pvals));
    valid = ~isnan(pvals);
    m = sum(valid);
    if m == 0, return; end
    pv = pvals(valid);
    [sorted, sI] = sort(pv);
    bh = sorted .* m ./ (1:m)';
    for i = m-1:-1:1, bh(i) = min(bh(i), bh(i+1)); end
    bh = min(bh, 1);
    adj = nan(m, 1); adj(sI) = bh;
    pAdj(valid) = adj;
end

function s = pToStars(p)
    if isnan(p), s = ''; return; end
    if p < 0.001, s = '***';
    elseif p < 0.01, s = '**';
    elseif p < 0.05, s = '*';
    else, s = '';
    end
end

function out = ternary(cond, a, b)
    if cond, out = a; else, out = b; end
end

function [mY, semY, nY, validYrs] = longitudinalStats(T, varName, yearVar, years, grp, grpCol)
    % Generic longitudinal summary: works with any group column.
    % grpCol: column vector of group labels (default: T.group4)
    nYrs = numel(years);
    mY = nan(nYrs,1); semY = nan(nYrs,1); nY = nan(nYrs,1);
    yrCol = safeCol(T, yearVar);
    if nargin < 6, grpCol = T.group4; end
    for yi = 1:nYrs
        idx = grpCol == grp & yrCol == years(yi);
        vals = safeCol(T, varName);
        vals = vals(idx); vals = vals(~isnan(vals));
        nY(yi) = numel(vals);
        if numel(vals) >= 5
            mY(yi) = mean(vals);
            semY(yi) = std(vals) / sqrt(numel(vals));
        end
    end
    validYrs = ~isnan(mY);
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

function addAttritionTable(ax, years, nMat, grpLabels, grpCols, yFrac)
    % Add a sample-size attrition table below the plot.
    % nMat: nGroups x nYears matrix of sample sizes
    % yFrac: vertical position as fraction of axes (default -0.22)
    if nargin < 6, yFrac = -0.22; end
    nG = size(nMat, 1); nY = numel(years);
    % Year headers
    for yi = 1:nY
        text(ax, years(yi), yFrac, sprintf('Yr %d',years(yi)), ...
            'Units','data','FontSize',4, 'Color',[0.4 0.4 0.4], ...
            'HorizontalAlignment','center','VerticalAlignment','top', ...
            'FontWeight','bold');
    end
end

function [slopes, grps] = computePatientSlopes(T, varName, grpColName)
    % Compute per-patient linear slopes and baseline group assignments.
    patients = unique(T.PATNO);
    slopes = nan(numel(patients), 1);
    grps   = nan(numel(patients), 1);
    for pi = 1:numel(patients)
        pIdx = T.PATNO == patients(pi);
        pT = T(pIdx, :);
        yrs = safeCol(pT, 'YEAR');
        vals = safeCol(pT, varName);
        ok = ~isnan(yrs) & ~isnan(vals);
        if sum(ok) >= 3
            b = polyfit(yrs(ok), vals(ok), 1);
            slopes(pi) = b(1);
        end
        blRow = yrs == 0 & ~isnan(safeCol(pT, grpColName));
        if any(blRow)
            grpVals = safeCol(pT, grpColName);
            grps(pi) = grpVals(find(blRow, 1));
        end
    end
end

function cmap = divergingColormap(n)
    if nargin < 1, n = 256; end
    half = floor(n/2);
    r1 = linspace(0.02,1,half)'; g1 = linspace(0.19,1,half)'; b1 = linspace(0.58,1,half)';
    r2 = linspace(1,0.70,n-half)'; g2 = linspace(1,0.10,n-half)'; b2 = linspace(1,0.12,n-half)';
    cmap = [r1,g1,b1; r2,g2,b2];
end

function setPublicationDefaults()
    set(0,'DefaultAxesFontName','Helvetica');
    set(0,'DefaultAxesFontSize',8);
    set(0,'DefaultTextFontName','Helvetica');
    set(0,'DefaultFigureColor',[1 1 1]);
    set(0,'DefaultAxesLineWidth',0.6);
    set(0,'DefaultAxesTickDir','out');
    set(0,'DefaultAxesTickLength',[0.015 0.015]);
    set(0,'DefaultAxesXColor',[0.2 0.2 0.2]);
    set(0,'DefaultAxesYColor',[0.2 0.2 0.2]);
    set(0,'DefaultLineLineWidth',1.2);
end

function exportFigure(fig, filepath)
    set(fig,'PaperPositionMode','auto','Renderer','painters');
    print(fig, [filepath '.pdf'], '-dpdf', '-bestfit');
    print(fig, [filepath '.png'], '-dpng', '-r600');
    close(fig);
end

function formatAxis(ax, COL)
    set(ax,'FontSize',7,'Box','off','TickDir','out','XGrid','on','YGrid','on', ...
        'GridColor',COL.grid,'GridAlpha',0.4);
end

function addRegressionLine(x, y, lineCol, txtCol)
    p = polyfit(x, y, 1);
    xf = linspace(min(x), max(x), 100);
    plot(xf, polyval(p, xf), '-', 'Color',lineCol, 'LineWidth',1.5);
    [rho, pv] = corr(x, y, 'Type','Spearman');
    pStr = ternary(pv<0.001, '< .001', sprintf('= %.3f',pv));
    text(0.05, 0.08, sprintf('\\rho = %.2f, p %s',rho,pStr), ...
        'Units','normalized','FontSize',6.5,'Color',txtCol);
end
% =========================================================================
%  FIGURE 1 HELPERS — three-panel orchestration + improved sub-panels
%  Follows the grammar of Figure 4: tiledlayout, formatAxis, COL palette,
%  Helvetica/bold axes, "(a) Title" panel-letter convention, exportFigure.
% =========================================================================
function [fig, sankeyData] = buildFigure1Panel(bl, T, isPD, COL, groupLabels, groupLabelsShort, varargin)
    % BUILDFIGURE1PANEL  Four-panel Figure 1 in a 2x2 grid (Fig 4 grammar).
    %   Row 1 (cross-sectional pair):
    %     (a) Baseline prevalence lollipop
    %     (b) Cleveland dot plot of clinical profiles (replaces radar)
    %         — z-scored group means on a shared x-axis, one dot per group
    %         per variable. Uses MATLAB native zscore().
    %   Row 2 (longitudinal companion pair):
    %     (c) Transition-probability heatmap (Y0->Y2 and Y2->Y4 blocks)
    %     (d) Sankey transitions across selected years
    %
    %   Name-value options:
    %     'cohort'              'all' | 'incident'   default 'all'
    %     'incidentDurationYrs' (1,1) double         default 1.0
    %     'sankeyYears'         (1,:) double         default [0 2 4]
    %     'minComplete'         (1,1) double         default 30
    %     'panelTitle'          char/string          default ''
    p = inputParser;
    p.addParameter('cohort','all', @(x) ismember(x,{'all','incident'}));
    p.addParameter('incidentDurationYrs', 1.0);
    p.addParameter('sankeyYears', [0 2 4]);
    p.addParameter('minComplete', 30);
    p.addParameter('panelTitle','');
    p.parse(varargin{:});
    opt = p.Results;

    % --- Figure setup (Figure 4 grammar) ---------------------------------
    fig = figure('Units','centimeters','Position',[1 1 26 18.5],'Color',COL.bg);   % was 26 x 16 (room for print-size text)
    set(fig, ...
        'DefaultAxesFontName','Helvetica','DefaultTextFontName','Helvetica', ...
        'DefaultAxesFontWeight','bold','DefaultTextFontWeight','normal', ...
        'DefaultAxesTitleFontWeight','bold');
    % explicit axes positions in PRINTED cm (17.0 x 12.1 at the
    % 170 mm inclusion width), same 2 x 2 arrangement. The former
    % tiledlayout(2,12,'compact') shrank panels a and c to make room for the
    % longer panel-b labels, so the print-size text in c collided.
    PW = 17.0; PH = 18.5 * 17 / 26;
    posF = @(l, b, w, h) [l/PW, b/PH, w/PW, h/PH];
    if ~isempty(opt.panelTitle)
        sgtitle(fig, opt.panelTitle, 'FontSize', f1pt(8.5), 'FontWeight','bold', 'Color', COL.txt, 'Interpreter','tex');
    end

    % --- Pre-compute the Sankey data ONCE (shared by panels c and d) --
    sankeyData = computeSankeyData(T, isPD, opt);

    % ---- Panel (a) Prevalence lollipop  (row 1, left) -----------------
    pA = posF(1.20, 7.45, 3.80, 4.00);
    ax_a = axes(fig, 'Position', pA); set(ax_a,'Color',COL.bg); hold(ax_a,'on');
    drawPrevalenceLollipop(ax_a, bl, COL, groupLabelsShort);
    formatAxis(ax_a, COL);
    set(ax_a,'XGrid','off','YGrid','off','FontSize',f1pt(6.5));
    ax_a.XLabel.FontSize = f1pt(7);

    % ---- Panel (b) Cleveland dot plot of clinical profiles (row 1, right) -
    % legend uses the short group names of panels a, c, d
    % (was groupLabels: "Apathy only", "ICB only", "Apathy + ICB").
    pB = posF(9.00, 7.45, 7.85, 4.00);
    ax_b = axes(fig, 'Position', pB); set(ax_b,'Color',COL.bg); hold(ax_b,'on');
    drawClevelandDotPlot(ax_b, bl, COL, groupLabelsShort);
    ax_b.Position = pB;

    % ---- Panel (c) Transition-probability heatmap  (row 2, left) ------
    pC = posF(1.45, 1.30, 5.45, 4.40);
    ax_c = axes(fig, 'Position', pC); set(ax_c,'Color',COL.bg); hold(ax_c,'on');
    drawTransitionHeatmap(ax_c, sankeyData.patGroupComplete, ...
        sankeyData.sankeyYears, COL, groupLabelsShort);
    cbC = findobj(fig, 'Type', 'colorbar');
    if ~isempty(cbC), cbC(1).Position = posF(7.00, 1.30, 0.17, 4.40); end
    ax_c.Position = pC;

    % ---- Panel (d) Sankey transitions  (row 2, right) -----------------
    % nodes are all participants classified at each visit
    % and flows all participants classified at both consecutive visits, so
    % the title says that (it said "n=550", the completers of panel c); the n
    % of every column and flow is printed inside the panel.
    pD = posF(8.60, 0.25, 8.30, 5.10);
    ax_d = axes(fig, 'Position', pD); set(ax_d,'Color',COL.bg); hold(ax_d,'on');
    drawSankeyPanelFromData(ax_d, sankeyData, COL, groupLabelsShort, opt);
    if strcmp(opt.cohort,'incident')
        ttlD = 'Incident Transitions'; sttl = 'incident participants classified at consecutive visits';
    else
        ttlD = 'Group Transitions'; sttl = 'all participants classified at consecutive visits';
    end

    % ---- Panel titles: overlay axes in printed cm, left-aligned at each
    % column start on a common baseline per row.
    axT = axes(fig, 'Position', [0 0 1 1], 'Visible', 'off', 'XLim', [0 PW], 'YLim', [0 PH], 'HitTest', 'off');
    tt = {0.10, 11.80, 'a', 'Baseline Prevalence'; ...
          5.40, 11.80, 'b', 'Clinical Profiles  (group mean z-score; right = worse)'; ...
          0.10,  6.00, 'c', sprintf('Transition Probabilities  (n = %d)', sankeyData.nComplete); ...   % the n of the completers used here (was on d)
          8.30,  6.00, 'd', ttlD};
    for k = 1:size(tt, 1)
        text(axT, tt{k,1}, tt{k,2}, ['\bf' tt{k,3} '\rm  ' tt{k,4}], 'FontSize', f1pt(8), 'Color', COL.txt, ...
            'HorizontalAlignment', 'left', 'VerticalAlignment', 'baseline');
    end
    text(axT, 8.30 + 0.33, 5.62, sttl, 'FontSize', f1pt(7), 'Color', COL.txt, ...
        'HorizontalAlignment', 'left', 'VerticalAlignment', 'baseline');   % panel d subtitle

    drawnow;
end


function drawPrevalenceLollipop(ax, bl, COL, groupLabelsShort)
    % Lollipop chart of baseline group prevalence — much cleaner than a bar.
    nBL    = height(bl);
    counts = arrayfun(@(g) sum(bl.group4 == g), 1:4);
    pcts   = 100 * counts / nBL;
    [~, ord] = sort(pcts, 'ascend');     % small at bottom, big at top

    yPos = 1:4;
    for i = 1:4
        g = ord(i);  y = yPos(i);  col = COL.groups(g,:);
        % Stem
        plot(ax, [0 pcts(g)], [y y], '-', 'Color', col, 'LineWidth', 2.4);
        % Lollipop head (white edge for crispness)
        plot(ax, pcts(g), y, 'o', 'MarkerSize', 9, ...
            'MarkerFaceColor', col, 'MarkerEdgeColor', 'w', 'LineWidth', 1);
        % Right-side numeric annotation
        text(ax, pcts(g) + max(pcts)*0.04, y, ...
            sprintf('%.1f%%   n=%d', pcts(g), counts(g)), ...
            'FontSize', f1pt(6.5), 'FontWeight','normal','Color', COL.txt, ...
            'VerticalAlignment','middle','HorizontalAlignment','left', ...
            'BackgroundColor',[1 1 1 0.75], 'Margin', 1);
        fprintf('    [Fig1a] %-8s %.1f%%   n=%d\n', groupLabelsShort{g}, pcts(g), counts(g));   % printed values to the log
    end
    set(ax, 'YTick', yPos, 'YTickLabel', groupLabelsShort(ord), ...
        'YDir','normal', 'YLim', [0.4 4.6], ...
        'XLim', [0 max(pcts)*1.65]);   % was 1.45; room for the print-size labels
    xlabel(ax, 'Prevalence (%)', 'FontSize', f1pt(7));
end


function drawClevelandDotPlot(ax, bl, COL, groupLabels)
    % Cleveland dot plot of clinical profiles by group.
    %   Replaces the radar — one dot per group per variable, on a shared
    %   z-score x-axis. Far less spaghetti than 4 overlapping polygons.
    %
    %   Standardisation: per-variable cohort z-score using MATLAB's native
    %   zscore(). Because zscore returns NaN if any element is NaN, we
    %   apply it column-wise on the non-NaN rows only, leaving NaN in
    %   place for missing observations.
    %
    %   Inverted variables (MoCA, ADL where higher = better) are sign-
    %   flipped after z-scoring so that "right of the reference line =
    %   worse" is true for every row, preserving the radar's directional
    %   convention without forcing a polygon.
    radarDefs = {
        'updrs3_score',      'Motor',       false;
        'moca',              'MoCA',        true;
        'apathy_composite',  'Apathy',      false;
        'depression_factor', 'Depression',  false;
        'stai_trait',        'Anxiety',     false;
        'quip',              'Impulsivity', false;
        'ess',               'Sleepiness',  false;
        'rem',               'REM Sleep',   false;
        'scopa',             'Autonomic',   false;
        'MSEADLG',           'ADL',         true;
    };
    nVar       = size(radarDefs, 1);
    varNames   = radarDefs(:,1);
    varLabels  = radarDefs(:,2)';
    isInverted = [radarDefs{:,3}];
    % display labels say which rows are sign-reversed
    % and that ADL is the Schwab-England scale (labels only).
    varLabels{strcmp(varNames, 'MSEADLG')} = ['Schwab' char(8211) 'England ADL (reversed)'];
    varLabels{strcmp(varNames, 'moca')}    = 'MoCA (reversed)';

    allData = bl{:, varNames};

    % --- Native MATLAB zscore() per column (NaN-safe) -----------------
    % zscore(X) returns NaN-laden output if X contains any NaN, so we
    % apply it only on the non-NaN rows of each column.
    zAll = nan(size(allData));
    for j = 1:nVar
        ok = ~isnan(allData(:, j));
        if sum(ok) >= 2
            zAll(ok, j) = zscore(allData(ok, j));
        end
    end

    % --- Per-group z-score profile (group mean of cohort z-scores) ----
    zPerGrp = nan(4, nVar);
    for g = 1:4
        idx = (bl.group4 == g);
        if ~any(idx), continue; end
        zg = mean(zAll(idx, :), 1, 'omitnan');
        zg(isInverted) = -zg(isInverted);    % flip so right = worse on every row
        zPerGrp(g, :) = zg;
    end

    % --- Layout setup --------------------------------------------------
    yPos = 1:nVar;
    xMag = max([abs(zPerGrp(:)); 0.5], [], 'omitnan') * 1.15;
    xMag = ceil(xMag * 4) / 4;            % round up to nearest 0.25

    hold(ax, 'on');

    % Subtle alternating row backgrounds for legibility
    for v = 1:nVar
        if mod(v, 2) == 0
            patch(ax, [-xMag xMag xMag -xMag], ...
                [v-0.45 v-0.45 v+0.45 v+0.45], ...
                [0.96 0.96 0.97], 'EdgeColor','none');
        end
    end

    % Reference line at z = 0 (cohort mean)
    plot(ax, [0 0], [0.4 nVar+0.6], '-', ...
        'Color', [0.45 0.45 0.45], 'LineWidth', 1.0);
    % Light dashed guides at z = +/- 1
    for vz = [-1, 1]
        plot(ax, [vz vz], [0.4 nVar+0.6], '--', ...
            'Color', [0.85 0.85 0.85], 'LineWidth', 0.5);
    end

    % Within-row connecting line spanning the group-mean range
    for v = 1:nVar
        zv = zPerGrp(:, v);
        ok = ~isnan(zv);
        if sum(ok) >= 2
            plot(ax, [min(zv(ok)) max(zv(ok))], [v v], '-', ...
                'Color', [0.78 0.78 0.82], 'LineWidth', 1.4);
        end
    end

    % --- Conditional vertical lanes: jitter only where dots overlap ------
    %   On most rows the four group means are well separated, so dots stay on
    %   the row centre line (yOff = 0). Only where two or more groups have
    %   (near-)identical z — e.g. Neither/Apathy on Impulsivity (both
    %   QUIP == 0), Apathy/Overlap on Motor, Apathy/ICB on MoCA — are the
    %   colliding groups fanned into separate vertical lanes so none is
    %   hidden. Groups not part of a collision on a given row stay centred.
    %
    %   Collisions are found per row by 1-D single-linkage clustering of the
    %   present groups along z: consecutive sorted means closer than
    %   collideThresh join the same cluster; clusters of size >= 2 are
    %   spread symmetrically about the centre line in steps of laneStep.
    collideThresh = 0.12;   % |Δz| below which two dots are treated as overlapping
    laneStep      = 0.22;   % vertical gap between lanes within a collision cluster
    yOff = zeros(4, nVar);
    for v = 1:nVar
        zv = zPerGrp(:, v);
        ok = find(~isnan(zv));
        if numel(ok) < 2, continue; end
        [zs, sIdx] = sort(zv(ok));                    % ascending z among present groups
        grpSorted  = ok(sIdx);
        clusterId  = cumsum([0; diff(zs) > collideThresh]) + 1;
        for c = 1:max(clusterId)
            members = grpSorted(clusterId == c);
            m = numel(members);
            if m >= 2
                lanes = (((0:m-1) - (m-1)/2) * laneStep)';   % symmetric about centre
                yOff(members, v) = lanes;
            end
        end
    end

    % Dots per group per variable (opaque; lane offset applied only on rows
    % where the group collides with another — yOff is 0 everywhere else).
    ph = gobjects(4, 1);
    for g = 1:4
        zv = zPerGrp(g, :);
        ok = ~isnan(zv);
        if any(ok)
            ph(g) = plot(ax, zv(ok), yPos(ok) + yOff(g, ok), 'o', ...
                'MarkerSize', 6.5, 'MarkerFaceColor', COL.groups(g,:), ...
                'MarkerEdgeColor', 'w', 'LineWidth', 0.8);
        end
    end

    % --- Axis cosmetics ------------------------------------------------
    set(ax, 'YTick', yPos, 'YTickLabel', varLabels, 'YDir','reverse', ...
        'XLim', [-xMag xMag], 'YLim', [0.4 nVar+0.6], ...
        'FontSize', f1pt(6.5), 'FontWeight','bold', 'Box','off', 'TickDir','out');
    xlabel(ax, 'Group mean (cohort z-score)', 'FontSize', f1pt(7), 'FontWeight', 'bold');

    % --- Legend (bottom, horizontal, 2 columns; original style) -------
    validIdx = find(arrayfun(@(h) isgraphics(h), ph));
    if ~isempty(validIdx)
        legNames = groupLabels(validIdx);     % no sample-size suffix (caller passes the short names)
        % legend in the empty lower-left corner of the plot (no
        % group mean lies below z = -0.7 in any row); it was southoutside,
        % where the print-size text would collide with the panel below.
        lg = legend(ax, ph(validIdx), legNames, ...
            'Location', 'southwest', 'Orientation', 'vertical', ...
            'FontSize', f1pt(6.5), 'Box', 'off', 'NumColumns', 1);
        lg.ItemTokenSize = [15 15];
    end
end


function drawRadarPanelV2(ax, bl, COL, groupLabels)  %#ok<DEFNU>
    % NOTE: superseded by drawClevelandDotPlot; kept for reference / fallback.
    % Polished radar of clinical profiles by group.
    %   *Z-SCORE normalisation* (centred at total cohort mean):
    %       z = (x - cohort_mu) / cohort_sigma
    %       r = 0.5 + z/4    (so z = -2 -> 0, z = 0 -> 0.5, z = +2 -> 1)
    %       For variables where higher = better (MoCA, ADL), the sign of z
    %       is flipped first so that 0 = best and 1 = worst on every axis,
    %       preserving the original directional convention.
    %   This puts the cohort mean on the middle ring and lets group means
    %   move *inward* (better than average) or *outward* (worse than average),
    %   greatly increasing inter-group visual contrast vs. min-max scaling.
    %
    %   Title + legend layout matches the ORIGINAL radarPlotPub:
    %       text(0, 1.35, 'Clinical Profiles')   - bold title
    %       text(0, 1.25, '(0 = Best, 1 = Worst)') - subtitle
    %       legend southoutside, horizontal, 2 cols, fontsize 8
    %       panel-letter '(b)' as a textbox at the axis top-left
    radarDefs = {
        'updrs3_score',      'Motor',       false;
        'moca',              'MoCA',        true;
        'apathy_composite',  'Apathy',      false;
        'depression_factor', 'Depression',  false;
        'stai_trait',        'Anxiety',     false;
        'quip',              'Impulsivity', false;
        'ess',               'Sleepiness',  false;
        'rem',               'REM Sleep',   false;
        'scopa',             'Autonomic',   false;
        'MSEADLG',           'ADL',         true;
    };
    nAxes      = size(radarDefs, 1);
    axLabels   = radarDefs(:,2)';
    isInverted = [radarDefs{:,3}];
    varNames   = radarDefs(:,1);

    allData = bl{:, varNames};

    % --- Cohort z-score parameters (per variable) ----------------------
    mu_c    = mean(allData, 1, 'omitnan');
    sigma_c = std(allData, 0, 1, 'omitnan');
    sigma_c(sigma_c == 0 | isnan(sigma_c)) = 1;   % guard against zero-var

    plotData = nan(4, nAxes);
    nPerGrp  = zeros(4,1);
    for g = 1:4
        idx = (bl.group4 == g);
        if ~any(idx), continue; end
        m  = mean(allData(idx,:), 1, 'omitnan');
        z  = (m - mu_c) ./ sigma_c;
        z(isInverted) = -z(isInverted);   % flip so higher z = worse on every axis
        r  = 0.5 + z / 4;                 % map z=[-2,+2] -> r=[0,1]
        plotData(g,:) = max(0, min(1, r));
        nPerGrp(g) = sum(idx);
    end
    colors = COL.groups(1:4,:);

    angles = linspace(pi/2, pi/2 - 2*pi, nAxes+1);
    angles = angles(1:end-1);
    rMax = 1; nRings = 4;
    rTicks = linspace(0, rMax, nRings+1); rTicks(1) = [];

    hold(ax, 'on');
    t = linspace(0, 2*pi, 200);
    fill(ax, rMax*cos(t), rMax*sin(t), [0.97 0.97 0.98], ...
        'EdgeColor', COL.grid, 'LineWidth', 0.5);
    for ri = 1:numel(rTicks)
        % Mid-ring (cohort mean, r=0.5) is drawn slightly stronger so the
        % reader can see the reference at a glance.
        if abs(rTicks(ri) - 0.5) < 1e-6
            plot(ax, rTicks(ri)*cos(t), rTicks(ri)*sin(t), '-', ...
                'Color', [0.55 0.55 0.55], 'LineWidth', 0.7);
        else
            plot(ax, rTicks(ri)*cos(t), rTicks(ri)*sin(t), '-', ...
                'Color', COL.grid, 'LineWidth', 0.4);
        end
    end
    for ai = 1:nAxes
        plot(ax, [0 rMax*cos(angles(ai))], [0 rMax*sin(angles(ai))], '-', ...
            'Color', COL.grid, 'LineWidth', 0.4);
    end
    for ri = 1:numel(rTicks)
        text(ax, 0.02, rTicks(ri), sprintf('%.2f', rTicks(ri)), ...
            'FontSize', 5.5, 'Color', [0.55 0.55 0.55], ...
            'HorizontalAlignment','left', 'FontWeight','normal');
    end

    % polygons (fill first, lines on top)
    ph = gobjects(4,1);
    for g = 1:4
        v = plotData(g,:);  if any(isnan(v)), continue; end
        vp = [v, v(1)]; ap = [angles, angles(1)];
        fill(ax, vp.*cos(ap), vp.*sin(ap), colors(g,:), ...
            'FaceAlpha', 0.10, 'EdgeColor','none');
    end
    for g = 1:4
        v = plotData(g,:); if any(isnan(v)), continue; end
        vp = [v, v(1)]; ap = [angles, angles(1)];
        xp = vp.*cos(ap);  yp = vp.*sin(ap);
        ph(g) = plot(ax, xp, yp, '-', 'Color', colors(g,:), 'LineWidth', 1.7);
        plot(ax, xp(1:end-1), yp(1:end-1), 'o', 'MarkerSize', 4.5, ...
            'MarkerFaceColor', colors(g,:), 'MarkerEdgeColor','w','LineWidth',0.5);
    end

    lOff = rMax * 1.20;
    for ai = 1:nAxes
        xL = lOff * cos(angles(ai));  yL = lOff * sin(angles(ai));
        if abs(cos(angles(ai))) < 0.1, ha = 'center';
        elseif cos(angles(ai)) > 0,    ha = 'left';
        else,                          ha = 'right';
        end
        if abs(sin(angles(ai))) < 0.1, va = 'middle';
        elseif sin(angles(ai)) > 0,    va = 'bottom';
        else,                          va = 'top';
        end
        text(ax, xL, yL, axLabels{ai}, 'HorizontalAlignment', ha, ...
            'VerticalAlignment', va, 'FontSize', 8, 'FontWeight','bold', ...
            'Color', COL.txt);
    end
    axis(ax, 'equal', 'off');
    pad = rMax * 0.40;
    xlim(ax, [-rMax-pad, rMax+pad]);
    ylim(ax, [-rMax-pad*0.7, rMax+pad*0.9]);

    % --- ORIGINAL title + subtitle + panel letter ---------------------
    text(ax, 0, 1.35, 'Clinical Profiles', ...
        'FontSize', 8.5, 'FontWeight', 'bold', 'HorizontalAlignment', 'center', ...
        'Color', COL.accent);
    text(ax, 0, 1.25, '(0 = Best, 1 = Worst)', ...
        'FontSize', 5.5, 'HorizontalAlignment', 'center', 'Color', [0.5 0.5 0.5]);
    text(ax, 0.02, 0.98, '(b)', ...
        'Units', 'normalized', 'FontSize', 10, 'FontWeight', 'bold', ...
        'Color', COL.accent, 'HorizontalAlignment', 'left', ...
        'VerticalAlignment', 'top');

    % --- ORIGINAL legend (southoutside, horizontal, 2 cols, fontsize 8) ----
    validIdx = find(arrayfun(@(h) isgraphics(h), ph));
    if ~isempty(validIdx)
        legNames = arrayfun(@(g) sprintf('%s (n=%d)', groupLabels{g}, nPerGrp(g)), ...
            validIdx', 'UniformOutput', false);
        lg = legend(ax, ph(validIdx), legNames, ...
            'Location', 'southoutside', 'Orientation', 'horizontal', ...
            'FontSize', 8, 'Box', 'off', 'NumColumns', min(numel(validIdx), 2));
        lg.ItemTokenSize = [15 15];
    end
end


function out = computeSankeyData(T, isPD, opt)
    % Compute the per-patient × per-timepoint group matrix used by both
    % the transition heatmap (panel b) and the Sankey diagram (panel d).
    % Cohort filter ('all' or 'incident') is applied here.
    Tpd   = T(isPD & ~isnan(T.group4), :);
    yrCol = safeCol(Tpd, 'YEAR');
    sankeyYears = opt.sankeyYears;
    nT = numel(sankeyYears);

    if strcmp(opt.cohort,'incident')
        durName = '';
        for cand = {'duration_yrs','disease_duration','PD_duration'}
            if ismember(cand{1}, Tpd.Properties.VariableNames)
                durName = cand{1}; break;
            end
        end
        if ~isempty(durName)
            blMask = (yrCol == 0);
            blPID  = Tpd.PATNO(blMask);
            blDur  = safeCol(Tpd(blMask,:), durName);
            incidentPID = blPID(blDur < opt.incidentDurationYrs);
            sankeyPatIDs = intersect(unique(Tpd.PATNO), incidentPID);
            fprintf('    Sankey filter: %d incident patients (<%g yr duration at baseline).\n', ...
                numel(sankeyPatIDs), opt.incidentDurationYrs);
        else
            warning('No disease-duration variable found; using all PD subjects.');
            sankeyPatIDs = unique(Tpd.PATNO);
        end
    else
        sankeyPatIDs = unique(Tpd.PATNO);
    end
    nSP = numel(sankeyPatIDs);

    patGroupMat = nan(nSP, nT);
    for ti = 1:nT
        for pi = 1:nSP
            idx = (Tpd.PATNO == sankeyPatIDs(pi)) & (yrCol == sankeyYears(ti));
            if any(idx)
                grps = Tpd.group4(idx);
                grps = grps(~isnan(grps));
                if ~isempty(grps), patGroupMat(pi, ti) = grps(1); end
            end
        end
    end
    complete         = all(~isnan(patGroupMat), 2);
    patGroupComplete = patGroupMat(complete,:);
    nComplete        = sum(complete);

    out = struct( ...
        'sankeyYears',     sankeyYears, ...
        'sankeyPatIDs',    sankeyPatIDs, ...
        'patGroupMat',     patGroupMat, ...
        'patGroupComplete',patGroupComplete, ...
        'complete',        complete, ...
        'nComplete',       nComplete);
end


function drawTransitionHeatmap(ax, patGroupComplete, sankeyYears, COL, groupLabelsShort)
    % Row-normalised transition probabilities P(group_{t+1} | group_t).
    %   One 4x4 matrix per interval, concatenated horizontally with a
    %   thin divider. Diagonals (stayers) outlined and bolded.
    %   Sequential gray-slate colormap (matches Fig 4 panel C).
    nG    = 4;
    nT    = numel(sankeyYears);
    nIntv = nT - 1;

    if size(patGroupComplete,1) < 1 || nIntv < 1
        text(ax, 0.5, 0.5, 'No transitions to display', ...
            'Units','normalized','HorizontalAlignment','center', ...
            'FontSize', 8, 'Color', [0.5 0.5 0.5]);
        axis(ax, 'off');
        return
    end

    % --- Build [nG x (nIntv*nG)] row-normalised probability matrix ----
    H = nan(nG, nIntv * nG);
    for ii = 1:nIntv
        cm = zeros(nG, nG);
        for g1 = 1:nG, for g2 = 1:nG
            cm(g1, g2) = sum(patGroupComplete(:, ii)   == g1 & ...
                             patGroupComplete(:, ii+1) == g2);
        end, end
        rs = sum(cm, 2); rs(rs == 0) = NaN;     % avoid 0/0
        prob = cm ./ rs;                         % row-normalised
        H(:, (ii-1)*nG + (1:nG)) = prob;
        % exact counts behind the printed cells, to the log only
        fprintf('    [Fig1c] Year %d->%d, n = %d, overall stability %d/%d = %.1f%%\n', sankeyYears(ii), sankeyYears(ii+1), ...
            sum(cm(:)), sum(diag(cm)), sum(cm(:)), 100*sum(diag(cm))/sum(cm(:)));
        for g1 = 1:nG
            fprintf('    [Fig1c]    from %-8s counts %s | row %% %s\n', groupLabelsShort{g1}, sprintf('%d ', cm(g1,:)), sprintf('%.1f ', 100*prob(g1,:)));
        end
    end

    % --- Render heatmap (clim [0,1]) ----------------------------------
    imagesc(ax, H, [0 1]);

    % Sequential gray-slate colormap (like Fig 4 §c)
    cmap = interp1(linspace(0, 1, 5), ...
        [0.96 0.96 0.97; 0.78 0.82 0.88; 0.55 0.62 0.74; 0.32 0.42 0.60; 0.14 0.22 0.42], ...
        linspace(0, 1, 256), 'pchip');
    colormap(ax, cmap);

    cb = colorbar(ax, 'FontSize', f1pt(6.5));
    cb.Label.String   = 'P(state_{t+1} | state_t)';
    cb.Label.FontSize = f1pt(7.5);   % subscripts (0.8x) print at 6 pt
    set(cb, 'FontWeight','bold');
    set(cb.Label, 'FontWeight','bold');
    cb.Ticks = 0:0.25:1;

    hold(ax, 'on');

    % --- Cell text + diagonal outlines --------------------------------
    for r = 1:nG
        for c = 1:size(H, 2)
            v = H(r, c);
            if isnan(v), txt = ''; v = 0;
            else,        txt = sprintf('%.2f', v);
            end
            tc = ternary(v > 0.55, 'w', 'k');
            withinC = c - floor((c-1)/nG)*nG;
            isDiag  = (r == withinC);
            fw = ternary(isDiag, 'bold', 'normal');
            if ~isempty(txt)
                text(ax, c, r, txt, 'HorizontalAlignment','center', ...
                    'VerticalAlignment','middle', 'FontSize', f1pt(6.5), ...
                    'FontWeight', fw, 'Color', tc);
            end
            if isDiag
                rectangle(ax, 'Position', [c-0.5, r-0.5, 1, 1], ...
                    'EdgeColor', [0.10 0.10 0.10], 'LineWidth', 1.4);
            end
        end
    end

    % --- Vertical dividers between intervals --------------------------
    for d = 1:(nIntv-1)
        plot(ax, [d*nG+0.5, d*nG+0.5], [0.5, nG+0.5], '-', ...
            'Color', [0.20 0.20 0.20], 'LineWidth', 1.6);
    end

    % --- Tick labels (group abbrev repeated per interval block) -------
    %   YLim is widened on the TOP side (negative y) so the year header
    %   can sit cleanly above the heatmap with breathing room. The X-tick
    %   labels are explicitly pinned to the bottom so they cannot drift
    %   up toward the year header.
    xLabels = repmat(groupLabelsShort(1:nG), 1, nIntv);
    headerPad = 0.85;       % data-units of vertical room reserved above row 1
    % labels rotated 45 degrees (at 0 they ran together
    % as "NeitherApathy" at the interval junction); print-size fonts.
    set(ax, 'XTick', 1:size(H,2), 'XTickLabel', xLabels, ...
        'YTick', 1:nG, 'YTickLabel', groupLabelsShort(1:nG), ...
        'FontSize', f1pt(6.5), 'FontWeight','bold', 'Box','off', 'TickDir','out', ...
        'XLim', [0.5, size(H,2)+0.5], 'YLim', [0.5 - headerPad, nG+0.5], ...
        'XAxisLocation','bottom', 'XTickLabelRotation', 45);

    % --- Interval headers above each block ----------------------------
    %   Placed near the top of the (widened) YLim, well clear of row 1.
    headerY = 0.5 - headerPad * 0.55;   % ~midway between top edge and row 1
    for ii = 1:nIntv
        cx = (ii-1)*nG + (nG+1)/2;
        text(ax, cx, headerY, ...
            sprintf('Year %d \\rightarrow Year %d', sankeyYears(ii), sankeyYears(ii+1)), ...
            'HorizontalAlignment','center', 'VerticalAlignment','middle', ...
            'FontSize', f1pt(7), 'FontWeight','bold', 'Color', COL.txt);
    end
    for r = 1:nG   % printed cell values to the log
        fprintf('    [Fig1c] from %-8s %s\n', groupLabelsShort{r}, sprintf('%.2f ', H(r, :)));
    end

    ylabel(ax, 'From state (t)', 'FontSize', f1pt(7), 'FontWeight', 'bold');
    xlabel(ax, 'To state (t+1)', 'FontSize', f1pt(7), 'FontWeight', 'bold');

    % Push the panel title slightly higher so it cannot crowd the year
    % header from above (the title of an axis sits at the YLim top edge
    % by default; we offset it via TitleHorizontalAlignment + nudging the
    % Title position after it is created by the caller).
    if ~isempty(ax.Title)
        ax.TitleFontSizeMultiplier = 1.0;
        ax.TitleHorizontalAlignment = 'center';
    end
end


function drawSankeyPanelFromData(ax, sankeyData, COL, groupLabelsShort, opt)
    % Sankey transition diagram (panel d), styled after the user's
    % plotApathySankey reference:
    %   - ABSOLUTE counts on Y (not normalised to 100%)
    %   - Each timepoint column is centred VERTICALLY against maxTotal,
    %     so attrition is read off directly as shrinking column height.
    %   - SOURCE-coloured flows (alpha 0.35); linear in X, smoothstep in Y.
    %   - Group labels to the LEFT of the first column.
    %   - Year labels at the TOP of each column.
    %   - Counts INSIDE each rectangle (luminance-aware contrast).
    %   - axis OFF.
    nG = 4;
    sankeyYears      = sankeyData.sankeyYears;
    patGroupMat      = sankeyData.patGroupMat;
    patGroupComplete = sankeyData.patGroupComplete;
    nComplete        = sankeyData.nComplete;
    nT               = numel(sankeyYears);

    % --- Per-timepoint counts (use ALL patients with data at that point) ---
    counts = zeros(nG, nT);
    for k = 1:nT
        for g = 1:nG
            counts(g, k) = sum(patGroupMat(:, k) == g);
        end
    end
    maxTotal = max(sum(counts, 1));

    if maxTotal < opt.minComplete
        text(ax, 0.5, 0.5, sprintf('Sankey not drawn: max n = %d < %d', maxTotal, opt.minComplete), ...
            'Units','normalized','HorizontalAlignment','center','FontSize',8, ...
            'Color',[0.5 0.5 0.5]);
        axis(ax, 'off');
        return
    end

    % --- Pair-wise flows (consecutive timepoints, paired data only) ---
    flows = zeros(nG, nG, nT-1);
    for k = 1:(nT-1)
        pairOK = ~isnan(patGroupMat(:,k)) & ~isnan(patGroupMat(:,k+1));
        for g1 = 1:nG
            for g2 = 1:nG
                flows(g1, g2, k) = sum(pairOK & ...
                    patGroupMat(:,k)==g1 & patGroupMat(:,k+1)==g2);
            end
        end
    end

    % --- Geometry: centred-stacking per timepoint ---------------------
    nodeWidth = 0.16;   % was 0.12, so the print-size counts fit inside the nodes
    x = 1:nT;
    nodeStart = zeros(nG, nT);
    nodeEnd   = zeros(nG, nT);
    for k = 1:nT
        cum = (maxTotal - sum(counts(:, k))) / 2;
        for g = 1:nG
            nodeStart(g, k) = cum;
            nodeEnd(g, k)   = cum + counts(g, k);
            cum = nodeEnd(g, k);
        end
    end

    hold(ax, 'on');

    % --- Flow bands (linear X, smoothstep Y, source-coloured) ---------
    for k = 1:(nT-1)
        targetCursor = zeros(nG, 1);
        for g1 = 1:nG
            sCursor = nodeStart(g1, k);
            for g2 = 1:nG
                f = flows(g1, g2, k);
                if f <= 0, continue; end
                s0 = sCursor;             s1 = sCursor + f;
                t0 = nodeStart(g2, k+1) + targetCursor(g2);
                t1 = t0 + f;
                drawSankeyFlowBand(ax, ...
                    x(k) + nodeWidth/2, x(k+1) - nodeWidth/2, ...
                    s0, s1, t0, t1, COL.groups(g1,:), 0.35);
                sCursor = s1;
                targetCursor(g2) = targetCursor(g2) + f;
            end
        end
    end

    % --- Bars + count labels inside each rectangle --------------------
    for k = 1:nT
        for g = 1:nG
            h = counts(g, k);
            if h <= 0, continue; end
            col = COL.groups(g, :);
            rectangle(ax, 'Position', ...
                [x(k) - nodeWidth/2, nodeStart(g,k), nodeWidth, h], ...
                'FaceColor', col, 'EdgeColor', 'w', 'LineWidth', 0.8);
            lum = 0.2126*col(1) + 0.7152*col(2) + 0.0722*col(3);
            tc  = ternary(lum > 0.55, 'k', 'w');
            if h >= maxTotal * 0.025
                text(ax, x(k), nodeStart(g,k) + h/2, sprintf('%d', h), ...
                    'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle', ...
                    'FontSize', f1pt(6.5), 'FontWeight', 'bold', 'Color', tc);
            end
        end
    end

    % --- Year labels above each column --------------------------------
    % each column states its own n (all participants
    % classified at that visit), since the columns are not the n = 550
    % completers of panel c.
    for k = 1:nT
        text(ax, x(k), maxTotal + maxTotal*0.04, ...
            sprintf('\\bfYear %d\\rm (n = %s)', sankeyYears(k), fmtThousands(sum(counts(:, k)))), ...
            'HorizontalAlignment','center', 'VerticalAlignment','bottom', ...
            'FontSize', f1pt(7), 'FontWeight','normal', 'Color', COL.txt);
    end
    % --- Paired n of each flow, below the diagram ---
    for k = 1:(nT-1)
        text(ax, (x(k) + x(k+1))/2, -maxTotal*0.02, ...
            sprintf('paired n = %s', fmtThousands(sum(sum(flows(:, :, k))))), ...
            'HorizontalAlignment','center', 'VerticalAlignment','top', ...
            'FontSize', f1pt(6.5), 'FontAngle','italic', 'Color', [0.30 0.30 0.30]);
    end
    fprintf('    [Fig1d] nodes (all classified per visit): %s\n', strjoin(arrayfun(@(k) sprintf('Year %d = %d (%s)', ...
        sankeyYears(k), sum(counts(:,k)), strjoin(arrayfun(@(g) sprintf('%d', counts(g,k)), 1:nG, 'UniformOutput', false), '/')), ...
        1:nT, 'UniformOutput', false), '; '));
    fprintf('    [Fig1d] flows (classified at both visits): %s\n', strjoin(arrayfun(@(k) sprintf('Year %d->%d = %d', ...
        sankeyYears(k), sankeyYears(k+1), sum(sum(flows(:,:,k)))), 1:(nT-1), 'UniformOutput', false), '; '));

    % --- Group labels to the LEFT of the first column -----------------
    % leftPad mixed units (maxTotal*0.005 is in count units, so
    % the labels were drawn about 7 x-units left of the axes, off-canvas).
    leftPad = nodeWidth/2 + 0.03;
    for g = 1:nG
        if counts(g, 1) > 0
            text(ax, x(1) - leftPad, ...
                (nodeStart(g,1) + nodeEnd(g,1))/2, ...
                groupLabelsShort{g}, ...
                'HorizontalAlignment','right', 'VerticalAlignment','middle', ...
                'FontSize', f1pt(7), 'FontWeight','bold', 'Color', COL.groups(g,:));
        end
    end

    % --- End-point stability annotation -------------------------------
    % same number, clearer wording (it is computed on the
    % completers only), no subscript, placed under the last column.
    if nComplete > 0
        transMat = zeros(nG);
        for g1 = 1:nG, for g2 = 1:nG
            transMat(g1,g2) = sum(patGroupComplete(:,1)==g1 & patGroupComplete(:,end)==g2);
        end, end
        stab = 100 * sum(diag(transMat)) / nComplete;
        nWords = {'one','two','three','four','five','six'};
        text(ax, nT + 0.45, maxTotal*0.02, ...
            {sprintf('Year %d \\rightarrow %d stability: %.1f%%', sankeyYears(1), sankeyYears(end), stab), ...
             sprintf('(n = %d at all %s visits)', nComplete, nWords{min(nT, numel(nWords))})}, ...
            'HorizontalAlignment','right','VerticalAlignment','bottom', ...
            'FontSize', f1pt(6.5), 'FontWeight','normal','Color', COL.txt);
        fprintf('    [Fig1d] end-point stability Year %d -> %d: %d / %d = %.1f%%\n', sankeyYears(1), sankeyYears(end), ...
            sum(diag(transMat)), nComplete, stab);
    end

    xlim(ax, [0.5, nT + 0.5]);
    ylim(ax, [-maxTotal*0.12, maxTotal * 1.12]);   % was -0.04 (room for the paired-n labels)
    axis(ax, 'off');
end


function drawSankeyFlowBand(ax, x0, x1, s0, s1, t0, t1, clr, alphaVal)
    % Cubic-smoothstep Sankey flow band:
    %   linear interpolation in X (preserves the time axis)
    %   smoothstep w(t) = t^2 (3 - 2 t) in Y
    if s1 <= s0 || t1 <= t0, return; end
    n  = 60;
    tt = linspace(0, 1, n);
    w  = tt .* tt .* (3 - 2*tt);     % cubic smoothstep
    x  = x0 + (x1 - x0) * tt;        % LINEAR X
    y0 = s0 + (t0 - s0) * w;         % smooth Y bottom
    y1 = s1 + (t1 - s1) * w;         % smooth Y top
    patch(ax, [x, fliplr(x)], [y0, fliplr(y1)], clr, ...
        'FaceAlpha', alphaVal, 'EdgeColor', 'none');
end

% =========================================================================
%  FIGURE 1 — DIMENSIONAL COMPANION FIGURE (continuous helpers)
%  Honest, threshold-free visualisation of apathy and ICB scores.
%  Layout (2x2):
%    (a) drawJointDensityHexbin  — integer-aware count map at baseline
%                                  + faint dashed 4-group reference cutoffs
%    (b) drawClinicalSurface     — Gaussian-weighted local mean of a
%                                  clinical variable over the (apathy, ICB) plane
%    (c) drawTimepointContours   — KDE contours at year 0 and year 4 overlaid
%    (d) drawVelocityField       — mean (deltaApathy, deltaICB) per starting bin
% =========================================================================
function fig = buildFigure1Continuous(bl, T, isPD, COL)
    apV  = safeCol(bl, 'apathy_composite');
    icbV = safeCol(bl, 'quip');

    % --- Infer the 4-group cutoffs from the existing labels -----------
    [apCut, icbCut] = inferGroupCutoffs(bl);

    fig = figure('Units','centimeters','Position',[1 1 26 18],'Color',COL.bg);
    set(fig, ...
        'DefaultAxesFontName','Helvetica','DefaultTextFontName','Helvetica', ...
        'DefaultAxesFontWeight','bold','DefaultTextFontWeight','normal', ...
        'DefaultAxesTitleFontWeight','bold');
    tl = tiledlayout(fig, 2, 2, 'TileSpacing','compact', 'Padding','compact');

    % ---- Panel (a) Joint density (baseline) -------------------------
    ax_a = nexttile(tl, 1); set(ax_a,'Color',COL.bg); hold(ax_a,'on');
    drawJointDensityHexbin(ax_a, apV, icbV, apCut, icbCut, COL);
    title(ax_a, '\bfa\rm  Baseline distribution (apathy \times ICB)', ...
        'FontSize', 9, 'Color', COL.txt);

    % ---- Panel (b) Smoothed clinical surface (depression) ----------
    ax_b = nexttile(tl, 2); set(ax_b,'Color',COL.bg); hold(ax_b,'on');
    depV = safeCol(bl, 'depression_factor');
    drawClinicalSurface(ax_b, apV, icbV, depV, ...
        'Depression factor', apCut, icbCut, COL);
    title(ax_b, '\bfb\rm  Depression burden over the apathy \times ICB plane', ...
        'FontSize', 9, 'Color', COL.txt);

    % ---- Panel (c) Year 0 vs Year 4 KDE contours -------------------
    ax_c = nexttile(tl, 3); set(ax_c,'Color',COL.bg); hold(ax_c,'on');
    drawTimepointContours(ax_c, T, isPD, [0, 4], apCut, icbCut, COL);
    title(ax_c, '\bfc\rm  Cohort drift: Year 0 vs Year 4', ...
        'FontSize', 9, 'Color', COL.txt);

    % ---- Panel (d) Mean velocity field -----------------------------
    ax_d = nexttile(tl, 4); set(ax_d,'Color',COL.bg); hold(ax_d,'on');
    drawVelocityField(ax_d, T, isPD, [0, 4], apCut, icbCut, COL);
    title(ax_d, '\bfd\rm  Mean velocity field (Year 0 \rightarrow Year 4)', ...
        'FontSize', 9, 'Color', COL.txt);

    drawnow;
end


function [apCut, icbCut] = inferGroupCutoffs(bl)
    % Infer apathy & ICB cutoffs from the existing 4-group labels:
    %   group 1 = Neither (low apathy AND low ICB)
    %   group 2 = Apathy-only (high apathy, low ICB)
    %   group 3 = ICB-only (low apathy, high ICB)
    %   group 4 = Both (high on both)
    apCut = NaN; icbCut = NaN;
    if ~ismember('group4', bl.Properties.VariableNames), return; end
    apV  = safeCol(bl, 'apathy_composite');
    icbV = safeCol(bl, 'quip');
    g    = bl.group4;
    apLo = apV(g == 1 | g == 3);
    apHi = apV(g == 2 | g == 4);
    if any(~isnan(apLo)) && any(~isnan(apHi))
        apCut = (max(apLo, [], 'omitnan') + min(apHi, [], 'omitnan')) / 2;
    end
    icbLo = icbV(g == 1 | g == 2);
    icbHi = icbV(g == 3 | g == 4);
    if any(~isnan(icbLo)) && any(~isnan(icbHi))
        icbCut = (max(icbLo, [], 'omitnan') + min(icbHi, [], 'omitnan')) / 2;
    end
end


function drawJointDensityHexbin(ax, x, y, apCut, icbCut, COL)
    % Integer-aware 2D count map of apathy x ICB.
    %   One cell per (apathy, ICB) integer pair; colour = subject count.
    %   Dashed reference lines mark the conventional 4-group cutoffs.
    ok = ~isnan(x) & ~isnan(y);
    x = x(ok); y = y(ok); n = numel(x);
    if n < 5
        text(ax, 0.5, 0.5, 'No data', 'Units','normalized', ...
            'HorizontalAlignment','center');
        return
    end
    xMin = floor(min(x)); xMax = ceil(max(x));
    yMin = floor(min(y)); yMax = ceil(max(y));
    xEdges = (xMin - 0.5):1:(xMax + 0.5);
    yEdges = (yMin - 0.5):1:(yMax + 0.5);
    N = histcounts2(x, y, xEdges, yEdges);

    xCenters = (xEdges(1:end-1) + xEdges(2:end)) / 2;
    yCenters = (yEdges(1:end-1) + yEdges(2:end)) / 2;
    Nplot = N'; Nplot(Nplot == 0) = NaN;
    imagesc(ax, xCenters, yCenters, Nplot, 'AlphaData', ~isnan(Nplot));
    set(ax, 'YDir', 'normal');

    cmap = interp1(linspace(0,1,5), ...
        [0.92 0.94 0.97; 0.65 0.74 0.85; 0.40 0.52 0.70; 0.22 0.35 0.55; 0.10 0.18 0.40], ...
        linspace(0,1,256), 'pchip');
    colormap(ax, cmap);

    % cell-count labels for the largest cells
    maxN = max(N(:));
    for ix = 1:numel(xCenters)
        for iy = 1:numel(yCenters)
            v = N(ix, iy);
            if v >= max(maxN * 0.04, 5)
                tc = 'k';
                if v > maxN * 0.55, tc = 'w'; end
                text(ax, xCenters(ix), yCenters(iy), num2str(v), ...
                    'HorizontalAlignment','center','VerticalAlignment','middle', ...
                    'FontSize', 6.5, 'FontWeight','bold','Color', tc);
            end
        end
    end

    if ~isnan(apCut)
        xline(ax, apCut, '--', 'Color', [0.45 0.45 0.45], 'LineWidth', 0.9);
    end
    if ~isnan(icbCut)
        yline(ax, icbCut, '--', 'Color', [0.45 0.45 0.45], 'LineWidth', 0.9);
    end

    text(ax, 0.02, 0.98, sprintf('n = %d (baseline)', n), ...
        'Units','normalized', 'VerticalAlignment','top', ...
        'FontSize', 7, 'Color', COL.txt, ...
        'BackgroundColor',[1 1 1 0.75], 'Margin', 1);
    text(ax, 0.98, 0.02, ...
        'dashed: 4-group reference cutoffs', ...
        'Units','normalized', 'HorizontalAlignment','right', 'VerticalAlignment','bottom', ...
        'FontSize', 6, 'Color', [0.4 0.4 0.4], 'FontAngle', 'italic', ...
        'BackgroundColor',[1 1 1 0.75], 'Margin', 1);

    cb = colorbar(ax, 'FontSize', 7);
    cb.Label.String = 'Subjects per cell';
    cb.Label.FontSize = 7;
    set(cb, 'FontWeight','bold');

    xlim(ax, [xMin-0.5, xMax+0.5]);
    ylim(ax, [yMin-0.5, yMax+0.5]);
    xlabel(ax, 'Apathy composite', 'FontSize', 8, 'FontWeight','bold');
    ylabel(ax, 'QUIP-Short total (ICB)', 'FontSize', 8, 'FontWeight','bold');
    set(ax, 'FontSize', 8, 'FontWeight','bold', 'Box','off', 'TickDir','out');
end


function drawClinicalSurface(ax, x, y, z, varLabel, apCut, icbCut, COL) %#ok<INUSD>
    % LOESS-like local Gaussian-weighted mean of z over the (x, y) plane.
    %   Bandwidths default to ~0.5 SD per axis. Areas with very low local
    %   data support are masked (NaN, transparent in imagesc).
    ok = ~isnan(x) & ~isnan(y) & ~isnan(z);
    x = x(ok); y = y(ok); z = z(ok);
    if numel(x) < 30
        text(ax, 0.5, 0.5, 'Insufficient data', 'Units','normalized', ...
            'HorizontalAlignment','center');
        return
    end
    xMin = floor(min(x)); xMax = ceil(max(x));
    yMin = floor(min(y)); yMax = ceil(max(y));
    nGrid = 60;
    [Xg, Yg] = meshgrid(linspace(xMin, xMax, nGrid), ...
                         linspace(yMin, yMax, nGrid));
    bwx = max(std(x)*0.5, 0.5);
    bwy = max(std(y)*0.5, 0.5);
    Z = nan(size(Xg)); W = nan(size(Xg));
    for i = 1:numel(Xg)
        d2 = ((x - Xg(i)) / bwx).^2 + ((y - Yg(i)) / bwy).^2;
        w  = exp(-0.5 * d2);
        sw = sum(w);
        W(i) = sw;
        if sw > 1.0
            Z(i) = sum(w .* z) / sw;
        end
    end
    Z(W < 2) = NaN;   % require >=2 effective subjects (sum of kernel weights)

    imagesc(ax, [xMin xMax], [yMin yMax], Z, 'AlphaData', ~isnan(Z));
    set(ax, 'YDir', 'normal');

    % Sequential colormap scaled to the data range (single-hue warm; "burden")
    z_min = min(Z(:),[],'omitnan');
    z_max = max(Z(:),[],'omitnan');
    if isnan(z_min) || isnan(z_max) || z_min >= z_max
        z_min = mean(z) - max(std(z), eps);
        z_max = mean(z) + max(std(z), eps);
    end
    set(ax, 'CLim', [z_min, z_max]);
    % White -> red (the red half of the previous diverging map)
    cmap_seq = interp1([0, 0.5, 1], ...
        [1 1 1; 0.93 0.65 0.65; 0.85 0.30 0.30], ...
        linspace(0, 1, 256), 'pchip');
    colormap(ax, cmap_seq);

    if ~isnan(apCut)
        xline(ax, apCut, '--', 'Color', [0.45 0.45 0.45], 'LineWidth', 0.9);
    end
    if ~isnan(icbCut)
        yline(ax, icbCut, '--', 'Color', [0.45 0.45 0.45], 'LineWidth', 0.9);
    end

    cb = colorbar(ax, 'FontSize', 7);
    cb.Label.String = varLabel;
    cb.Label.FontSize = 7;
    set(cb, 'FontWeight','bold');

    xlim(ax, [xMin xMax]); ylim(ax, [yMin yMax]);
    xlabel(ax, 'Apathy composite', 'FontSize', 8, 'FontWeight','bold');
    ylabel(ax, 'QUIP-Short total (ICB)', 'FontSize', 8, 'FontWeight','bold');
    set(ax, 'FontSize', 8, 'FontWeight','bold', 'Box','off', 'TickDir','out');
end


function drawTimepointContours(ax, T, isPD, years, apCut, icbCut, COL) %#ok<INUSD>
    % Bivariate KDE contours for each requested timepoint, overlaid on
    % the same axes for direct comparison of cohort drift.
    Tpd = T(isPD, :);
    yrCol = safeCol(Tpd, 'YEAR');
    apCol = safeCol(Tpd, 'apathy_composite');
    icCol = safeCol(Tpd, 'quip');
    okAll = ~isnan(apCol) & ~isnan(icCol) & ismember(yrCol, years);
    if sum(okAll) < 30
        text(ax, 0.5, 0.5, 'Insufficient longitudinal data', ...
            'Units','normalized', 'HorizontalAlignment','center');
        return
    end
    xLim  = [floor(min(apCol(okAll))), ceil(max(apCol(okAll)))];
    yLim  = [floor(min(icCol(okAll))), ceil(max(icCol(okAll)))];
    nGrid = 60;
    [Xg, Yg]  = meshgrid(linspace(xLim(1), xLim(2), nGrid), ...
                          linspace(yLim(1), yLim(2), nGrid));
    grid_pts  = [Xg(:), Yg(:)];

    cols = [0.20 0.40 0.65; 0.85 0.30 0.30];   % blue Y0, red Y4
    h        = gobjects(numel(years), 1);
    legNames = cell(numel(years), 1);
    for i = 1:numel(years)
        yr = years(i);
        ok = ~isnan(apCol) & ~isnan(icCol) & yrCol == yr;
        if sum(ok) < 20, continue; end
        x = apCol(ok); y = icCol(ok);
        bw = [max(std(x)*0.6, 0.4), max(std(y)*0.6, 0.4)];
        d  = mvksdensity([x y], grid_pts, 'Bandwidth', bw);
        F  = reshape(d, size(Xg));
        levs = max(F(:)) * [0.20, 0.45, 0.75];
        [~, hC] = contour(ax, Xg, Yg, F, levs, ...
            'Color', cols(i,:), 'LineWidth', 1.5);
        h(i) = hC;
        legNames{i} = sprintf('Year %d (n=%d)', yr, sum(ok));
    end

    if ~isnan(apCut), xline(ax, apCut, '--', 'Color', [0.45 0.45 0.45], 'LineWidth', 0.7); end
    if ~isnan(icbCut), yline(ax, icbCut, '--', 'Color', [0.45 0.45 0.45], 'LineWidth', 0.7); end

    valid = arrayfun(@(o) isgraphics(o), h);
    if any(valid)
        lg = legend(ax, h(valid), legNames(valid), ...
            'Location','northeast', 'FontSize', 7, 'Box','off');
        lg.ItemTokenSize = [12 12];
    end

    xlim(ax, xLim); ylim(ax, yLim);
    xlabel(ax, 'Apathy composite', 'FontSize', 8, 'FontWeight','bold');
    ylabel(ax, 'QUIP-Short total (ICB)', 'FontSize', 8, 'FontWeight','bold');
    set(ax, 'FontSize', 8, 'FontWeight','bold', 'Box','off', 'TickDir','out');
end


function drawVelocityField(ax, T, isPD, years, apCut, icbCut, COL)
    % Per-bin mean velocity vector between years(1) and years(end).
    %   Each subject contributes one (deltaApathy, deltaICB) vector;
    %   vectors are averaged within bins of starting position.
    %   Faint background trajectories drawn first for context.
    Tpd   = T(isPD, :);
    yrCol = safeCol(Tpd, 'YEAR');
    apCol = safeCol(Tpd, 'apathy_composite');
    icCol = safeCol(Tpd, 'quip');
    yr1 = years(1); yr2 = years(end);

    pids = unique(Tpd.PATNO);
    apY1 = nan(numel(pids), 1); apY2 = nan(numel(pids), 1);
    icY1 = nan(numel(pids), 1); icY2 = nan(numel(pids), 1);
    for pi = 1:numel(pids)
        sel1 = Tpd.PATNO == pids(pi) & yrCol == yr1;
        sel2 = Tpd.PATNO == pids(pi) & yrCol == yr2;
        if any(sel1)
            apY1(pi) = mean(apCol(sel1), 'omitnan');
            icY1(pi) = mean(icCol(sel1), 'omitnan');
        end
        if any(sel2)
            apY2(pi) = mean(apCol(sel2), 'omitnan');
            icY2(pi) = mean(icCol(sel2), 'omitnan');
        end
    end
    ok = ~isnan(apY1) & ~isnan(apY2) & ~isnan(icY1) & ~isnan(icY2);
    apY1 = apY1(ok); apY2 = apY2(ok);
    icY1 = icY1(ok); icY2 = icY2(ok);
    n = sum(ok);

    if n < 30
        text(ax, 0.5, 0.5, sprintf('Insufficient paired data (n=%d)', n), ...
            'Units','normalized', 'HorizontalAlignment','center');
        return
    end

    % --- Faint background trajectories (with tiny jitter for visibility) -
    rng(42);   % reproducible jitter
    jit = 0.06;
    for i = 1:numel(apY1)
        plot(ax, [apY1(i), apY2(i)] + jit*randn(1,2), ...
                  [icY1(i), icY2(i)] + jit*randn(1,2), '-', ...
            'Color', [0.55 0.55 0.55 0.10], 'LineWidth', 0.4);
    end

    % --- Bin starting positions ---------------------------------------
    xLim = [floor(min(apY1)), ceil(max(apY1))];
    yLim = [floor(min(icY1)), ceil(max(icY1))];
    nB   = 6;
    xEdges = linspace(xLim(1), xLim(2), nB+1);
    yEdges = linspace(yLim(1), yLim(2), nB+1);
    [Xg, Yg] = meshgrid((xEdges(1:end-1)+xEdges(2:end))/2, ...
                         (yEdges(1:end-1)+yEdges(2:end))/2);
    U = nan(size(Xg)); V = nan(size(Xg));
    for ix = 1:nB
        for iy = 1:nB
            inBin = apY1 >= xEdges(ix) & apY1 < xEdges(ix+1) & ...
                    icY1 >= yEdges(iy) & icY1 < yEdges(iy+1);
            if sum(inBin) >= 4
                U(iy, ix) = mean(apY2(inBin) - apY1(inBin));
                V(iy, ix) = mean(icY2(inBin) - icY1(inBin));
            end
        end
    end

    % --- Aggregate (cohort-mean) drift -------------------------------------
    mu_u = mean(apY2 - apY1);
    mu_v = mean(icY2 - icY1);
    x0   = mean(apY1);
    y0   = mean(icY1);

    % --- Common scale so per-bin AND aggregate are directly comparable -----
    %   Reference = MEDIAN per-bin magnitude (a typical bin fills ~70% of a
    %   bin width). Outlier bins overflow into neighbours — the correct visual
    %   signal that those bins are unusually large. The aggregate is drawn on
    %   the same scale, so its length is directly readable relative to a
    %   "typical bin".
    bin_dx  = mean(diff(xEdges));
    bin_dy  = mean(diff(yEdges));
    mags    = sqrt(U(:).^2 + V(:).^2);
    ref_mag = median(mags(~isnan(mags)));
    sc      = 0.7 * min(bin_dx, bin_dy) / max(ref_mag, eps);

    % --- Per-bin vector field ----------------------------------------------
    quiver(ax, Xg, Yg, sc*U, sc*V, 'AutoScale','off', ...
        'Color', [0.10 0.20 0.50], 'LineWidth', 1.8, 'MaxHeadSize', 1.0);

    if ~isnan(apCut), xline(ax, apCut, '--', 'Color', [0.45 0.45 0.45], 'LineWidth', 0.7); end
    if ~isnan(icbCut), yline(ax, icbCut, '--', 'Color', [0.45 0.45 0.45], 'LineWidth', 0.7); end

    % --- Aggregate arrow (same scale as per-bin) ---------------------------
    plot(ax, x0, y0, 'o', 'MarkerSize', 5, ...
        'MarkerFaceColor', [0.85 0.20 0.20], ...
        'MarkerEdgeColor', 'w', 'LineWidth', 0.8);
    quiver(ax, x0, y0, sc*mu_u, sc*mu_v, 'AutoScale','off', ...
        'Color', [0.85 0.20 0.20], 'LineWidth', 2.6, 'MaxHeadSize', 1.5);

    text(ax, 0.02, 0.98, ...
        sprintf('n = %d (year %d \\rightarrow %d)\nmean drift = (%.2f, %.2f)', ...
                n, yr1, yr2, mu_u, mu_v), ...
        'Units','normalized', 'VerticalAlignment','top', ...
        'FontSize', 7, 'Color', COL.txt, ...
        'BackgroundColor',[1 1 1 0.75], 'Margin', 1);

    xlim(ax, xLim); ylim(ax, yLim);
    xlabel(ax, 'Apathy composite (baseline)', 'FontSize', 8, 'FontWeight','bold');
    ylabel(ax, 'QUIP-Short total (baseline)', 'FontSize', 8, 'FontWeight','bold');
    set(ax, 'FontSize', 8, 'FontWeight','bold', 'Box','off', 'TickDir','out');
end

function radarPlotPub(data, axisLabels, groupNames, colors, COL)
    % Publication-quality radar/spider chart
    % data: nGroups x nAxes matrix of values
    
    nAxes = size(data, 2);
    nGroups = size(data, 1);
    
    % Angles: start from top (pi/2) going clockwise
    angles = linspace(pi/2, pi/2 - 2*pi, nAxes+1);
    angles = angles(1:end-1);
    
    ax = gca; 
    hold(ax, 'on');
    
    % Determine plotting range
    allVals = data(:); 
    allVals = allVals(~isnan(allVals));
    if isempty(allVals), return; end
    
    absMax = max(abs(allVals)) + 0.1; % Add slight padding
    rMax = ceil(absMax * 5) / 5;      % Round to neat intervals
    
    rMin = -rMax;  
    if min(allVals) >= -0.05
        rMin = 0;  % If all positive (like our 0-1 data), start from 0
    end
    
    % Grid rings
    nRings = 5;
    rTicks = linspace(rMin, rMax, nRings);
    t = linspace(0, 2*pi, 200);
    
    % Light background fill
    fill(ax, rMax*cos(t), rMax*sin(t), [0.97 0.97 0.97], 'EdgeColor', COL.grid, 'LineWidth', 0.5);
    
    for ri = 1:nRings
        if rTicks(ri) == 0
            plot(ax, rTicks(ri)*cos(t), rTicks(ri)*sin(t), '-', 'Color', [0.6 0.6 0.6], 'LineWidth', 0.6);
        else
            plot(ax, abs(rTicks(ri))*cos(t), abs(rTicks(ri))*sin(t), '-', 'Color', COL.grid, 'LineWidth', 0.3);
        end
        
        % Numeric label on first spoke
        rAbs = rTicks(ri);
        if rAbs ~= 0 || rMin < 0
            text(ax, 0.05, rAbs + 0.05, sprintf('%.1f', rAbs), ...
                'FontSize', 6, 'Color', [0.5 0.5 0.5], 'HorizontalAlignment', 'left');
        end
    end
    
    % Spoke lines
    for ai = 1:nAxes
        plot(ax, [0 rMax*cos(angles(ai))], [0 rMax*sin(angles(ai))], '-', ...
            'Color', COL.grid, 'LineWidth', 0.4);
    end
    
    % Plot data polygons: fill first
    ph = gobjects(nGroups, 1);
    for g = 1:nGroups
        vals = data(g,:);
        if any(isnan(vals)), continue; end
        
        vp = [vals, vals(1)]; 
        ap = [angles, angles(1)];
        xp = vp .* cos(ap); 
        yp = vp .* sin(ap);
        
        fill(ax, xp, yp, colors(g,:), 'FaceAlpha', 0.1, ...
            'EdgeColor', colors(g,:), 'LineWidth', 0.5, 'EdgeAlpha', 0.3);
    end
    
    % Lines and markers on top (separate loop to avoid overlap issues)
    for g = 1:nGroups
        vals = data(g,:);
        if any(isnan(vals)), continue; end
        
        vp = [vals, vals(1)]; 
        ap = [angles, angles(1)];
        xp = vp .* cos(ap); 
        yp = vp .* sin(ap);
        
        ph(g) = plot(ax, xp, yp, '-', 'Color', colors(g,:), 'LineWidth', 1.6, ...
            'MarkerFaceColor', colors(g,:), 'MarkerSize', 5, 'MarkerEdgeColor', 'w');
    end
    
    % Axis labels with positioning logic
    lOff = rMax * 1.15;
    for ai = 1:nAxes
        xL = lOff * cos(angles(ai));
        yL = lOff * sin(angles(ai));
        
        if abs(cos(angles(ai))) < 0.1
            ha = 'center';
        elseif cos(angles(ai)) > 0
            ha = 'left';
        else
            ha = 'right';
        end
        
        if abs(sin(angles(ai))) < 0.1
            va = 'middle';
        elseif sin(angles(ai)) > 0
            va = 'bottom';
        else
            va = 'top';
        end
        
        text(ax, xL, yL, axisLabels{ai}, 'HorizontalAlignment', ha, ...
            'VerticalAlignment', va, 'FontSize', 8, 'FontWeight', 'bold', 'Color', COL.txt);
    end
    
    axis(ax, 'equal', 'off');
    pad = rMax * 0.35;
    xlim(ax, [-rMax-pad, rMax+pad]);
    ylim(ax, [-rMax-pad, rMax+pad]);
    
    % Legend
    validH = ph(arrayfun(@(h) isgraphics(h) && isvalid(h), ph));
    validNames = groupNames(arrayfun(@(h) isgraphics(h) && isvalid(h), ph));
    if numel(validH) > 0
        lg = legend(ax, validH, validNames, ...
            'Location', 'southoutside', 'Orientation', 'horizontal', ...
            'FontSize', 8, 'Box', 'off', 'NumColumns', min(numel(validH), 2));
        lg.ItemTokenSize = [15 15];
    end
    
    hold(ax, 'off');
end

function drawEllipse(data, col, confidence)
    mu = mean(data); C = cov(data);
    [V, D] = eig(C);
    chi2val = chi2inv(confidence, 2);
    t = linspace(0,2*pi,100);
    a = sqrt(chi2val*D(1,1)); b = sqrt(chi2val*D(2,2));
    ell = V * [a*cos(t); b*sin(t)];
    plot(ell(1,:)+mu(1), ell(2,:)+mu(2), '-', 'Color',[col 0.4],'LineWidth',1);
end
function v = sCol(T, name)
    if ismember(name, T.Properties.VariableNames), v = safeNum(T.(name));
    else, v = nan(height(T),1);
    end
end
function v = safeNum(x)
    if isnumeric(x), v = double(x);
    elseif iscell(x)||isstring(x), v = str2double(string(x));
    elseif iscategorical(x), v = double(x);
    elseif islogical(x), v = double(x);
    else, v = nan(numel(x),1);
    end
end
function fmtAx(ax, COL)
    set(ax,'FontSize',7,'Box','off','TickDir','out','XGrid','on','YGrid','on',...
        'GridColor',COL.grid,'GridAlpha',0.4);
end

function c = f1pt(p)
%F1PT printed point size -> canvas point size for Figure 1:
%the 26 cm wide canvas is included at 170 mm in main_CLD_revised.tex.
    c = p * 26 / 17;
end

function s = fmtThousands(n)
%FMTTHOUSANDS integer with thousands separators, e.g. 1,502.
    s = sprintf('%d', round(n));
    s = regexprep(s, '(\d)(?=(\d{3})+$)', '$1,');
end

function auditFigureText(fig, scale, csvPath)
%AUDITFIGURETEXT legibility QA: lists every visible text
%element (tick labels, text objects incl. titles/labels, legends, colorbar
%ticks and labels) with its printed size = FontSize x scale, where scale =
%inclusion width / canvas width. TeX sub/superscripts are drawn at 0.8x the
%base size, reported separately. Writes a CSV and prints the minima.
    rows = cell(0, 6);
    hasScript = @(str, interp) strcmpi(interp, 'tex') && ~isempty(regexp(str, '[_^]', 'once'));
    axs = findall(fig, 'Type', 'axes');
    for k = 1:numel(axs)
        a = axs(k);
        if ~strcmp(a.Visible, 'on'), continue; end
        tagA = sprintf('axes@[%.3f %.3f]', a.Position(1), a.Position(2));
        rulers = {'x', a.XAxis, a.XColor, a.XTickLabel; 'y', a.YAxis(1), a.YColor, a.YTickLabel};
        for r = 1:2
            col = rulers{r,3}; lbl = rulers{r,4};
            if (ischar(col) && strcmp(col, 'none')) || isempty(lbl), continue; end
            lbl = cellstr(lbl); lblStr = strjoin(lbl(:)', ' | ');
            fs = rulers{r,2}.FontSize;
            sc = hasScript(lblStr, a.TickLabelInterpreter);
            rows(end+1,:) = {[rulers{r,1} 'tick'], tagA, lblStr, fs, fs*scale, ternary(sc, 0.8*fs*scale, NaN)}; %#ok<AGROW>
        end
    end
    txt = findall(fig, 'Type', 'text');
    for k = 1:numel(txt)
        t = txt(k); s = t.String;
        if iscell(s), s = strjoin(s(:)', ' / '); end
        s = char(s); if size(s,1) > 1, s = strjoin(cellstr(s)', ' / '); end
        if isempty(strtrim(s)) || ~strcmp(t.Visible, 'on'), continue; end
        par = t.Parent;
        if isa(par, 'matlab.graphics.axis.Axes') && strcmp(par.Visible, 'off') && (isequal(t, par.XLabel) || isequal(t, par.YLabel)), continue; end
        fs = t.FontSize;
        rows(end+1,:) = {'text', class(par), s, fs, fs*scale, ternary(hasScript(s, t.Interpreter), 0.8*fs*scale, NaN)}; %#ok<AGROW>
    end
    lgs = findall(fig, 'Type', 'legend');
    for k = 1:numel(lgs)
        L = lgs(k); if ~strcmp(L.Visible, 'on'), continue; end
        s = strjoin(cellstr(L.String), ' | '); fs = L.FontSize;
        rows(end+1,:) = {'legend', 'legend', s, fs, fs*scale, ternary(hasScript(s, L.Interpreter), 0.8*fs*scale, NaN)}; %#ok<AGROW>
    end
    cbs = findall(fig, 'Type', 'colorbar');
    for k = 1:numel(cbs)
        C = cbs(k); if ~strcmp(C.Visible, 'on'), continue; end
        s = strjoin(cellstr(C.TickLabels), ' | '); fs = C.FontSize;
        rows(end+1,:) = {'cbtick', 'colorbar', s, fs, fs*scale, NaN}; %#ok<AGROW>
        s = C.Label.String; if iscell(s), s = strjoin(s(:)', ' / '); end
        if ~isempty(strtrim(char(s)))
            fs = C.Label.FontSize;
            rows(end+1,:) = {'cblabel', 'colorbar', char(s), fs, fs*scale, ternary(hasScript(char(s), C.Label.Interpreter), 0.8*fs*scale, NaN)}; %#ok<AGROW>
        end
    end
    T = cell2table(rows, 'VariableNames', {'kind', 'parent', 'text', 'canvas_pt', 'printed_pt', 'printed_script_pt'});
    writetable(T, csvPath);
    fprintf('    [font audit] %d text elements; min printed %.2f pt (%s: "%s"); min printed sub/superscript %.2f pt -> %s\n', ...
        height(T), min(T.printed_pt), T.kind{find(T.printed_pt == min(T.printed_pt), 1)}, ...
        T.text{find(T.printed_pt == min(T.printed_pt), 1)}, min(T.printed_script_pt, [], 'omitnan'), csvPath);
end

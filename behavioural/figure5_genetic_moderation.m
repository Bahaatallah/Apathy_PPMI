% Role: Genetic moderation LMEs (35-term FDR), apathy/ICB trajectory tests by genotype, per-subgroup
%     correlations and the genetics figure, which is manuscript Figure 5 (the file keeps the name of
%     the original figure5_genetic_moderation.m).
% Run: on the workspace left by behavioural/cohort_phenotyping.m, loaded in a fresh MATLAB
%     session, with APATHY_MODE_OVERRIDE absent (default modes: APATHY_MODE 'continuous',
%     STRATA_MODE 'visit'); outputs go to outputDir. The baseline-apathy refit of Supplementary
%     Table S17 sets APATHY_MODE_OVERRIDE = 'baseline' after copying the baseline apathy and ICB
%     z-scores of lmeT into Tpd_lme as bl_apathy_z and bl_icb_z.
% Manuscript item: Figure 5; Results (GBA1 section); Supplementary Tables S16 and S17.
%% Generates manuscript Figure 5 (Fig3_Genetics_Motivation_Cognition; the PNG keeps the historical
%  file name). Text is sized for print at the 140 mm inclusion width. Optional: FIG_FONT_AUDIT = true
%  writes a CSV of printed font sizes.

%% =========================================================================
%  FIG 3: GENETICS, MOTIVATION, AND COGNITIVE DECLINE
%  =========================================================================
%  Narrative:
%    Row 1 — Genetic landscape
%    Row 2 — Apathy & ICB trajectories by genetics
%    Row 3 — Genetics x motivation interaction on cognitive decline
%
%  Statistical:
%    - DVs z-scored for comparable betas
%    - Full covariates: age, sex, education, UPDRS-III, LEDD, depression
%    - LME: DV ~ Year * Gene * Apathy + covariates + (Year|PATNO)
%    - FDR correction across interaction terms
%    - Trajectory visualisation uses within-subject error correction
%      (Loftus & Masson 1994, Morey 2008), with Gaussian smoothing
%
%  Dependencies: safeCol, exportFigure, formatAxis, COL, THRESH
%  =========================================================================

fprintf('\n%s\n', repmat('=',1,72));
fprintf('  FIG 3: GENETICS, MOTIVATION & COGNITIVE DECLINE\n');
fprintf('%s\n\n', repmat('=',1,72));

%% ====  CONFIG: COGNITIVE VARIABLE  ======================================
%  Single knob for the cognitive DV used throughout this script.
%  Alternatives to test:
%     'cfa_g'            - CFA global factor (correlated-factors model)
%     'bifactor_g'       - Bifactor g-factor
%     'moca'             - MoCA total (0-30)
%     'obs_g'            - Observed composite
%     'COG_COMPOSITE_INT'- Integer cognitive composite
%
%  The selected variable is z-scored and stored in Tg.cog_z as the DV
%  for all LME models and trajectory plots below.
COG_VAR       = 'bifactor_g';
COG_VAR_LABEL = sprintf('%s (z)', COG_VAR);

% Apathy parameterisation in LME models:
%   'continuous' - z-scored apathy_composite (apathy_z); recommended.
%                  Beta = effect per 1 SD of apathy. Higher power, no
%                  threshold sensitivity.
%   'binary'     - apathy_composite >= THRESH.apathy_comp (apathy_bin).
%                  Beta = effect of being above clinical threshold.
%
% NOTE: visualisation strata (panel i) always use the
% binary split for interpretability — only the LME models are affected.
if exist('APATHY_MODE_OVERRIDE', 'var'), APATHY_MODE = APATHY_MODE_OVERRIDE; else, APATHY_MODE = 'continuous'; end   % 'baseline' = Year 0 apathy and ICB z (bl_apathy_z, bl_icb_z), the primary-model design

fprintf('  Cognitive DV : %s  (z-scored -> Tg.cog_z)\n', COG_VAR);
fprintf('  Apathy mode  : %s\n\n', APATHY_MODE);

%% ====  PALETTE + LEGEND STYLE  ==========================================
%  Deterministic colour-blind-safe palette, identical across the revised
%  figure scripts. Overrides the workspace COL from the master script so
%  run order cannot reintroduce the old colours.
COL.neither = [0.00 0.60 0.50];   % Okabe-Ito bluish green
COL.apathy  = [0.90 0.60 0.00];   % Okabe-Ito orange
COL.icb     = [0.00 0.45 0.70];   % Okabe-Ito blue
COL.overlap = [0.45 0.20 0.50];   % deep purple
COL.groups  = [COL.neither; COL.apathy; COL.icb; COL.overlap];
COL.bg      = [1 1 1];
COL.grid    = [0.85 0.85 0.85];
COL.txt     = [0.15 0.15 0.15];
COL.accent  = [0.10 0.10 0.10];

%  One legend style for every panel (was 3.5 to 5 pt with mixed token sizes).
%  sizes are now PRINT sizes at the 140 mm inclusion width.
PRINT_W_CM = 14.0;                    % \includegraphics width in main_CLD_revised.tex
FIG_W_CM   = 20.0;                    % authored canvas width (unchanged)
FIG_H_CM   = 21.0;                    % authored canvas height (was 22; prints 14.7 cm so float + caption fit the 22 cm text block)
PRINT_H_CM = FIG_H_CM * PRINT_W_CM / FIG_W_CM;
PT = @(p) p * FIG_W_CM / PRINT_W_CM;  % printed pt -> canvas pt
posCM = @(l, b, w, h) [l/PRINT_W_CM, b/PRINT_H_CM, w/PRINT_W_CM, h/PRINT_H_CM];  % printed cm -> normalized
FS_TITLE  = PT(7.5);   % panel titles
FS_LABEL  = PT(7.5);   % axis labels (7.5 so beta_z / p_FDR subscripts print at 6)
FS_TICK   = PT(6.5);   % tick labels
FS_ANNOT  = PT(6.0);   % in-panel numbers (cell values, n labels)
FS_SCRIPT = PT(7.5);   % annotations that contain a sub/superscript
LEG_FS    = PT(6.5);   % legend font size, all panels
LEG_TOKEN = [14 8];    % legend ItemTokenSize (points), all panels except i
EPS_ = char(949); MINUS_ = char(8722); TIMES_ = char(215);   % Greek epsilon, minus sign, multiplication sign

%% ====  0. DATA PREPARATION  =============================================
fprintf('[0] Preparing data...\n');

assert(exist('Tpd_lme','var')==1, 'Tpd_lme required.');
Tg = Tpd_lme;
fprintf('    Longitudinal table: %d rows, %d subjects\n', height(Tg), numel(unique(Tg.PATNO)));

[~, firstIdx] = unique(Tg.PATNO, 'first');
bl_g = Tg(firstIdx, :);
fprintf('    Baseline table:     %d subjects\n', height(bl_g));

% --- Clean genetic subgroup labels ---
fprintf('    Cleaning genetic subgroup labels...\n');
sgRaw = string(Tg.subgroup);

rawUni = unique(sgRaw);
fprintf('    Raw subgroup labels (%d unique):\n', numel(rawUni));
for ri = 1:numel(rawUni)
    fprintf('      "%s"  (n=%d)\n', rawUni(ri), sum(sgRaw == rawUni(ri)));
end

sgClean = repmat("", height(Tg), 1);
sgClean(contains(sgRaw,'LRRK2','IgnoreCase',true) & contains(sgRaw,'GBA','IgnoreCase',true)) = "LRRK2+GBA";
sgClean(contains(sgRaw,'LRRK2','IgnoreCase',true) & sgClean == "") = "LRRK2";
sgClean(contains(sgRaw,'GBA','IgnoreCase',true)   & sgClean == "") = "GBA";
sgClean(contains(sgRaw,'SNCA','IgnoreCase',true)   & sgClean == "") = "SNCA";
sgClean((contains(sgRaw,'PRKN','IgnoreCase',true) | contains(sgRaw,'PINK1','IgnoreCase',true)) & sgClean == "") = "PRKN/PINK1";

isNormo = contains(lower(sgRaw), 'normosomic') | contains(lower(sgRaw), 'normosmic') ...
        | contains(lower(sgRaw), 'sporadic')   | contains(lower(sgRaw), 'idiopathic');
sgClean(isNormo & sgClean == "") = "Sporadic";
sgClean(sgClean == "" & ~ismissing(sgRaw) & strlength(sgRaw) > 0) = "Sporadic";
sgClean(sgClean == "") = "Sporadic";

Tg.gene_subgroup = sgClean;
bl_g.gene_subgroup = sgClean(firstIdx);

cleanUni = unique(sgClean);
fprintf('    Cleaned subgroup labels (%d unique):\n', numel(cleanUni));
for ri = 1:numel(cleanUni)
    fprintf('      %-15s  n_subjects = %d\n', cleanUni(ri), sum(bl_g.gene_subgroup == cleanUni(ri)));
end

% --- APOE, apathy, group4 ---
apoe4_raw = double(safeCol(Tg, 'APOE_e4'));
Tg.apoe4_bin = double(apoe4_raw >= 1); Tg.apoe4_bin(isnan(apoe4_raw)) = NaN;
bl_g.apoe4_bin = Tg.apoe4_bin(firstIdx);

ap_cont = double(safeCol(Tg, 'apathy_composite'));
Tg.apathy_bin = double(ap_cont >= THRESH.apathy_comp); Tg.apathy_bin(isnan(ap_cont)) = NaN;
bl_g.apathy_bin = Tg.apathy_bin(firstIdx);

% --- trajectory strata (panel i)
%     by BASELINE (Year 0) apathy status, broadcast to every visit row, so
%     that participants are counted once and the caption ('baseline apathy
%     status') is true. STRATA_MODE = 'visit' restores the per-visit split.
STRATA_MODE = 'visit';      % default reproduces the published figure (per-visit strata); 'baseline' available
yr0 = double(safeCol(Tg, 'YEAR'));
bl0 = Tg(yr0 == 0, :); [~, u0] = unique(bl0.PATNO, 'first');
apBL = ap_cont(yr0 == 0); apBL = apBL(u0); ids0 = string(bl0.PATNO(u0));
[tfB, locB] = ismember(string(Tg.PATNO), ids0);
apBLrow = nan(height(Tg), 1); apBLrow(tfB) = apBL(locB(tfB));
Tg.apathy_bin_bl = double(apBLrow >= THRESH.apathy_comp); Tg.apathy_bin_bl(isnan(apBLrow)) = NaN;
fprintf('    Baseline strata flag: %d participants with a Year 0 composite; first-row-not-Year-0 = %d\n', ...
    sum(~isnan(apBL)), sum(double(safeCol(bl_g, 'YEAR')) ~= 0));
if strcmpi(STRATA_MODE, 'baseline'), Tg.apathy_bin_disp = Tg.apathy_bin_bl; else, Tg.apathy_bin_disp = Tg.apathy_bin; end

if ismember('group4', Tg.Properties.VariableNames)
    bl_g.group4 = double(safeCol(bl_g, 'group4'));
    Tg.group4_num = double(safeCol(Tg, 'group4'));
    hasGroup4 = true;
else
    hasGroup4 = false;
    warning('group4 not found.');
end

Tg.year_num = double(safeCol(Tg, 'YEAR'));

% --- Z-score DVs ---
fprintf('    Z-scoring continuous variables...\n');
zscore_safe = @(x) (x - nanmean(x)) ./ nanstd(x);
fitLME_safe  = @(tbl,f1,f2) tryFitLME(tbl,f1,f2);   % must exist before first use (panel d)
fitGLME_safe = @(tbl,f1,f2,distName,linkName) tryFitGLME(tbl,f1,f2,distName,linkName);
Tg.cog_z    = zscore_safe(double(safeCol(Tg, COG_VAR)));
Tg.apathy_z = zscore_safe(double(safeCol(Tg, 'apathy_composite')));
fprintf('    DV (%s) non-missing: %d / %d rows\n', COG_VAR, sum(~isnan(Tg.cog_z)), height(Tg));

% Pick which apathy column the LME models use (set above by APATHY_MODE).
% apathy_bin is always built (the panel i strata use it for display).
switch lower(APATHY_MODE)
    case 'continuous', apathy_col_lme = 'apathy_z';
    case 'binary',     apathy_col_lme = 'apathy_bin';
    case 'baseline',   apathy_col_lme = 'bl_apathy_z';   % baseline z from the master workspace (same as the primary cognition LMEs)
    otherwise, error('Unknown APATHY_MODE "%s" (use ''continuous'', ''binary'' or ''baseline'')', APATHY_MODE);
end
fprintf('    LME apathy column: Tg.%s\n', apathy_col_lme);

hasQUIP = ismember('quip', Tg.Properties.VariableNames);
if hasQUIP, Tg.quip_z = zscore_safe(double(safeCol(Tg, 'quip'))); end

% --- ICB columns (parallel to apathy) ---
%   icb_z   : continuous z-scored measure (preferred: QUIP; fallback: group4)
%   icb_bin : 1 if any ICB symptom (QUIP > 0, or group4 in {3,4})
% Used for the gene x icb x time analysis added below.
if hasQUIP
    icb_cont   = double(safeCol(Tg, 'quip'));
    Tg.icb_z   = Tg.quip_z;
    if isfield(THRESH,'icb_quip'), Tg.icb_bin = double(icb_cont >= THRESH.icb_quip);
    else,                          Tg.icb_bin = double(icb_cont > 0); end
    Tg.icb_bin(isnan(icb_cont)) = NaN;
elseif hasGroup4
    warning('No QUIP found; deriving ICB from group4 only (no continuous version).');
    Tg.icb_bin = double(Tg.group4_num == 3 | Tg.group4_num == 4);
    Tg.icb_bin(isnan(Tg.group4_num)) = NaN;
    Tg.icb_z   = Tg.icb_bin;
else
    warning('No ICB measure available (no quip, no group4) — ICB models will be skipped.');
    Tg.icb_z   = nan(height(Tg),1);
    Tg.icb_bin = nan(height(Tg),1);
end
bl_g.icb_bin = Tg.icb_bin(firstIdx);

% Pick LME ICB column following the same APATHY_MODE rule (parallel parameterisation)
switch lower(APATHY_MODE)
    case 'continuous', icb_col_lme = 'icb_z';
    case 'binary',     icb_col_lme = 'icb_bin';
end
if strcmpi(APATHY_MODE, 'baseline'), icb_col_lme = 'bl_icb_z'; fprintf('    baseline mode: ICB column -> Tg.bl_icb_z\n'); end
hasICB = any(~isnan(Tg.(icb_col_lme)));
fprintf('    LME ICB column   : Tg.%s  (n=%d non-missing)\n', icb_col_lme, sum(~isnan(Tg.(icb_col_lme))));

% --- Full covariate set ---
%   _raw : primary column name expected in Tg
%   _alt : fallback column name tried if _raw is missing
%   _lme : name used inside the LME formula (centered for continuous vars,
%          raw for sex — see loop below)
covNames_raw  = {'age_at_visit','SEX','EDUCYRS','updrs3_score','LEDD','depression_factor'};
covNames_alt  = {'AGE',         'sex','educ',   'UPDRS3_SCORE','ledd','depression'};
covNames_lme  = {'age_c',       'sex_num','edu_c','updrs3_c',  'ledd_c','depression_c'};

covPresent = false(numel(covNames_raw), 1);

for ci = 1:numel(covNames_raw)
    vn = covNames_raw{ci}; an = covNames_alt{ci}; ln = covNames_lme{ci};
    if ismember(vn, Tg.Properties.VariableNames), rawVals = double(safeCol(Tg, vn));
    elseif ismember(an, Tg.Properties.VariableNames), rawVals = double(safeCol(Tg, an));
    else, fprintf('      Covariate "%s" not found.\n', vn); continue;
    end
    covPresent(ci) = true;
    if strcmp(ln,'sex_num'), Tg.(ln) = rawVals;
    else, Tg.(ln) = rawVals - nanmean(rawVals);
    end
    fprintf('      %-20s -> %s  (n=%d)\n', vn, ln, sum(~isnan(rawVals)));
end
covStr = strjoin(covNames_lme(covPresent), ' + ');
hasCovariates = any(covPresent);
fprintf('    Covariate string: %s\n', covStr);

% --- Valid genetic groups ---
MIN_N = 20;
allGenes = unique(bl_g.gene_subgroup);
validGenes = strings(0,1);
for gi = 1:numel(allGenes)
    if sum(bl_g.gene_subgroup == allGenes(gi)) >= MIN_N
        validGenes(end+1) = allGenes(gi); %#ok<SAGROW>
    end
end
isSpo = validGenes == "Sporadic";
validGenes = ["Sporadic", sort(validGenes(~isSpo))];
validGenes(validGenes == "") = [];
nGG = numel(validGenes);
fprintf('    Valid groups (n>=%d): %s\n', MIN_N, strjoin(validGenes, ', '));

years = 0:max(Tg.year_num);
nYrs = numel(years);


% --- Gaussian smoothing setup ---
GAUSS_BW = 0.20;
gauss_sigma = GAUSS_BW * max(years);
smoothTraj = @(yrs, vals) gaussSmooth(yrs, vals, gauss_sigma);

% --- Colours ---
% Okabe-Ito subgroup palette. Positional mapping for the
% five groups present in the submitted figure (validGenes = Sporadic then
% alphabetical): Sporadic, GBA, LRRK2, PRKN/PINK1, SNCA. The old palette
% paired orange-red (GBA), yellow-green (PRKN/PINK1), and dark red (SNCA):
% simulated-deuteranopia dE for GBA1 vs PRKN/PINK1 is 61.9.
% Rows 6-8 are spares for any additional subgroup (e.g. LRRK2+GBA).
% PRKN/PINK1 bluish green -> Okabe-Ito sky blue.
% Bluish green and the Sporadic gray both render gray under deuteranopia
% (Machado 2009 simulation, CIE76 dE 11.3); sky blue is 26 L* units lighter
% than the gray (dE 44.0). Minimum over all ten subgroup pairs: normal 25.9,
% deuteranopia 24.7, protanopia 24.3, tritanopia 23.8 (LRRK2 vs PRKN/PINK1).
colPalette = [0.40 0.40 0.40;    % Sporadic   : gray (unchanged)
              0.80 0.40 0.00;    % GBA        : Okabe-Ito vermillion
              0.00 0.45 0.70;    % LRRK2      : Okabe-Ito blue
              0.35 0.70 0.90;    % PRKN/PINK1 : Okabe-Ito sky blue (was bluish green 0.00 0.60 0.50)
              0.10 0.10 0.10;    % SNCA       : near-black
              0.80 0.60 0.70;    % spare      : Okabe-Ito reddish purple
              0.00 0.60 0.50;    % spare      : Okabe-Ito bluish green (was sky blue, now used by PRKN/PINK1)
              0.95 0.90 0.25];   % spare      : Okabe-Ito yellow
geneColours = zeros(nGG, 3);
for gi = 1:nGG
    geneColours(gi,:) = colPalette(mod(gi-1, size(colPalette,1))+1, :);
end
% display names (italic gene symbols, GBA1, APOE e4). The data
% strings (validGenes, geneFactors_tested, table GeneFactor) are unchanged.
validGenesDisp = arrayfun(@(g) geneDisplayName(g), validGenes, 'UniformOutput', false);


%% =========================================================================
%  FIG 3  (3 x 3)
%  =========================================================================
fprintf('\n[1] Building Fig 3...\n');
fig3 = figure('Units','centimeters','Position',[1 1 FIG_W_CM FIG_H_CM],'Color',COL.bg);   % was [1 1 20 22]

% Layout: explicit axes positions in PRINTED cm (14.0 x 14.7),
% same 3 x 3 arrangement. Rows a-c share the subgroup rows (y = 1..5).
LAY.a  = posCM( 1.70, 11.08, 3.25, 3.00);
LAY.b  = posCM( 5.25, 11.08, 2.95, 3.00);   % rows shared with a (no y tick labels)
LAY.c  = posCM( 9.95, 11.08, 3.00, 3.00);   LAY.cbc = posCM(13.03, 11.08, 0.17, 3.00);
LAY.d  = posCM( 0.95,  6.63, 3.30, 3.00);
LAY.e  = posCM( 5.35,  6.63, 3.30, 3.00);
LAY.f  = posCM( 9.95,  6.63, 3.75, 3.00);
% row 3 in the columns of rows 1-2 (g under a, h under b/e, i under
% c/f); h shares g's rows, so g's colour bar is horizontal below g (LAY.cbg).
LAY.g  = posCM( 1.70,  2.03, 3.25, 3.00);   LAY.cbg = posCM( 0.40,  0.52, 1.20, 0.14);
LAY.h  = posCM( 5.30,  2.03, 3.35, 3.00);
LAY.i  = posCM( 9.95,  2.03, 3.75, 3.00);
TTL = struct('x', {0.10, 5.00, 8.60, 0.10, 4.60, 9.10, 0.10, 5.00, 8.60}, ...
             'y', {14.37, 14.37, 14.37, 9.91, 9.91, 9.91, 5.31, 5.31, 5.31});   % title baselines (printed cm), panels a..i

% =====================================================================
%  ROW 1: GENETIC LANDSCAPE  (panels a, b, c — unchanged)
% =====================================================================

% ----- (a) Subgroup frequencies + APOE overlay -----
ax_a = axes(fig3, 'Position', LAY.a); set(ax_a,'Color',COL.bg); hold(ax_a, 'on');   % was subplot(3,3,1)
nPerGroup = zeros(nGG,1); nAPOE_pos = zeros(nGG,1);
for gi = 1:nGG
    idx = bl_g.gene_subgroup == validGenes(gi);
    nPerGroup(gi) = sum(idx);
    nAPOE_pos(gi) = sum(bl_g.apoe4_bin(idx) == 1);
end

% --- Stats: is APOE e4 prevalence comparable across gene subgroups? -------
%   Omnibus chi-square of e4 carrier status (apoe4_bin) x gene subgroup.
%   If any expected cell < 5 the chi-square approximation is unreliable, so a
%   label-permutation p (Freeman-Halton-style exact test for the 2xk table)
%   is added. RNG is saved/restored so no other seeded analysis is perturbed.
okAE   = ~isnan(bl_g.apoe4_bin) & ismember(bl_g.gene_subgroup, validGenes);
aeVec  = double(bl_g.apoe4_bin(okAE));
sgVec  = removecats(categorical(bl_g.gene_subgroup(okAE)));
[obsAE, chi2AE, pAE] = crosstab(aeVec, sgVec);
dfAE   = (size(obsAE,1)-1) * (size(obsAE,2)-1);
expAE  = sum(obsAE,2) * sum(obsAE,1) / sum(obsAE(:));
isSparseAE = any(expAE(:) < 5);
fprintf('\n  [APOE e4 x subgroup] e4+ prevalence (APOE-known denominator):\n');
for gi = 1:nGG
    idxg = bl_g.gene_subgroup == validGenes(gi);
    nk   = sum(idxg & ~isnan(bl_g.apoe4_bin));
    np   = sum(bl_g.apoe4_bin(idxg) == 1);
    fprintf('      %-10s %5.1f%% (%d/%d)\n', char(validGenes(gi)), 100*np/max(nk,1), np, nk);
end
fprintf('      Omnibus chi2(%d) = %.2f, p = %.3f%s\n', dfAE, chi2AE, pAE, ...
    ternary(isSparseAE, '   [low expected counts -> use permutation p]', ''));
if isSparseAE
    sPrev = rng; rng(1);   % local + reproducible; restored below
    [~,~,gCode] = unique(sgVec);
    nGsz = accumarray(gCode,1); Nae = numel(aeVec); mPos = sum(aeVec);
    if mPos > 0 && mPos < Nae
        ePos = nGsz*mPos/Nae; eNeg = nGsz*(Nae-mPos)/Nae;
        chi2fun = @(xp) sum((xp-ePos).^2./ePos) + sum(((nGsz-xp)-eNeg).^2./eNeg);
        statObs = chi2fun(accumarray(gCode, aeVec));
        nPerm = 10000; geCount = 0;
        for pp = 1:nPerm
            if chi2fun(accumarray(gCode, aeVec(randperm(Nae)))) >= statObs - 1e-9, geCount = geCount + 1; end
        end
        fprintf('      Permutation (2xk exact) p = %.4f  (%d perms)\n', (geCount+1)/(nPerm+1), nPerm);
    end
    rng(sPrev);
end

% the e4+ percentage uses the GENOTYPED denominator
% (APOE-known participants), exactly as the log block above and the Results
% text; the bar is split e4+ / e4- / not genotyped so the drawn proportion
% matches the printed one. Counts only; no test is changed.
nAPOE_known = zeros(nGG,1);
for gi = 1:nGG
    nAPOE_known(gi) = sum(bl_g.gene_subgroup == validGenes(gi) & ~isnan(bl_g.apoe4_bin));
end
nAPOE_neg  = nAPOE_known - nAPOE_pos;
nAPOE_miss = nPerGroup - nAPOE_known;
hBarA = barh(ax_a, (1:nGG)', [nAPOE_pos, nAPOE_neg, nAPOE_miss], 0.6, 'stacked', 'EdgeColor', 'none');
hBarA(1).FaceColor = [0.80 0.40 0.00]; hBarA(1).FaceAlpha = 0.85;          % e4+ : Okabe-Ito vermillion (as before)
hBarA(2).FaceColor = [0.75 0.75 0.75];                                  % e4- : gray (as before)
hBarA(3).FaceColor = [1 1 1]; hBarA(3).EdgeColor = [0.55 0.55 0.55]; hBarA(3).LineWidth = 0.4;   % not genotyped
xLabA = 1250;   % label column (data units) right of the longest bar
for gi = 1:nGG
    text(ax_a, xLabA, gi, {sprintf('n = %d', nPerGroup(gi)), ...
        sprintf('%s4+ %d%% of %d', EPS_, round(100*nAPOE_pos(gi)/max(nAPOE_known(gi),1)), nAPOE_known(gi))}, ...
        'FontSize', FS_ANNOT, 'VerticalAlignment','middle', 'Color', COL.txt);
    fprintf('    panel a label: %-10s n = %d; e4+ %d/%d genotyped = %.1f%% (printed %d%%)\n', char(validGenes(gi)), ...
        nPerGroup(gi), nAPOE_pos(gi), nAPOE_known(gi), 100*nAPOE_pos(gi)/max(nAPOE_known(gi),1), ...
        round(100*nAPOE_pos(gi)/max(nAPOE_known(gi),1)));
end
set(ax_a,'YTick',1:nGG,'YTickLabel',validGenesDisp,'FontSize',FS_TICK,'YDir','reverse',...
    'Box','off','TickDir','out','XGrid','off','GridColor',COL.grid, 'YLim',[0.5 nGG+0.5], ...
    'XColor','none');   % n is printed, so the count axis is not drawn
xlim(ax_a, [0, 2800]);   % bars use the left 45%; labels sit in the right-hand column
lg = legend(ax_a, hBarA, {['\itAPOE\rm ' EPS_ '4+'], [EPS_ '4' MINUS_], 'not genotyped'}, 'FontSize', LEG_FS, ...
    'Box', 'off', 'Orientation', 'horizontal');
lg.ItemTokenSize = [8 6];
drawnow; lgP = lg.Position; lg.Position = [0.30/PRINT_W_CM, LAY.a(2) - lgP(4) - 0.05/PRINT_H_CM, lgP(3), lgP(4)];   % below panel a, off the data
aNoteY = (LAY.a(2) - lgP(4)) * PRINT_H_CM - 0.08;   % printed cm; the denominator note is drawn under the legend (overlay, below)
TTL(1).s = 'Subgroup Composition';

% ----- (b) Apathy & ICB prevalence -----
ax_b = axes(fig3, 'Position', LAY.b); set(ax_b,'Color',COL.bg); hold(ax_b, 'on');   % was subplot(3,3,2)
prevAp=nan(nGG,1); prevICB=nan(nGG,1);
ciAp_lo=nan(nGG,1); ciAp_hi=nan(nGG,1); ciICB_lo=nan(nGG,1); ciICB_hi=nan(nGG,1);
for gi = 1:nGG
    if hasGroup4
        g4 = bl_g.group4(bl_g.gene_subgroup==validGenes(gi)); g4=g4(~isnan(g4)); nv=numel(g4);
        if nv>0
            pA=sum(g4==2|g4==4)/nv; pI=sum(g4==3|g4==4)/nv;
            prevAp(gi)=pA*100; prevICB(gi)=pI*100;
            [ciAp_lo(gi),ciAp_hi(gi)]=wilsonCI(pA,nv);
            [ciICB_lo(gi),ciICB_hi(gi)]=wilsonCI(pI,nv);
        end
    end
end
bw=0.28;
b1=barh(ax_b,(1:nGG)'-bw/2,prevAp,bw,'FaceColor',COL.apathy,'EdgeColor','none','FaceAlpha',0.8);
b2=barh(ax_b,(1:nGG)'+bw/2,prevICB,bw,'FaceColor',COL.icb,'EdgeColor','none','FaceAlpha',0.8);
for gi=1:nGG
    if ~isnan(ciAp_lo(gi)), plot(ax_b,[ciAp_lo(gi) ciAp_hi(gi)]*100,[gi-bw/2 gi-bw/2],'-','Color',COL.apathy*0.6,'LineWidth',0.8); end
    if ~isnan(ciICB_lo(gi)), plot(ax_b,[ciICB_lo(gi) ciICB_hi(gi)]*100,[gi+bw/2 gi+bw/2],'-','Color',COL.icb*0.6,'LineWidth',0.8); end
end
% rows shared with panel a (same y positions), so no y tick
% labels; x axis 0-100% leaves the top-right corner empty for the legend
% and the test annotation (was southeast, over the SNCA bars).
set(ax_b,'YTick',1:nGG,'YTickLabel',[],'FontSize',FS_TICK,'YDir','reverse','YLim',[0.5 nGG+0.5],...
    'Box','off','TickDir','out','XGrid','on','GridColor',COL.grid, 'XLim',[0 100], 'XTick',0:20:100, ...
    'XTickLabelRotation',0);
xlabel(ax_b, 'Prevalence (%)','FontSize',FS_LABEL,'FontWeight','bold');
lg=legend(ax_b,[b1 b2],{'Apathy','ICB'},'FontSize',LEG_FS,'Box','off','Location','northeast'); lg.ItemTokenSize=LEG_TOKEN;
TTL(2).s = 'Apathy & ICB Prevalence';
if hasGroup4
    g4_bl=bl_g.group4(ismember(bl_g.gene_subgroup,validGenes));
    sg_bl=bl_g.gene_subgroup(ismember(bl_g.gene_subgroup,validGenes));
    ok=~isnan(g4_bl);
    % REVISED (bugfix): crosstab's 2nd output is the chi-square STATISTIC,
    % not the p-value; the submitted figure therefore printed the impossible
    % "chi^2 p=5.249". Capture the 3rd output (p) and print both values.
    % the 1st output (the table) is kept only to print df, and
    % the annotation names the test (apathy prevalence across subgroups).
    % the test is printed to the console only; the in-panel
    % annotation was removed at the author's request.
    if sum(ok)>20, [tblAp,chi2Ap,pChi]=crosstab(double(g4_bl(ok)==2|g4_bl(ok)==4),sg_bl(ok));
        dfAp = (size(tblAp,1)-1) * (size(tblAp,2)-1);
        fprintf('    panel b: apathy x subgroup chi2(%d) = %.2f, p = %.3f (not drawn)\n', dfAp, chi2Ap, pChi);
    end
end

% ----- (c) Co-occurrence heatmap -----
ax_c = axes(fig3, 'Position', LAY.c); set(ax_c,'Color',COL.bg);   % was subplot(3,3,3)
coMatrix=nan(nGG,4);
for gi=1:nGG
    idx=bl_g.gene_subgroup==validGenes(gi); n=sum(idx);
    if n>0
        ae=bl_g.apoe4_bin(idx); ab=bl_g.apathy_bin(idx);
        coMatrix(gi,:) = [sum(ae==0&ab==0,'omitnan'), sum(ae==0&ab==1,'omitnan'), ...
                          sum(ae==1&ab==0,'omitnan'), sum(ae==1&ab==1,'omitnan')] * 100 / n;
    end
end
% participants with APOE genotype or baseline
% apathy unavailable are shown in a fifth "missing" column (percent of the
% same whole-subgroup denominator), so each row sums to 100%. The four
% existing cells are computed exactly as before.
missPct = nan(nGG,1);
for gi=1:nGG
    idx=bl_g.gene_subgroup==validGenes(gi); n=sum(idx);
    if n>0
        missPct(gi) = 100 * sum(isnan(bl_g.apoe4_bin(idx)) | isnan(bl_g.apathy_bin(idx))) / n;
        fprintf('    panel c row %-10s: %.1f / %.1f / %.1f / %.1f / missing %.1f (sum %.1f%%)\n', char(validGenes(gi)), ...
            coMatrix(gi,1), coMatrix(gi,2), coMatrix(gi,3), coMatrix(gi,4), missPct(gi), sum(coMatrix(gi,:)) + missPct(gi));
    end
end
imagesc(ax_c, coMatrix);
colormap(ax_c,[linspace(1,COL.accent(1),64)',linspace(1,COL.accent(2),64)',linspace(1,COL.accent(3),64)']);
cb=colorbar(ax_c); cb.Label.String='% of subgroup'; cb.FontSize=FS_TICK; cb.Label.FontSize=PT(7);
cb.Position = LAY.cbc; ax_c.Position = LAY.c;
xMiss = 5.12;   % "missing" column set slightly apart from the colour-scaled cells
set(ax_c,'XLim',[0.5 xMiss+0.5],'XTick',[1:4 xMiss],'XTickLabel',{[EPS_ '4' MINUS_ '/apathy' MINUS_], [EPS_ '4' MINUS_ '/apathy+'], ...
    [EPS_ '4+/apathy' MINUS_], [EPS_ '4+/apathy+'], 'missing'},...
    'YTick',1:nGG,'YTickLabel',validGenesDisp,'FontSize',FS_TICK,'Box','off','TickDir','out');
xtickangle(ax_c,30);   % 30 degrees (less height than 45)
hold(ax_c, 'on');
for gi=1:nGG
    rectangle(ax_c, 'Position', [xMiss-0.5, gi-0.5, 1, 1], 'FaceColor', [1 1 1], 'EdgeColor', [0.70 0.70 0.70], 'LineWidth', 0.5);
    if ~isnan(missPct(gi))
        text(ax_c, xMiss, gi, sprintf('%.0f%%', missPct(gi)), 'HorizontalAlignment','center', 'FontSize', FS_ANNOT, ...
            'Color', [0.35 0.35 0.35], 'FontAngle', 'italic');
    end
end
for gi=1:nGG, for ci=1:4
    if ~isnan(coMatrix(gi,ci))
        tc=ternary(coMatrix(gi,ci)>35,[1 1 1],COL.txt);
        text(ax_c,ci,gi,sprintf('%.0f%%',coMatrix(gi,ci)),'HorizontalAlignment','center','FontSize',FS_ANNOT,'Color',tc,'FontWeight','bold');
    end
end, end
TTL(3).s = ['Co-occurrence: Gene ' TIMES_ ' \itAPOE\rm ' TIMES_ ' Apathy'];


% =====================================================================
%  ROW 2: APATHY & ICB TRAJECTORIES (within-subject SE + Gaussian smooth)
% =====================================================================

% ----- (d) Apathy trajectory by subgroup -----
ax_d = axes(fig3, 'Position', LAY.d); set(ax_d,'Color',COL.bg); hold(ax_d, 'on');   % was subplot(3,3,4)
phG_d = gobjects(nGG,1);
for gi = 1:nGG
    rowIdx = Tg.gene_subgroup == validGenes(gi);
    phG_d(gi) = plotCPTraj(ax_d, Tg, rowIdx, 'year_num', 'apathy_composite', ...
        geneColours(gi,:), nYrs, gauss_sigma);
end
% no per-panel legend (it sat on the GBA1 band and SNCA curve,
% F5-18); one subgroup key for d and e is drawn below the row (after e).
formatAxis(ax_d,COL); set(ax_d,'FontSize',FS_TICK); xlim(ax_d,[-0.3 max(years)+0.3]);
xlabel(ax_d,'Year','FontSize',FS_LABEL,'FontWeight','bold'); ylabel(ax_d,'Apathy Composite','FontSize',FS_LABEL,'FontWeight','bold');
TTL(4).s = 'Apathy Trajectory by Subgroup';

% LME for apathy trajectories
fprintf('    Testing apathy trajectories (LME, z-scored, full covariates)...\n');
testGenes_noSpo = validGenes(validGenes ~= "Sporadic");
for gi = 1:numel(testGenes_noSpo)
    idx = Tg.gene_subgroup == testGenes_noSpo(gi) | Tg.gene_subgroup == "Sporadic";
    T_sub = Tg(idx,:); T_sub.gene_flag = double(T_sub.gene_subgroup == testGenes_noSpo(gi));
    ok = ~isnan(T_sub.apathy_z) & ~isnan(T_sub.year_num); T_sub = T_sub(ok,:);
    if height(T_sub) < 50, continue; end
        if hasCovariates
            f1_ap = ['apathy_z ~ year_num * gene_flag + ' covStr ' + (year_num|PATNO)'];
            f2_ap = ['apathy_z ~ year_num * gene_flag + ' covStr ' + (1|PATNO)'];
        else
            f1_ap = 'apathy_z ~ year_num * gene_flag + (year_num|PATNO)';
            f2_ap = 'apathy_z ~ year_num * gene_flag + (1|PATNO)';
        end
        [lme_tmp,~] = fitLME_safe(T_sub, f1_ap, f2_ap);
        if ~isempty(lme_tmp)
            c = lme_tmp.Coefficients; ii = strcmp(c.Name,'year_num:gene_flag');
            if any(ii), fprintf('      %s: beta=%.4f, p=%.4f\n', testGenes_noSpo(gi), c.Estimate(ii), c.pValue(ii)); end
        end
end

% LME for ICB (QUIP) trajectories — Poisson on raw counts.
% QUIP is a count of endorsed items (>=0 integer); a Poisson GLME with log link
% is the natural family when QUIP is the outcome (Gaussian LME on z-scored QUIP
% would assume normality and homoscedasticity that count data do not satisfy).
if hasQUIP
    fprintf('    Testing ICB trajectories (GLME Poisson, raw QUIP counts, full covariates)...\n');
    for gi = 1:numel(testGenes_noSpo)
        idx = Tg.gene_subgroup == testGenes_noSpo(gi) | Tg.gene_subgroup == "Sporadic";
        T_sub = Tg(idx,:);
        T_sub.gene_flag = double(T_sub.gene_subgroup == testGenes_noSpo(gi));
        T_sub.quip = round(double(safeCol(T_sub,'quip')));   % enforce integer counts
        T_sub.quip(T_sub.quip < 0) = NaN;
        ok = ~isnan(T_sub.quip) & ~isnan(T_sub.year_num);
        T_sub = T_sub(ok,:);
        if height(T_sub) < 50, continue; end
        if hasCovariates
            f1_ic = ['quip ~ year_num * gene_flag + ' covStr ' + (year_num|PATNO)'];
            f2_ic = ['quip ~ year_num * gene_flag + ' covStr ' + (1|PATNO)'];
        else
            f1_ic = 'quip ~ year_num * gene_flag + (year_num|PATNO)';
            f2_ic = 'quip ~ year_num * gene_flag + (1|PATNO)';
        end
        [glme_tmp,~] = fitGLME_safe(T_sub, f1_ic, f2_ic, 'Poisson', 'log');
        if ~isempty(glme_tmp)
            c = glme_tmp.Coefficients;
            ii = strcmp(c.Name,'year_num:gene_flag');
            if any(ii)
                fprintf('      %s: beta=%.4f, p=%.4f (Poisson)\n', ...
                    testGenes_noSpo(gi), c.Estimate(ii), c.pValue(ii));
            end
        end
    end
end

% ----- (e) ICB trajectory by subgroup -----
ax_e = axes(fig3, 'Position', LAY.e); set(ax_e,'Color',COL.bg); hold(ax_e, 'on');   % was subplot(3,3,5)
phG_e = gobjects(nGG,1);
if hasQUIP
    for gi = 1:nGG
        rowIdx = Tg.gene_subgroup == validGenes(gi);
        phG_e(gi) = plotCPTraj(ax_e, Tg, rowIdx, 'year_num', 'quip', ...
            geneColours(gi,:), nYrs, gauss_sigma);
    end
    eYLabel = 'QUIP Total';
elseif hasGroup4
    % For proportion-based ICB, within-subject SE doesn't directly apply
    % (it's a rate not a continuous DV per subject), so use standard approach
    for gi = 1:nGG
        tSub = Tg(Tg.gene_subgroup == validGenes(gi),:);
        if numel(unique(tSub.PATNO)) < 10, continue; end
        mY=nan(nYrs,1);
        for yi=1:nYrs
            g4=tSub.group4_num(tSub.year_num==years(yi)); g4=g4(~isnan(g4));
            if numel(g4)>=5, mY(yi)=100*sum(g4==3|g4==4)/numel(g4); end
        end
        vYrs=~isnan(mY);
        if sum(vYrs)>=2
            yrsV=years(vYrs); mV=mY(vYrs);
            yrsD=linspace(min(yrsV),max(yrsV),200);
            mSm=smoothTraj(yrsD, interp1(yrsV,mV,yrsD,'pchip'));
            phG_e(gi)=plot(ax_e,yrsD,mSm,'-','Color',geneColours(gi,:),'LineWidth',1.5);
            plot(ax_e,yrsV,mV,'.','Color',geneColours(gi,:),'MarkerSize',6);
        end
    end
    eYLabel = 'ICB Prevalence (%)';
else
    eYLabel = '';
end
formatAxis(ax_e,COL); set(ax_e,'FontSize',FS_TICK); xlim(ax_e,[-0.3 max(years)+0.3]);
xlabel(ax_e,'Year','FontSize',FS_LABEL,'FontWeight','bold'); ylabel(ax_e,eYLabel,'FontSize',FS_LABEL,'FontWeight','bold');
TTL(5).s = 'ICB Trajectory by Subgroup';
% one subgroup key for panels d and e, centred below the row
% (the per-panel legends sat on the SNCA band and curves).
vph=phG_d(isgraphics(phG_d)); vn=validGenesDisp(isgraphics(phG_d));
if ~isempty(vph)
    lgDE = legend(ax_d, vph, vn, 'Orientation','horizontal', 'FontSize', LEG_FS, 'Box', 'off');
    lgDE.ItemTokenSize = LEG_TOKEN;
    drawnow; lgP = lgDE.Position; xcDE = (LAY.d(1) + LAY.e(1) + LAY.e(3)) / 2;
    lgDE.Position = [xcDE - lgP(3)/2, LAY.d(2) - 0.70/PRINT_H_CM - lgP(4), lgP(3), lgP(4)];
end

% ----- (f) Apathy–ICB relationship by subgroup -----
% Layout:
%   Main: scatter + per-subgroup regression lines
%   Top margin:   apathy KDE per subgroup + mean±SE bars
%   Right margin: ICB KDE per subgroup + mean±SE bars
%   Single legend in bottom-left of scatter (away from densities)
ax_f = axes(fig3, 'Position', LAY.f); set(ax_f,'Color',COL.bg); hold(ax_f, 'on');   % was subplot(3,3,6)

if hasQUIP
    icbVar = 'quip';
    icbLabel = 'QUIP Total';
else
    icbVar = '';
    icbLabel = '';
end

allAp = double(safeCol(bl_g, 'apathy_composite'));

if ~isempty(icbVar)
    allIcb = double(safeCol(bl_g, icbVar));
    
    % --- Compute per-group stats first (needed for layout) ---
    grpMeanAp = nan(nGG,1); grpSEAp = nan(nGG,1);
    grpMeanIcb= nan(nGG,1); grpSEIcb= nan(nGG,1);
    rhoPerGroup = nan(nGG,1); pPerGroup = nan(nGG,1);
    regCoefs = nan(nGG, 2);  % [intercept, slope]
    
    for gi = 1:nGG
        idx = bl_g.gene_subgroup == validGenes(gi);
        xv = allAp(idx); yv = allIcb(idx);
        ok = ~isnan(xv) & ~isnan(yv); xv = xv(ok); yv = yv(ok);
        if numel(xv) < 10, continue; end
        grpMeanAp(gi) = mean(xv); grpSEAp(gi) = std(xv)/sqrt(numel(xv));
        grpMeanIcb(gi)= mean(yv); grpSEIcb(gi)= std(yv)/sqrt(numel(yv));
        [rhoPerGroup(gi), pPerGroup(gi)] = corr(xv, yv, 'Type', 'Spearman');
        regCoefs(gi,:) = [ones(numel(xv),1), xv] \ yv;
    end
    
        % FDR-correct per-subgroup correlations
    validRho = ~isnan(pPerGroup);
    pPerGroup_FDR = nan(nGG, 1);
    pPerGroup_FDR(validRho) = fdr_bh(pPerGroup(validRho));
    fprintf('    Apathy-ICB Spearman correlations (FDR-corrected):\n');
    for gi = 1:nGG
        if ~isnan(rhoPerGroup(gi))
            fprintf('      %-15s  rho=%.3f  p_unc=%.4f  p_FDR=%.4f\n', ...
                validGenes(gi), rhoPerGroup(gi), pPerGroup(gi), pPerGroup_FDR(gi));
        end
    end
    % --- Layout geometry ---
    % X: scatter from 0 to data max, then gap, then ICB density + bars
    % Y: scatter from 0 to data max, then gap, then apathy density + bars
    xDataMax = max(allAp(~isnan(allAp))) * 1.05;
    yDataMax = max(allIcb(~isnan(allIcb))) * 1.05;
    xDataMin = 0;
    yDataMin = min(0, min(allIcb(~isnan(allIcb))));
    
    xRange = xDataMax - xDataMin;
    yRange = yDataMax - yDataMin;
    
    % Band sizes as fractions of data range
    gapFrac     = 0.03;
    barFrac     = 0.08;   % mean±SE bars
    densFrac    = 0.14;   % density curves
    
    % Y layout (top margin for apathy density)
    yGap     = yRange * gapFrac;
    yBarBase = yDataMax + yGap;
    yBarH    = yRange * barFrac;
    yDensBase= yBarBase + yBarH + yGap * 0.5;
    yDensH   = yRange * densFrac;
    
    % X layout (right margin for ICB density)
    xGap     = xRange * gapFrac;
    xBarBase = xDataMax + xGap;
    xBarW    = xRange * barFrac;
    xDensBase= xBarBase + xBarW + xGap * 0.5;
    xDensW   = xRange * densFrac;
    
    % --- Separator lines ---
    % plot(ax_f, [xDataMin xDataMax], [yDataMax yDataMax], '-', 'Color', COL.grid, 'LineWidth', 0.3);
    % plot(ax_f, [xDataMax xDataMax], [yDataMin yDataMax], '-', 'Color', COL.grid, 'LineWidth', 0.3);
    % 
    % --- Main scatter + regression lines ---
    phG_f = gobjects(nGG, 1);
    for gi = 1:nGG
        idx = bl_g.gene_subgroup == validGenes(gi);
        xv = allAp(idx); yv = allIcb(idx);
        ok = ~isnan(xv) & ~isnan(yv); xv = xv(ok); yv = yv(ok);
        if numel(xv) < 10, continue; end
        
        % Scatter with jitter
        jitX = 0.015 * randn(numel(xv),1) * xRange;
        jitY = 0.015 * randn(numel(yv),1) * yRange;
        scatter(ax_f, xv+jitX, yv+jitY, 6, geneColours(gi,:), ...
            'filled', 'MarkerFaceAlpha', 0.15, 'MarkerEdgeColor', 'none');
        
        % Regression line (within data range)
        xFit = linspace(max(xDataMin, min(xv)), min(xDataMax, max(xv)), 50);
        phG_f(gi) = plot(ax_f, xFit, regCoefs(gi,1) + regCoefs(gi,2)*xFit, '-', ...
            'Color', geneColours(gi,:), 'LineWidth', 1.6);
    end
    
    % Apathy threshold
    xline(ax_f, THRESH.apathy_comp, ':', 'Color', COL.apathy, 'LineWidth', 0.8, 'Alpha', 0.4);
    
    % --- TOP MARGIN: apathy density + mean±SE bars ---
    
    % Mean±SE bars (one per group, vertically stacked in bar band)
% Mean±SE bars (one per group, vertically stacked in top bar band)
    dodgeAmountY = yBarH * 1.5; 
    barYPositions = linspace(yBarBase - dodgeAmountY, yBarBase + dodgeAmountY, nGG);
    
    % Set a fixed cap size (in points) for the error bars
    capSize = 4; 
    
    for gi = 1:nGG
        if isnan(grpMeanAp(gi)), continue; end
        
        yp = barYPositions(gi);
        col = geneColours(gi,:);
        
        % Plot the horizontal error bar, caps, and mean dot simultaneously
        % Note the 'horizontal' argument to specify the error is on the X-axis
        errorbar(ax_f, grpMeanAp(gi), yp, grpSEAp(gi), 'horizontal', 'o', ...
            'Color', col, ...
            'LineWidth', 1, ...
            'CapSize', capSize, ...
            'MarkerSize', 5, ...
            'MarkerFaceColor', col, ...
            'MarkerEdgeColor', 'w');
    end
    
    % Density curves (overlaid above bars)
    % for gi = 1:nGG
    %     idx = bl_g.gene_subgroup == validGenes(gi);
    %     xv = allAp(idx); xv = xv(~isnan(xv));
    %     if numel(xv) < 10, continue; end
    % 
    %     [dens, xDens] = ksdensity(xv, 'NumPoints', 100);
    %     dens = dens / max(dens) * yDensH;  % normalise height
    % 
    %     fill(ax_f, [xDens fliplr(xDens)], ...
    %         [yDensBase + dens, repmat(yDensBase, 1, numel(xDens))], ...
    %         geneColours(gi,:), 'FaceAlpha', 0.20, 'EdgeColor', geneColours(gi,:), ...
    %         'LineWidth', 0.5, 'EdgeAlpha', 0.6);
    % end
    
    % --- RIGHT MARGIN: ICB density + mean±SE bars ---
    
% Mean±SE bars (one per group, horizontally stacked in bar band)
dodgeAmount = xBarW * 3; 
barXPositions = linspace(xBarBase - dodgeAmount, xBarBase + dodgeAmount, nGG);
% Set a fixed cap size (in points) for the error bars
capSize = 4; 

for gi = 1:nGG
    if isnan(grpMeanIcb(gi)), continue; end
    
    xp = barXPositions(gi);
    col = geneColours(gi,:);
    
    % Plot the error bar, caps, and mean dot simultaneously
    errorbar(ax_f, xp, grpMeanIcb(gi), grpSEIcb(gi), 'o', ...
        'Color', col, ...
        'LineWidth', 1, ...
        'CapSize', capSize, ...
        'MarkerSize', 5, ...
        'MarkerFaceColor', col, ...
        'MarkerEdgeColor', 'w');
end
    
    % Density curves (rotated 90°, plotted to the right)
    % for gi = 1:nGG
    %     idx = bl_g.gene_subgroup == validGenes(gi);
    %     yv = allIcb(idx); yv = yv(~isnan(yv));
    %     if numel(yv) < 10, continue; end
    % 
    %     [dens, yDens] = ksdensity(yv, 'NumPoints', 100);
    %     dens = dens / max(dens) * xDensW;  % normalise width
    % 
    %     % Rotated: y-values stay on y-axis, density extends rightward on x-axis
    %     fill(ax_f, [xDensBase + dens, repmat(xDensBase, 1, numel(yDens))], ...
    %         [yDens fliplr(yDens)], ...
    %         geneColours(gi,:), 'FaceAlpha', 0.20, 'EdgeColor', geneColours(gi,:), ...
    %         'LineWidth', 0.5, 'EdgeAlpha', 0.6);
    % end

    
    
    % --- Axis limits ---
    % extra headroom (0.35 x the QUIP range) above the marginal
    % band so the enlarged legend sits in empty space; y ticks limited to
    % the observed QUIP range (the band above holds the mean +/- SE markers).
    xlim(ax_f, [xDataMin, xDensBase + xDensW * 1.05]);
    ylim(ax_f, [yDataMin, yDensBase + yDensH * 1.05 + 0.35 * yRange]);
    set(ax_f, 'YTick', 0:floor(max(allIcb(~isnan(allIcb)))));

    % --- Single legend: regression lines only, with rho ---
    vph = phG_f(isgraphics(phG_f));
    vIdx = find(isgraphics(phG_f));
    legStr = cell(numel(vIdx), 1);
    for li = 1:numel(vIdx)
        gi = vIdx(li);
        if ~isnan(rhoPerGroup(gi))
            legStr{li} = sprintf('%s \\rho = %.2f%s', validGenesDisp{gi}, ...
                rhoPerGroup(gi), ternary(pPerGroup_FDR(gi)<0.05, '*', ''));
        else
            legStr{li} = validGenesDisp{gi};
        end
    end
    if ~isempty(vph)
        lg = legend(ax_f, vph, legStr, 'FontSize', LEG_FS, 'Box', 'off', 'Location', 'northeast');
        lg.ItemTokenSize = LEG_TOKEN;
        % nudge 0.15 cm right so it clears the marginal markers
        drawnow; lgP = lg.Position; lg.Position = [LAY.f(1) + LAY.f(3) - lgP(3) + 0.15/PRINT_W_CM, lgP(2), lgP(3), lgP(4)];
    end
    
else
    % Fallback: strip plots of apathy_composite by subgroup
    for gi = 1:nGG
        idx = bl_g.gene_subgroup == validGenes(gi);
        vals = allAp(idx); vals = vals(~isnan(vals));
        if isempty(vals), continue; end
        jit = 0.15 * (rand(numel(vals),1) - 0.5);
        scatter(ax_f, gi+jit, vals, 8, geneColours(gi,:), 'filled', 'MarkerFaceAlpha', 0.3);
        q = quantile(vals, [0.25 0.5 0.75]);
        plot(ax_f, [gi-0.15 gi+0.15], [q(2) q(2)], '-', 'Color', geneColours(gi,:), 'LineWidth', 2);
        rectangle('Position', [gi-0.18, q(1), 0.36, q(3)-q(1)], ...
            'EdgeColor', geneColours(gi,:), 'LineWidth', 0.8, 'Curvature', 0.1);
    end
    set(ax_f, 'XTick', 1:nGG, 'XTickLabel', validGenesDisp, 'FontSize', FS_TICK);
    xtickangle(ax_f, 30);
    yline(ax_f, THRESH.apathy_comp, ':', 'Color', COL.apathy, 'LineWidth', 0.8);
end

set(ax_f, 'FontSize', 6, 'Box', 'off', 'TickDir', 'out');
formatAxis(ax_f, COL); set(ax_f, 'FontSize', FS_TICK);
if ~isempty(icbVar)
    xlabel(ax_f, 'Apathy Composite', 'FontSize', FS_LABEL, 'FontWeight', 'bold');
    ylabel(ax_f, icbLabel, 'FontSize', FS_LABEL, 'FontWeight', 'bold');
else
    ylabel(ax_f, 'Apathy Composite', 'FontSize', FS_LABEL, 'FontWeight', 'bold');
end
TTL(6).s = 'Apathy–ICB by Subgroup';


% =====================================================================
%  ROW 3: GENETICS x MOTIVATION -> COGNITIVE DECLINE
% =====================================================================

% ----- LME Models -----
%  For each gene factor we fit THREE models on the SAME N (so AICs are
%  directly comparable):
%    apathy_only : cog_z ~ Year * Gene * Apathy           + cov + (Year|PATNO)
%    icb_only    : cog_z ~ Year * Gene * ICB              + cov + (Year|PATNO)
%    combined    : cog_z ~ Year * Gene * Apathy
%                       + Year * Gene * ICB               + cov + (Year|PATNO)
%  The "combined" model gives the gene x year x ICB effect *adjusted for*
%  gene x year x apathy (and vice versa) — this is the test the user
%  asked for ("ICB term in the current model").
fprintf('\n    Running LME models (z-scored, full covariates)...\n');
resultsCell = {}; geneFactors_tested = {};
modelTypes  = {'apathy_only','icb_only','combined'};
if ~hasICB, modelTypes = {'apathy_only'}; warning('ICB unavailable — running apathy_only only.'); end

% APOE model
fprintf('    Gene factor: APOE e4...\n');
idx_a = ~isnan(Tg.apoe4_bin) & ~isnan(Tg.(apathy_col_lme)) & ~isnan(Tg.cog_z);
if hasICB, idx_a = idx_a & ~isnan(Tg.(icb_col_lme)); end
Ta = Tg(idx_a,:);
Ta.gene_flag   = Ta.apoe4_bin;
Ta.apathy_flag = Ta.(apathy_col_lme);
if hasICB, Ta.icb_flag = Ta.(icb_col_lme); end
geneFactors_tested{end+1} = 'APOE_e4';
for mi = 1:numel(modelTypes)
    mt = modelTypes{mi};
    [f1,f2] = buildModelFormula(mt, false, hasCovariates, covStr);
    [lme,~] = fitLME_safe(Ta, f1, f2);
    if ~isempty(lme)
        resultsCell = appendLMEResults(resultsCell, lme, 'APOE_e4', mt, Ta);
        fprintf('      [%-11s] AIC=%.1f  N=%d  %d subj\n', mt, ...
            lme.ModelCriterion.AIC, height(Ta), numel(unique(Ta.PATNO)));
    else
        fprintf('      [%-11s] FIT FAILED.\n', mt);
    end
end

% Subgroup models
testGenes = validGenes(validGenes~="Sporadic");
for gi = 1:numel(testGenes)
    gN = testGenes(gi); fprintf('    Gene factor: %s...\n',gN);
    isV = Tg.gene_subgroup==gN; isS = Tg.gene_subgroup=="Sporadic";
    idx = (isV|isS) & ~isnan(Tg.(apathy_col_lme)) & ~isnan(Tg.cog_z) & ~isnan(Tg.apoe4_bin);
    if hasICB, idx = idx & ~isnan(Tg.(icb_col_lme)); end
    Ts = Tg(idx,:);
    Ts.gene_flag   = double(Ts.gene_subgroup==gN);
    Ts.apathy_flag = Ts.(apathy_col_lme);
    if hasICB, Ts.icb_flag = Ts.(icb_col_lme); end
    if numel(unique(Ts.PATNO(Ts.gene_flag==1))) < MIN_N
        fprintf('      Skipped (n<%d carriers after filtering).\n', MIN_N); continue;
    end
    geneFactors_tested{end+1} = char(gN); %#ok<SAGROW>
    for mi = 1:numel(modelTypes)
        mt = modelTypes{mi};
        [f1,f2] = buildModelFormula(mt, true, hasCovariates, covStr);
        [lme,~] = fitLME_safe(Ts, f1, f2);
        if ~isempty(lme)
            resultsCell = appendLMEResults(resultsCell, lme, char(gN), mt, Ts);
            fprintf('      [%-11s] AIC=%.1f  N=%d  %d subj\n', mt, ...
                lme.ModelCriterion.AIC, height(Ts), numel(unique(Ts.PATNO)));
        else
            fprintf('      [%-11s] FIT FAILED.\n', mt);
        end
    end
end

% =====================================================================
% Sporadic-only: Apathy x ICB x Time on cognitive decline
% =====================================================================
% Replicates the Figure 4 primary LME framework, restricted to the
% sporadic PD subgroup. Tests whether apathy and ICB jointly moderate
% cognitive decline (3-way interaction).
%   g ~ apathy*icb*time + age + sex + educ + nmiss + duration
%       + ledd + updrs3_score + depfactor + (1 + time | PATNO) + (1 | time)
fprintf('\n  Sporadic-only Apathy x ICB x Time LME...\n');

% --- Detect optional covariates (duration, nmiss); centre if present ---
durCands   = {'duration_yrs','duration','disease_duration','PD_duration','pd_duration'};
nmissCands = {'nmiss','n_miss','missing_n','NMISS'};
extraCovs  = {};
for ci = 1:numel(durCands)
    if ismember(durCands{ci}, Tg.Properties.VariableNames)
        v = double(safeCol(Tg, durCands{ci}));
        Tg.duration_c = v - nanmean(v);
        extraCovs{end+1} = 'duration_c'; %#ok<SAGROW>
        fprintf('    duration ("%s") found, centred as duration_c.\n', durCands{ci});
        break;
    end
end
for ci = 1:numel(nmissCands)
    if ismember(nmissCands{ci}, Tg.Properties.VariableNames)
        v = double(safeCol(Tg, nmissCands{ci}));
        Tg.nmiss_c = v - nanmean(v);
        extraCovs{end+1} = 'nmiss_c'; %#ok<SAGROW>
        fprintf('    nmiss ("%s") found, centred as nmiss_c.\n', nmissCands{ci});
        break;
    end
end
sporadicCovStr = covStr;
if ~isempty(extraCovs), sporadicCovStr = [sporadicCovStr ' + ' strjoin(extraCovs,' + ')]; end

% --- Subset to sporadic & non-missing predictors ---
idxSpo = Tg.gene_subgroup == "Sporadic" ...
       & ~isnan(Tg.(apathy_col_lme)) & ~isnan(Tg.(icb_col_lme)) ...
       & ~isnan(Tg.cog_z)    & ~isnan(Tg.year_num);
Tspo = Tg(idxSpo,:);
Tspo.time   = Tspo.year_num;
Tspo.apathy = Tspo.(apathy_col_lme);
Tspo.icb    = Tspo.(icb_col_lme);

if hasCovariates
    f1_spo = ['cog_z ~ apathy*icb*time + ' sporadicCovStr ...
              ' + (1 + time | PATNO) + (1 | time)'];
    f2_spo = ['cog_z ~ apathy*icb*time + ' sporadicCovStr ...
              ' + (1 | PATNO)'];
else
    f1_spo = 'cog_z ~ apathy*icb*time + (1 + time | PATNO) + (1 | time)';
    f2_spo = 'cog_z ~ apathy*icb*time + (1 | PATNO)';
end

fprintf('    Formula: %s\n', f1_spo);
fprintf('    N=%d observations, %d subjects\n', height(Tspo), numel(unique(Tspo.PATNO)));
[lme_spo,~] = fitLME_safe(Tspo, f1_spo, f2_spo);

if ~isempty(lme_spo)
    c = lme_spo.Coefficients; ci_spo = coefCI(lme_spo);
    keyTerms = {'apathy','icb','time','apathy:icb','apathy:time','icb:time','apathy:icb:time'};
    fprintf('    [Sporadic LME, AIC=%.1f]\n', lme_spo.ModelCriterion.AIC);
    fprintf('    %-20s  %+9s  %8s  %18s  %s\n','Term','beta','SE','95%% CI','p');
    for ti = 1:numel(keyTerms)
        ii = strcmp(c.Name, keyTerms{ti});
        if any(ii)
            fprintf('    %-20s  %+9.4f  %8.4f  [%+.3f, %+.3f]  %.4g\n', ...
                keyTerms{ti}, c.Estimate(ii), c.SE(ii), ci_spo(ii,1), ci_spo(ii,2), c.pValue(ii));
        end
    end

    % Save full coefficients table
    sporadicTable = table(c.Name, c.Estimate, c.SE, c.tStat, c.DF, c.pValue, ...
                         ci_spo(:,1), ci_spo(:,2), ...
                         repmat(height(Tspo),               height(c),1), ...
                         repmat(numel(unique(Tspo.PATNO)),  height(c),1), ...
                         repmat(lme_spo.ModelCriterion.AIC, height(c),1), ...
        'VariableNames', {'Term','Beta','SE','tStat','DF','pValue', ...
                          'CI_lo','CI_hi','N_obs','N_subj','AIC'});
    writetable(sporadicTable, fullfile(outputDir, ...
        sprintf('Table_LME_SporadicApathyICB_%s.csv', COG_VAR)));
    fprintf('    Saved: Table_LME_SporadicApathyICB_%s.csv\n', COG_VAR);
end

% Build table + FDR (now with Model column)
colNames={'GeneFactor','Model','Term','Beta','SE','tStat','DF','pValue','CI_lo','CI_hi','N_obs','N_subj','AIC'};
resultsTable=cell2table(resultsCell,'VariableNames',colNames);
resultsTable.Beta=round(resultsTable.Beta,4); resultsTable.SE=round(resultsTable.SE,4);
resultsTable.tStat=round(resultsTable.tStat,3); resultsTable.pValue=round(resultsTable.pValue,6);
resultsTable.CI_lo=round(resultsTable.CI_lo,4); resultsTable.CI_hi=round(resultsTable.CI_hi,4);
interactionTerms = contains(resultsTable.Term,':');
% FDR-correct *within each model family* (apathy_only / icb_only / combined).
% These three model families test overlapping hypotheses on the same data,
% so pooling them in one FDR pass would inflate the correction.
resultsTable.pFDR = nan(height(resultsTable),1);
mt_list = unique(resultsTable.Model);
for k = 1:numel(mt_list)
    mask = strcmp(resultsTable.Model, mt_list{k}) & interactionTerms;
    if any(mask)
        resultsTable.pFDR(mask) = fdr_bh(resultsTable.pValue(mask));
    end
end

nGFt = numel(geneFactors_tested);

% Significance matrix is built from the COMBINED model when available
% (so apathy and ICB 3-ways are mutually adjusted); else apathy_only.
matrixModel = ternary(hasICB, 'combined', 'apathy_only');

if hasICB
    intLabels={'gene_flag', 'year_num:gene_flag', ...
               'gene_flag:apathy_flag','year_num:gene_flag:apathy_flag', ...
               'gene_flag:icb_flag',   'year_num:gene_flag:icb_flag'};
    intDisplay={'Gene (main)','Gene x Year', ...
                'Gene x Apathy','Gene x Year x Apathy', ...
                'Gene x ICB','Gene x Year x ICB'};
    tw_apathy = 4;   % column for Gene x Year x Apathy
    tw_icb    = 6;   % column for Gene x Year x ICB
else
    intLabels={'gene_flag','year_num:gene_flag','gene_flag:apathy_flag','year_num:gene_flag:apathy_flag'};
    intDisplay={'Gene (main)','Gene x Year','Gene x Apathy','Gene x Year x Apathy'};
    tw_apathy = 4; tw_icb = NaN;
end
nIT = numel(intLabels);
betaMatrix=nan(nGFt,nIT); pMatrix=nan(nGFt,nIT); pFDRMatrix=nan(nGFt,nIT);
ciLoMatrix=nan(nGFt,nIT); ciHiMatrix=nan(nGFt,nIT);
for gi=1:nGFt, for ti=1:nIT
    idx = strcmp(resultsTable.GeneFactor,geneFactors_tested{gi}) & ...
          strcmp(resultsTable.Model,matrixModel) & ...
          strcmp(resultsTable.Term,intLabels{ti});
    if any(idx), r=find(idx,1);
        betaMatrix(gi,ti)=resultsTable.Beta(r); pMatrix(gi,ti)=resultsTable.pValue(r);
        pFDRMatrix(gi,ti)=resultsTable.pFDR(r);
        ciLoMatrix(gi,ti)=resultsTable.CI_lo(r); ciHiMatrix(gi,ti)=resultsTable.CI_hi(r);
    end
end, end

% ----- (g) Significance matrix: diverging beta -----
ax_g = axes(fig3, 'Position', LAY.g); set(ax_g,'Color',COL.bg);   % was subplot(3,3,7)
bp = betaMatrix; bp(isnan(bp)) = 0;

% Tight, robust colour range: 90th percentile of |beta|, capped to a
% sensible band so a single outlier doesn't desaturate the rest.
absVals = abs(bp(:)); absVals = absVals(absVals > 0);
if isempty(absVals)
    clim_max = 0.05;
else
    clim_max = quantile(absVals, 0.90);
    clim_max = max(0.02, min(0.10, clim_max));   % keep range in [0.02, 0.10]
end
clims = [-1 1] * clim_max;
imagesc(ax_g, bp, clims);
nCM=256; cmap_div=zeros(nCM,3);
for ci=1:nCM, t=(ci-1)/(nCM-1);
    if t<0.5, s=t/0.5; cmap_div(ci,:)=[0.12+0.88*s,0.30+0.70*s,0.85+0.15*s];
    else, s=(t-0.5)/0.5; cmap_div(ci,:)=[1.0,1.0-0.70*s,1.0-0.85*s]; end
end
colormap(ax_g, cmap_div);
% horizontal colour bar below the panel (h now sits where the
% vertical bar was); its label and the star key are drawn on the overlay axes.
cb = colorbar(ax_g, 'Location', 'southoutside'); cb.FontSize = FS_TICK;
cb.Ticks = [];   % the three tick labels (clims(1), 0, clims(2)) are drawn as text below the bar on the
                 % overlay axes; the colour bar's own labels inherited the 30-degree rotation of the x tick labels
cb.Position = LAY.cbg; ax_g.Position = LAY.g;
hold(ax_g, 'on');

% Cell annotations: print every beta value (small font), put significance
% stars in the upper-right corner so they don't collide with the number.
starsUsed = {};   % star symbols drawn, for the key under the colour bar
for gi = 1:nGFt
    for ti = 1:nIT
        if isnan(betaMatrix(gi,ti)), continue; end
        pF = pFDRMatrix(gi,ti);
        % Pick stars
        if     isnan(pF),       st='';
        elseif pF < 0.001,      st='***';
        elseif pF < 0.01,       st='**';
        elseif pF < 0.05,       st='*';
        elseif pF < 0.10,       st=char(8224);   % dagger for trend
        else,                   st='';
        end
        if ~isempty(st), starsUsed{end+1} = st; end %#ok<AGROW>
        % Text colour vs background (dark text on light cells, white on saturated)
        nB = (bp(gi,ti)-clims(1))/(clims(2)-clims(1));
        tc = ternary(nB > 0.30 && nB < 0.70, [0.15 0.15 0.15], [1 1 1]);
        % Beta number — every cell (nudged down when the cell carries stars)
        text(ax_g, ti, gi + ternary(isempty(st), 0, 0.12), sprintf('%.2f', betaMatrix(gi,ti)), ...
            'HorizontalAlignment','center','VerticalAlignment','middle', ...
            'FontSize', FS_ANNOT, 'Color', tc);
        % Stars in upper-right corner of cell
        if ~isempty(st)
            text(ax_g, ti+0.46, gi-0.46, st, ...
                'HorizontalAlignment','right','VerticalAlignment','top', ...
                'FontSize', FS_ANNOT, 'Color', tc, 'FontWeight','bold');
        end
        % Outline for FDR<0.05 (solid) or <0.10 (dashed)
        if ~isnan(pF) && pF < 0.05
            rectangle(ax_g, 'Position',[ti-0.5,gi-0.5,1,1], ...
                'EdgeColor',[0.1 0.1 0.1], 'LineWidth', 1.6);
        elseif ~isnan(pF) && pF < 0.10
            rectangle(ax_g, 'Position',[ti-0.5,gi-0.5,1,1], ...
                'EdgeColor',[0.4 0.4 0.4], 'LineWidth', 0.8, 'LineStyle','--');
        end
    end
end
% the title names the FDR family the code uses
% (pFDR computed per model type over all its interaction terms; for the
% combined model 7 terms x 5 gene factors = 35), not "within model".
nFDRfamily = sum(strcmp(resultsTable.Model, matrixModel) & interactionTerms);
fprintf('    panel g: BH-FDR family = %d interaction terms of the %s models\n', nFDRfamily, matrixModel);
geneFactorsDisp = cellfun(@geneDisplayName, geneFactors_tested, 'UniformOutput', false);
intDisplayFig = strrep(intDisplay, ' x ', [' ' TIMES_ ' ']);
set(ax_g, 'XTick', 1:nIT, 'XTickLabel', intDisplayFig, ...
          'YTick', 1:nGFt, 'YTickLabel', geneFactorsDisp, ...
          'FontSize', FS_TICK, 'Box', 'off', 'TickDir', 'out');
xtickangle(ax_g, 30);   % was 40 (less height at print size)
% the FDR family moved from the title to the key under the colour
% bar (the title must end before panel h's title, which now starts at 5.0 cm).
TTL(7).s = sprintf('LME \\beta_z (%s model)', strrep(matrixModel,'_','\_'));
starDefs = {'***', '< 0.001'; '**', '< 0.01'; '*', '< 0.05'; char(8224), '< 0.10'};
keyParts = {};
for kk = 1:size(starDefs, 1)
    if any(strcmp(starsUsed, starDefs{kk,1}))
        keyParts{end+1} = sprintf('%s %s%s', starDefs{kk,1}, ternary(isempty(keyParts), 'p_{FDR} ', ''), starDefs{kk,2}); %#ok<AGROW>
    end
end
gKeyStr = strjoin(keyParts, ', ');   % the FDR family (nFDRfamily terms) is stated in the caption and Methods
fprintf('    panel g key: %s (FDR across %d terms)\n', gKeyStr, nFDRfamily);

% ----- (h) Forest plot: gene x year x Apathy (squares) and gene x year x ICB (circles) -----
ax_h = axes(fig3, 'Position', LAY.h); set(ax_h,'Color',COL.bg); hold(ax_h, 'on');   % was subplot(3,3,8)
% rows shared with panel g (same gene order, as b shares the rows
% of a), so the genes are no longer sorted by their smallest p.
si = (1:nGFt)';
gnSorted = geneFactors_tested(si);

% Helper: pick colour from p
pickCol = @(p) ternary(~isnan(p)&&p<0.05, COL.accent, ...
                ternary(~isnan(p)&&p<0.10, [0.85 0.55 0.25], [0.55 0.55 0.55]));

dy = 0.18;   % vertical offset between Apathy and ICB markers
hApathyMarker = []; hIcbMarker = [];   % handles for the legend
pLab = cell(0, 4);   % p labels {y, text, colour, gene row}, drawn after the loop at the right end of the row
for gi = 1:nGFt
    g = si(gi);
    % Apathy 3-way: square marker, plotted slightly above the row centre
    bA = betaMatrix(g,tw_apathy); lA = ciLoMatrix(g,tw_apathy); hA = ciHiMatrix(g,tw_apathy); pA = pFDRMatrix(g,tw_apathy);
    if ~isnan(bA)
        c=pickCol(pA);
        plot(ax_h,[lA hA],[gi-dy gi-dy],'-','Color',c,'LineWidth',1);
        hMk = plot(ax_h,bA,gi-dy,'s','MarkerSize',6,'MarkerFaceColor',c,'MarkerEdgeColor','w','LineWidth',0.4);
        if isempty(hApathyMarker), hApathyMarker = hMk; end
        if ~isnan(pA) && pA<0.10
            % labelled, "< 0.001" below 0.001 (was an
            % unlabelled "0.000"). drawn after the loop at the
            % right end of the panel on this row (the left-end placement met
            % the marker of the row above once h shared the rows of g).
            pLab(end+1,:) = {gi-dy-0.06, ['p_{FDR} ' fmtP(pA)], c, g}; %#ok<AGROW>
        end
    end
    if hasICB
        bI = betaMatrix(g,tw_icb); lI = ciLoMatrix(g,tw_icb); hI = ciHiMatrix(g,tw_icb); pI = pFDRMatrix(g,tw_icb);
        if ~isnan(bI)
            c=pickCol(pI);
            plot(ax_h,[lI hI],[gi+dy gi+dy],'-','Color',c,'LineWidth',1);
            hMk = plot(ax_h,bI,gi+dy,'o','MarkerSize',6,'MarkerFaceColor',c,'MarkerEdgeColor','w','LineWidth',0.4);
            if isempty(hIcbMarker), hIcbMarker = hMk; end
            if ~isnan(pI) && pI<0.10
                pLab(end+1,:) = {gi+dy+0.06, ['p_{FDR} ' fmtP(pI)], c, g}; %#ok<AGROW>   % as the apathy label
            end
        end
    end
end
% x-limits from the CIs, widened on the right so that a p label
% (about 1.5 cm at print size) fits right of the markers of its row.
ciAll = [ciLoMatrix(:,tw_apathy); ciHiMatrix(:,tw_apathy)];
if hasICB, ciAll = [ciAll; ciLoMatrix(:,tw_icb); ciHiMatrix(:,tw_icb)]; end
ciAll = ciAll(~isnan(ciAll));
xlH = [min(ciAll) - 0.004, max(ciAll) + 0.004];
if ~isempty(pLab)
    lblFrac = 1.70 / (LAY.h(3) * PRINT_W_CM);            % label plus margins as a fraction of the axis width
    rowsLab = unique(cell2mat(pLab(:,4)));
    mk = betaMatrix(rowsLab, tw_apathy); if hasICB, mk = [mk; betaMatrix(rowsLab, tw_icb)]; end
    need = max(mk(~isnan(mk))) + 0.006;                    % right edge of the markers on the labelled rows
    xlH(2) = max(xlH(2), (need - lblFrac * xlH(1)) / (1 - lblFrac));
end
xlim(ax_h, xlH);
for k = 1:size(pLab, 1)
    text(ax_h, xlH(2) - 0.015 * diff(xlH), pLab{k,1}, pLab{k,2}, 'FontSize', FS_SCRIPT, 'Color', pLab{k,3}, ...
        'HorizontalAlignment', 'right', 'VerticalAlignment', 'middle');
end
fprintf('    panel h: x-limits [%.3f %.3f]; %d p label(s) at the right end of the row\n', xlH(1), xlH(2), size(pLab, 1));
xline(ax_h,0,'-','Color',[0.5 0.5 0.5],'LineWidth',0.8,'Alpha',0.5);
gnSortedDisp = cellfun(@geneDisplayName, gnSorted, 'UniformOutput', false);
set(ax_h,'YTick',1:nGFt,'YTickLabel',[],'FontSize',6,'Box','off','TickDir','out','YDir','reverse', ...
    'YLim',[0.5 nGFt+0.5]);   % rows and row labels shared with panel g
fprintf('    panel h rows (top to bottom, as panel g): %s\n', strjoin(gnSorted, ', '));
TTL(8).s = 'Three-way interactions';

% Legend: marker shape -> motivation factor (colour-neutral grey for clarity)
% moved below the axis label (it sat on the PRKN/PINK1 CI).
if hasICB && ~isempty(hApathyMarker) && ~isempty(hIcbMarker)
    % Build neutral-grey proxy markers so the legend reads shape, not significance
    hPxA = plot(ax_h, NaN, NaN, 's', 'MarkerSize',6, 'MarkerFaceColor',[0.45 0.45 0.45], 'MarkerEdgeColor','w');
    hPxI = plot(ax_h, NaN, NaN, 'o', 'MarkerSize',6, 'MarkerFaceColor',[0.45 0.45 0.45], 'MarkerEdgeColor','w');
    lg = legend(ax_h, [hPxA hPxI], {'Apathy','ICB'}, 'Orientation','horizontal', ...
        'FontSize', LEG_FS, 'Box','off');
    lg.ItemTokenSize = LEG_TOKEN;
    drawnow; lgP = lg.Position;
    lg.Position = [LAY.h(1) + LAY.h(3)/2 - lgP(3)/2, LAY.h(2) - 1.00/PRINT_H_CM - lgP(4), lgP(3), lgP(4)];   % centred under h
end
formatAxis(ax_h,COL); set(ax_h,'FontSize',FS_TICK);
xlabel(ax_h,'\beta_z  (Gene \times Year \times Motivation)','FontSize',FS_LABEL,'FontWeight','bold');

% ----- (i) Cognitive trajectories: strongest 3-way (apathy or ICB) -----
ax_i = axes(fig3, 'Position', LAY.i); set(ax_i,'Color',COL.bg); hold(ax_i, 'on');   % was subplot(3,3,9)

% Find smallest p across the available 3-ways and remember which one
[minP_apathy, mIdxA] = min(pFDRMatrix(:,tw_apathy));
if hasICB
    [minP_icb, mIdxI] = min(pFDRMatrix(:,tw_icb));
else
    minP_icb = NaN; mIdxI = NaN;
end
% Choose whichever motivation shows the smaller p (treating NaN as worst)
useICB = false;
if hasICB && (~isnan(minP_icb)) && (isnan(minP_apathy) || minP_icb < minP_apathy)
    useICB = true; minP = minP_icb; minIdx = mIdxI;
else
    minP = minP_apathy; minIdx = mIdxA;
end
% Fallback to gene x year main effect if no 3-way was significant
if isnan(minP) || minP >= 0.05
    [mp2, mi2] = min(pFDRMatrix(:,2));
    if ~isnan(mp2) && (isnan(minP) || mp2 < minP), minP = mp2; minIdx = mi2; useICB = false; end
end
bestGene = geneFactors_tested{minIdx};
motivLabel = ternary(useICB, 'ICB', 'Apathy');
fprintf('    Panel (i): %s x %s (p_FDR=%.4f)\n', bestGene, motivLabel, minP);

% Stratify visualisation by binary of whichever motivation is selected
% display labels spelled out ("no apathy"/"apathy"; was NoAp/Ap)
if useICB, motivBin = 'icb_bin';   shortNo = 'no ICB'; shortYes = 'ICB';
else,      motivBin = 'apathy_bin_disp'; shortNo = 'no apathy'; shortYes = 'apathy'; end

if strcmp(bestGene,'APOE_e4')
    isC=Tg.apoe4_bin==1; isNC=Tg.apoe4_bin==0;
    vO=~isnan(Tg.apoe4_bin)&~isnan(Tg.(motivBin))&~isnan(Tg.cog_z);
    cL=[EPS_ '4+']; ncL=[EPS_ '4' MINUS_];
    colNC = [0.40 0.40 0.40]; colC = [0.80 0.40 0.00];   % as panel a (e4- gray, e4+ vermillion)
else
    isC=Tg.gene_subgroup==string(bestGene); isNC=Tg.gene_subgroup=="Sporadic";
    vO=(isC|isNC)&~isnan(Tg.(motivBin))&~isnan(Tg.cog_z);
    cL=geneDisplayName(bestGene); ncL='Sporadic';
    colNC = geneColours(validGenes=="Sporadic",:); colC = geneColours(validGenes==string(bestGene),:);
end

strata={isNC&Tg.(motivBin)==0&vO, isNC&Tg.(motivBin)==1&vO, ...
        isC &Tg.(motivBin)==0&vO, isC &Tg.(motivBin)==1&vO};
% strata take the subgroup colours
% (was gray / orange / blue / purple, which reused the symptom-group colours
% and gave GBA1 the LRRK2 blue). Light tint + solid = no apathy; full
% subgroup colour + dashed = apathy. Simulated deuteranopia/protanopia:
% every pair CIE76 dE >= 29.
tintC = @(c) c + 0.5 * (1 - c);
sColors=[tintC(colNC); colNC; tintC(colC); colC];
sLines={'-','--','-','--'};
sNames={[ncL ', ' shortNo],[ncL ', ' shortYes],[cL ', ' shortNo],[cL ', ' shortYes]};
if strcmpi(STRATA_MODE, 'baseline') && ~useICB, sNames = cellfun(@(x) [x ' at baseline'], sNames, 'UniformOutput', false); end
phS=gobjects(4,1);

for si_=1:4
    rowIdx = strata{si_};
    nSub = numel(unique(Tg.PATNO(rowIdx)));
    if nSub < 5, continue; end
    phS(si_) = plotCPTraj(ax_i, Tg, rowIdx, 'year_num', 'cog_z', ...
        sColors(si_,:), nYrs, gauss_sigma, sLines{si_});
    % legend without the n; the n is logged here.
    fprintf('    panel i stratum %d: %s (n = %d)\n', si_, strrep(strrep(sNames{si_}, '\it', ''), '\rm', ''), nSub);
end

yline(ax_i,0,'--','Color',[0.7 0.7 0.7],'LineWidth',0.5);
formatAxis(ax_i,COL); set(ax_i,'FontSize',FS_TICK); xlim(ax_i,[-0.3 max(years)+0.5]);
if strcmp(COG_VAR, 'bifactor_g'), cogYLabel = 'Global cognition (\itg\rm, z)';   % was 'bifactor_g (z)'; 'z' (was 'z-score') so the label ends below the panel title
else, cogYLabel = strrep(COG_VAR_LABEL, '_', '\_'); end
xlabel(ax_i,'Year','FontSize',FS_LABEL,'FontWeight','bold'); ylabel(ax_i,cogYLabel,'FontSize',FS_LABEL,'FontWeight','bold');
TTL(9).s = sprintf('%s %s %s: p_{FDR} %s', geneDisplayName(bestGene), TIMES_, motivLabel, fmtP(minP));
% legend below the axis label (it sat on the GBA1-apathy curve,
% F5-18); tokens long enough to show the dashed strata.
vph=phS(isgraphics(phS)); vn=sNames(isgraphics(phS));
if ~isempty(vph)
    lg=legend(ax_i,vph,vn,'FontSize',LEG_FS,'Box','off'); lg.ItemTokenSize=[26 8];
    drawnow; lgP = lg.Position;
    lg.Position = [LAY.i(1) + LAY.i(3)/2 - lgP(3)/2, LAY.i(2) - 0.70/PRINT_H_CM - lgP(4), lgP(3), lgP(4)];   % centred under i
end

% sgtitle('Fig 3: Genetics, Motivation & Cognitive Decline','FontSize',10,'FontWeight','bold','Color',COL.accent);
% panel titles drawn on a transparent overlay axes in printed cm,
% left-aligned at each column start, so neighbouring titles cannot collide.
axTTL = axes(fig3, 'Position', [0 0 1 1], 'Visible', 'off', 'XLim', [0 PRINT_W_CM], 'YLim', [0 PRINT_H_CM], 'HitTest', 'off');
panelLetters = 'abcdefghi';
for k = 1:numel(TTL)
    if isempty(TTL(k).s), continue; end
    text(axTTL, TTL(k).x, TTL(k).y, ['\bf' panelLetters(k) '\rm  ' TTL(k).s], 'FontSize', FS_TITLE, 'Color', COL.txt, ...
        'HorizontalAlignment', 'left', 'VerticalAlignment', 'baseline');
end
text(axTTL, 0.30, aNoteY, ['\itAPOE\rm ' EPS_ '4+ %: carriers / genotyped'], 'FontSize', FS_ANNOT, 'Color', COL.txt, ...
    'HorizontalAlignment', 'left', 'VerticalAlignment', 'top');   % panel a denominator
% panel g colour-bar label above the bar and the star key on the
% bottom line, right of the bar's tick labels (printed cm; the rotated x tick
% labels of g reach down to about 0.65 cm at 1.7 to 3.3 cm).
text(axTTL, (LAY.cbg(1) + LAY.cbg(3)/2)*PRINT_W_CM, (LAY.cbg(2) + LAY.cbg(4))*PRINT_H_CM + 0.05, '\beta_z', ...
    'FontSize', FS_SCRIPT, 'Color', COL.txt, 'HorizontalAlignment', 'center', 'VerticalAlignment', 'bottom');
cbTickX = LAY.cbg(1)*PRINT_W_CM + [0 0.5 1] * LAY.cbg(3)*PRINT_W_CM;
cbTickS = {sprintf('%g', clims(1)), '0', sprintf('%g', clims(2))};
for kk = 1:3
    text(axTTL, cbTickX(kk), LAY.cbg(2)*PRINT_H_CM - 0.05, cbTickS{kk}, 'FontSize', FS_TICK, 'Color', COL.txt, ...
        'HorizontalAlignment', 'center', 'VerticalAlignment', 'top');
end
if ~isempty(gKeyStr)
    text(axTTL, 2.10, 0.02, gKeyStr, 'FontSize', FS_SCRIPT, 'Color', COL.txt, ...
        'HorizontalAlignment', 'left', 'VerticalAlignment', 'bottom');
end
if exist('FIG_FONT_AUDIT', 'var') && isequal(FIG_FONT_AUDIT, true)
    auditFigureText(fig3, PRINT_W_CM / FIG_W_CM, fullfile(outputDir, sprintf('FontAudit_Fig3_Genetics_Motivation_Cognition_%s.csv', COG_VAR)));
end
exportFigure(fig3, fullfile(outputDir,sprintf('Fig3_Genetics_Motivation_Cognition_%s',COG_VAR)));
fprintf('\n    -> Fig3_Genetics_Motivation_Cognition_%s saved\n', COG_VAR);


%% ====  TABLES  ===========================================================
fprintf('\n[Tables] Saving...\n');
writetable(resultsTable, fullfile(outputDir,sprintf('Table_LME_GeneticsCognition_%s.csv',COG_VAR)));
summaryTable=resultsTable(interactionTerms,{'GeneFactor','Model','Term','Beta','SE','CI_lo','CI_hi','pValue','pFDR'});
summaryTable.pFDR=round(summaryTable.pFDR,4);
writetable(summaryTable, fullfile(outputDir,sprintf('Table_LME_InteractionsSummary_%s.csv',COG_VAR)));

demoVars={'age_at_visit','SEX','EDUCYRS','apathy_composite',COG_VAR,'APOE_e4','updrs3_score','depression_factor'};
demoNames={'Age','Sex(%F)','Education','Apathy',sprintf('Cognition(%s)',COG_VAR),'APOE(%)','UPDRS-III','Depression'};
demoCell={};
for gi=1:nGG
    idx=bl_g.gene_subgroup==validGenes(gi); n=sum(idx);
    row={char(validGenes(gi)),n};
    for vi=1:numel(demoVars)
        if ~ismember(demoVars{vi},bl_g.Properties.VariableNames), row{end+1}='N/A'; continue; end
        vals=double(safeCol(bl_g,demoVars{vi})); vals=vals(idx); vals=vals(~isnan(vals));
        if strcmpi(demoVars{vi},'sex')||strcmp(demoVars{vi},'APOE_e4'), row{end+1}=sprintf('%.1f',100*mean(vals>=1));
        else, row{end+1}=sprintf('%.2f +/- %.2f',mean(vals),std(vals)); end
    end
    demoCell(end+1,:)=row;
end
writetable(cell2table(demoCell,'VariableNames',['Subgroup','N',demoNames]), fullfile(outputDir,'Table_Genetic_Demographics.csv'));
fprintf('    All tables saved.\n');


%% ====  CONSOLE SUMMARY  =================================================
fprintf('\n%s\n',repmat('=',1,72));
fprintf('  DV         : %s (z-scored -> cog_z)\n', COG_VAR);
fprintf('  Apathy col : %s    ICB col: %s\n', apathy_col_lme, ternary(hasICB,icb_col_lme,'(unavailable)'));
fprintf('  Models per gene factor: %s\n', strjoin(modelTypes, ', '));
fprintf('  Covariates : %s\n', covStr);
fprintf('  Error bars : within-subject (Loftus-Masson 1994, Morey 2008)\n\n');
fprintf('  %-12s  %-12s  %-32s  %8s  %8s  %8s\n','Gene','Model','Term','Beta_z','p_unc','p_FDR');
fprintf('  %s\n',repmat('-',1,82));
for ri=1:height(summaryTable)
    sig='';
    if ~isnan(summaryTable.pFDR(ri))
        if summaryTable.pFDR(ri)<0.001, sig=' ***';
        elseif summaryTable.pFDR(ri)<0.01, sig=' **';
        elseif summaryTable.pFDR(ri)<0.05, sig=' *'; end
    end
    fprintf('  %-12s  %-12s  %-32s  %8.4f  %8.4f  %8.4f%s\n',...
        summaryTable.GeneFactor{ri},summaryTable.Model{ri},summaryTable.Term{ri},...
        summaryTable.Beta(ri),summaryTable.pValue(ri),...
        ternary(isnan(summaryTable.pFDR(ri)),NaN,summaryTable.pFDR(ri)),sig);
end
fprintf('\n');


%% =========================================================================
%  LOCAL HELPER FUNCTIONS
%  =========================================================================

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

function [lme, converged] = tryFitLME(tbl, f1, f2)
    converged = true;
    try, lme = fitlme(tbl, f1, 'FitMethod', 'REML');
    catch, try, lme = fitlme(tbl, f2, 'FitMethod', 'REML'); converged = false;
           catch, lme = []; converged = false; end
    end
end

function [glme, converged] = tryFitGLME(tbl, f1, f2, distName, linkName)
% Poisson/log fallback companion to tryFitLME for count outcomes (e.g. QUIP).
    converged = true;
    try
        glme = fitglme(tbl, f1, 'Distribution', distName, 'Link', linkName);
    catch
        try
            glme = fitglme(tbl, f2, 'Distribution', distName, 'Link', linkName);
            converged = false;
        catch
            glme = []; converged = false;
        end
    end
end

function rc = appendLMEResults(rc, lme, gn, modelType, tbl)
    c = lme.Coefficients; ns = numel(unique(tbl.PATNO));
    no = height(tbl); a = lme.ModelCriterion.AIC; ci = coefCI(lme);
    for ri = 1:height(c)
        rc(end+1,:) = {gn, modelType, c.Name{ri}, c.Estimate(ri), c.SE(ri), c.tStat(ri), ...
            c.DF(ri), c.pValue(ri), ci(ri,1), ci(ri,2), no, ns, a}; %#ok<AGROW>
    end
end

function [f1, f2] = buildModelFormula(modelType, includeAPOE, hasCovariates, covStr)
%BUILDMODELFORMULA  Construct the two LME formulas (random-slope, random-intercept).
    switch modelType
        case 'apathy_only'
            core = 'cog_z ~ year_num*gene_flag*apathy_flag';
        case 'icb_only'
            core = 'cog_z ~ year_num*gene_flag*icb_flag';
        case 'combined'
            core = 'cog_z ~ year_num*gene_flag*apathy_flag + year_num*gene_flag*icb_flag';
        otherwise
            error('buildModelFormula: unknown modelType "%s"', modelType);
    end
    if includeAPOE,    core = [core ' + apoe4_bin']; end
    if hasCovariates,  core = [core ' + ' covStr];  end
    f1 = [core ' + (year_num|PATNO)'];
    f2 = [core ' + (1|PATNO)'];
end

function pFDR = fdr_bh(pV)
    m=numel(pV); [ps,si]=sort(pV); pa=nan(m,1); pa(m)=ps(m);
    for i=m-1:-1:1, pa(i)=min(pa(i+1),ps(i)*m/i); end
    pa=min(pa,1); pFDR=nan(m,1); pFDR(si)=pa;
end

function [lo,hi] = wilsonCI(p,n)
    z=1.96; d=1+z^2/n; c=(p+z^2/(2*n))/d;
    m=z*sqrt((p*(1-p)+z^2/(4*n))/n)/d;
    lo=max(0,c-m); hi=min(1,c+m);
end

function yS = gaussSmooth(x, y, sigma)
    yS = nan(size(y));
    for i = 1:numel(x)
        w = exp(-0.5 * ((x - x(i)) / sigma).^2);
        w = w / sum(w);
        yS(i) = sum(w .* y);
    end
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





function out = ternary(cond, a, b)
    if cond, out = a; else, out = b; end
end

function exportFigure(fig, filepath)
    set(fig,'PaperPositionMode','auto','Renderer','painters');
    print(fig, [filepath '.pdf'], '-dpdf', '-bestfit');
    print(fig, [filepath '.png'], '-dpng', '-r600');
    close(fig);
end

function formatAxis(ax, COL)
    set(ax,'FontSize',7,'Box','off','TickDir','out','XGrid','off','YGrid','off', ...
        'GridColor',COL.grid,'GridAlpha',0.4);
end

function s = geneDisplayName(g)
%GENEDISPLAYNAME figure label for a subgroup or gene factor:
%gene symbols in italics (TeX), GBA -> GBA1, APOE_e4 -> APOE e4 (Greek
%epsilon). Display only; the data strings are never changed.
    g = char(g);
    switch g
        case 'Sporadic',   s = 'Sporadic';
        case 'GBA',        s = '\itGBA1\rm';
        case 'LRRK2',      s = '\itLRRK2\rm';
        case 'SNCA',       s = '\itSNCA\rm';
        case 'PRKN/PINK1', s = '\itPRKN\rm/\itPINK1\rm';
        case 'LRRK2+GBA',  s = '\itLRRK2\rm+\itGBA1\rm';
        case 'APOE_e4',    s = ['\itAPOE\rm ' char(949) '4'];
        otherwise,         s = strrep(g, '_', '\_');
    end
end

function s = fmtP(p)
%FMTP p-value text for figure annotations: '< 0.001' below
%0.001, otherwise '= 0.xxx'.
    if isnan(p),        s = '= n/a';
    elseif p < 0.001,   s = '< 0.001';
    else,               s = sprintf('= %.3f', p);
    end
end

function auditFigureText(fig, scale, csvPath)
%AUDITFIGURETEXT legibility QA: lists every visible text
%element (tick labels, text objects incl. titles/labels, legends, colorbar
%ticks) with its printed size = FontSize x scale, where scale = inclusion
%width / canvas width. TeX sub/superscripts are drawn at 0.8x the base
%size, reported separately. Writes a CSV and prints the minima.
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

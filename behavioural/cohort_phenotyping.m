% Cohort phenotyping for the PPMI analyses of apathy and impulse-control behaviours (ICB) in
%     Parkinson's disease. Reads the curated PPMI data cut and the GDS-15 and MoCA item files,
%     builds the apathy composite (MDS-UPDRS Part I item 1.5 plus GDS-15 items 2, 9 and 13, range
%     0 to 7) and the depression factor (GDS-15 without the apathy items and item 10), defines the
%     four apathy by ICB groups (apathy composite of 2 or more; QUIP-Current-Short with at least
%     one positive domain) and the baseline PD cohort, reports the internal consistency of the
%     apathy composite and builds the long-format tables with the bifactor g scores merged.
% Inputs are read from $PPMI_CLINICAL_ROOT/ppmi_data (the curated cut, and the item CSVs under
%     download/) and from cfa_factor_scores.csv, written by behavioural/cognition_bifactor_cfa.R
%     and copied into outputDir ($PPMI_CLINICAL_ROOT/revision/figures/regenerated) before running.
% The workspace it leaves (T, isPD, bl, groupLabels, groupLabelsShort, THRESH, COL, outputDir,
%     Tpd_lme, lmeT, okBase and coeffToTable) is the input of
%     behavioural/figure1_cohort_characterisation.m, figure2_cognitive_trajectories.m and
%     figure5_genetic_moderation.m and of imaging/connectivity_outcomes/fc_cognitive_slope_mediation.m.
% Manuscript items are the cohort and four-group counts (Results, Figure 1a) and the reliability
%     of the apathy composite (console output).

clear; clc; close all;

%% ========================  0. CONFIGURATION  ============================

% --- Paths ---------------------------------------------------------------
basePath    = fullfile(getenv('PPMI_CLINICAL_ROOT'), 'ppmi_data');
dlDir     = fullfile(basePath, 'download');
curatedFile = fullfile(basePath, 'PPMI_Curated_Data_Cut_Public_20260223.xlsx');
sheetName   = '20260105';

% Raw instrument CSV paths (only needed for item-level data)
downloadPath = fullfile(basePath, 'download');
gdsFile   = fullfile(downloadPath, 'Geriatric_Depression_Scale__Short_Version__16Mar2026.csv');
mocaFile  = fullfile(downloadPath, 'Montreal_Cognitive_Assessment__MoCA__16Mar2026.csv');

% --- Output directory ----------------------------------------------------
outputDir = fullfile(getenv('PPMI_CLINICAL_ROOT'), 'revision/figures/regenerated');
if ~exist(outputDir, 'dir'), mkdir(outputDir); end

% --- Clinical thresholds -------------------------------------------------
THRESH.apathy_comp = 2;  % Apathy composite (NP1APAT + 3 GDS items, range 0-7) >= 2
THRESH.quip_icb    = 1;  % QUIP >= 1 positive domain

% --- Publication colour palette (Okabe-Ito, colourblind-safe) ------------
COL.neither  = [0.40 0.76 0.65];   % teal
COL.apathy   = [0.90 0.56 0.13];   % amber
COL.icb      = [0.34 0.47 0.73];   % steel blue
COL.overlap  = [0.80 0.33 0.38];   % rose
COL.groups   = [COL.neither; COL.apathy; COL.icb; COL.overlap];
COL.bg       = [1 1 1];            % white
COL.grid     = [0.85 0.85 0.85];
COL.txt      = [0.15 0.15 0.15];
COL.accent   = [0.10 0.10 0.10];

groupLabels      = {'Neither','Apathy only','ICB only','Apathy + ICB'};
groupLabelsShort = {'Neither','Apathy','ICB','Both'};
% --- Typography defaults -------------------------------------------------
setPublicationDefaults();

fprintf('==========================================================\n');
fprintf('  PPMI: Apathy & Impulsivity in Parkinson''s Disease\n');
fprintf('==========================================================\n\n');

%% ====================  1. DATA LOADING  =================================
fprintf('[1] Loading curated dataset + raw instrument CSVs...\n');

% --- 1A. Curated dataset (primary: demographics, scores, imaging, genetics)
T = readtable(curatedFile, 'Sheet', sheetName, 'VariableNamingRule', 'preserve');
T.PATNO = string(T.PATNO);
fprintf('    Curated: %d rows x %d cols\n', height(T), width(T));

% Curated file column name reference (verified for 20260223 cut):
%   Demographics:  PATNO, COHORT, subgroup, analytic_subgroup, SEX, age_at_visit, EDUCYRS, YEAR, EVENT_ID
%   Motor:         updrs3_score, updrs2_score, updrs1_score, NHY, hy, MSEADLG, LEDD
%   Apathy item:   NP1APAT
%   Cognitive:     moca, bjlot, lns, SDMTOTAL, VLTANIM, TMT_A, TMT_B, hvlt_immediaterecall,
%                  hvlt_discrimination, hvlt_retention, DVT_TOTAL_RECALL, DVT_RECOG_DISC_INDEX
%   Neurobehav:    quip, quip_any, ess, rem, gds, stai_trait, stai_state, scopa
%   Imaging:       MIA_CAUDATE_BILAT, MIA_PUTAMEN_BILAT, MIA_CAUDATE_R/L, MIA_PUTAMEN_R/L
%   Genetics:      APOE, APOE_e4, subgroup
%   Staging:       NSD_STAGE, Stage_S, Stage_D, Stage_G, cogstate
%   New:           Clinical_Stage, COG_COMPOSITE_EXT, COG_COMPOSITE_INT

% --- 1B. GDS raw items (for apathy composite & depression factor) --------
if isfile(gdsFile)
    gdsRaw = readtable(gdsFile, 'VariableNamingRule', 'preserve');
    gdsRaw.PATNO = string(gdsRaw.PATNO);
    fprintf('    GDS raw items: %d rows\n', height(gdsRaw));
else
    gdsRaw = [];
    fprintf('    WARNING: GDS raw CSV not found\n');
end

% --- 1C. MoCA raw items (for subdomain analysis) ------------------------
if isfile(mocaFile)
    mocaRaw = readtable(mocaFile, 'VariableNamingRule', 'preserve');
    mocaRaw.PATNO = string(mocaRaw.PATNO);
    fprintf('    MoCA raw items: %d rows\n', height(mocaRaw));
else
    mocaRaw = [];
    fprintf('    WARNING: MoCA raw CSV not found\n');
end
fprintf('\n');

%% ====================  2. COMPOSITE SCORES  =============================
fprintf('[2] Computing composite scores from raw items...\n');

% --- 2A. GDS Apathy items & Depression factor ----------------------------
%  PPMI GDS coding: 0 = No, 1 = Yes (raw response).
%  Reverse-scored items (positive wording, Yes = NOT depressed):
%    Item 1 (GDSSATIS), 5 (GDSGSPIR), 7 (GDSHAPPY), 11 (GDSALIVE), 13 (GDSENRGY)
%  Direct-scored items (negative wording, Yes = depressed): all others
%
%  GDS item number → PPMI column:
%    1-GDSSATIS  2-GDSDROPD  3-GDSEMPTY  4-GDSBORED   5-GDSGSPIR
%    6-GDSAFRAD  7-GDSHAPPY  8-GDSHLPLS  9-GDSHOME   10-GDSMEMRY
%   11-GDSALIVE 12-GDSWRTLS 13-GDSENRGY 14-GDSHOPLS  15-GDSBETER

if ~isempty(gdsRaw)
    gdsItemCols = {'GDSSATIS','GDSDROPD','GDSEMPTY','GDSBORED','GDSGSPIR', ...
                   'GDSAFRAD','GDSHAPPY','GDSHLPLS','GDSHOME','GDSMEMRY', ...
                   'GDSALIVE','GDSWRTLS','GDSENRGY','GDSHOPLS','GDSBETER'};
    reverseIdx = [1, 5, 7, 11, 13];  % Items where Yes=0 (not depressed)

    % Score all 15 items in the depression direction
    gdsScored = nan(height(gdsRaw), 15);
    for it = 1:15
        if ismember(gdsItemCols{it}, gdsRaw.Properties.VariableNames)
            rawV = safeNumeric(gdsRaw.(gdsItemCols{it}));
            if ismember(it, reverseIdx)
                gdsScored(:, it) = 1 - rawV;   % Reverse: Yes(1)→0, No(0)→1
            else
                gdsScored(:, it) = rawV;        % Direct: Yes(1)→1
            end
        end
    end

    % Apathy 3 GDS items: items 2 (GDSDROPD), 9 (GDSHOME), 13 (GDSENRGY reversed)
    ap3 = gdsScored(:, [2, 9, 13]);
    gdsRaw.gds_apathy_3item = nanSumIfComplete(ap3, 2);

    % Depression factor: items 1,3,4,5,6,7,8,11,14,15 (excl. apathy 2,9,13 & memory 10)
    depIdx = [1, 3, 4, 5, 6, 7, 8, 11, 14, 15];
    gdsRaw.depression_factor = nanSumIfComplete(gdsScored(:, depIdx), 7);

    % Merge into main table
    gds_merge = gdsRaw(:, intersect(gdsRaw.Properties.VariableNames, ...
        {'PATNO','EVENT_ID','gds_apathy_3item','depression_factor'}, 'stable'));
    [~, uIdx] = unique(gds_merge(:, {'PATNO','EVENT_ID'}), 'rows', 'first');
    gds_merge = gds_merge(uIdx, :);
    T = outerjoin(T, gds_merge, 'Keys', {'PATNO','EVENT_ID'}, 'MergeKeys', true, 'Type', 'left');
    fprintf('    GDS: apathy 3-item + depression factor (10-item, reversed) merged\n');
end

% --- 2B. MoCA subdomains from raw items ---------------------------------
%  Standard MoCA subdomain groupings with CORRECT item assignments:
%    Visuospatial/Executive (5 pts): MCAALTTM + MCACUBE + MCACLCKC + MCACLCKN + MCACLCKH
%    Naming (3 pts):                 MCALION + MCARHINO + MCACAMEL
%    Attention (6 pts):              MCAFDS + MCABDS + MCAVIGIL + MCASER7(0-3)
%    Language (3 pts):               MCASNTNC(0-2) + MCAVF(0-1)  ← NOT MCAVFNUM!
%    Abstraction (2 pts):            MCAABSTR(0-2)
%    Delayed Recall (5 pts):         MCAREC1 + MCAREC2 + MCAREC3 + MCAREC4 + MCAREC5
%    Orientation (6 pts):            MCADATE + MCAMONTH + MCAYR + MCADAY + MCAPLACE + MCACITY

if ~isempty(mocaRaw)
    mc = mocaRaw;
    mc.moca_visuospatial = safeNumeric(mc.MCAALTTM) + safeNumeric(mc.MCACUBE) + ...
                           safeNumeric(mc.MCACLCKC) + safeNumeric(mc.MCACLCKN) + safeNumeric(mc.MCACLCKH);
    mc.moca_naming       = safeNumeric(mc.MCALION) + safeNumeric(mc.MCARHINO) + safeNumeric(mc.MCACAMEL);
    mc.moca_attention    = safeNumeric(mc.MCAFDS) + safeNumeric(mc.MCABDS) + ...
                           safeNumeric(mc.MCAVIGIL) + safeNumeric(mc.MCASER7);
    mc.moca_language     = safeNumeric(mc.MCASNTNC) + safeNumeric(mc.MCAVF);  % MCAVF not MCAVFNUM!
    mc.moca_abstraction  = safeNumeric(mc.MCAABSTR);
    mc.moca_memory       = safeNumeric(mc.MCAREC1) + safeNumeric(mc.MCAREC2) + ...
                           safeNumeric(mc.MCAREC3) + safeNumeric(mc.MCAREC4) + safeNumeric(mc.MCAREC5);
    mc.moca_orientation  = safeNumeric(mc.MCADATE) + safeNumeric(mc.MCAMONTH) + ...
                           safeNumeric(mc.MCAYR) + safeNumeric(mc.MCADAY) + ...
                           safeNumeric(mc.MCAPLACE) + safeNumeric(mc.MCACITY);

    moca_merge = mc(:, intersect(mc.Properties.VariableNames, ...
        {'PATNO','EVENT_ID','moca_visuospatial','moca_naming','moca_attention', ...
         'moca_language','moca_abstraction','moca_memory','moca_orientation'}, 'stable'));
    [~, uIdx] = unique(moca_merge(:, {'PATNO','EVENT_ID'}), 'rows', 'first');
    moca_merge = moca_merge(uIdx, :);
    T = outerjoin(T, moca_merge, 'Keys', {'PATNO','EVENT_ID'}, 'MergeKeys', true, 'Type', 'left');
    fprintf('    MoCA subdomains: 7 domains computed (language uses MCAVF, not MCAVFNUM)\n');
end

% --- 2C. Apathy/Anhedonia Composite (NP1APAT + 3 GDS items) -------------
np1apat = safeCol(T, 'NP1APAT');
gds_ap3 = safeCol(T, 'gds_apathy_3item');
T.apathy_composite = np1apat + gds_ap3;
T.apathy_composite(isnan(np1apat) | isnan(gds_ap3)) = NaN;
fprintf('    Apathy composite (0-7): %d valid rows\n', sum(~isnan(T.apathy_composite)));
fprintf('\n');

%% ============  3. DEFINE APATHY × ICB GROUPS  ==========================
fprintf('[3] Defining clinical groups...\n');

apathyFlag = safeCol(T, 'apathy_composite') >= THRESH.apathy_comp;
quipVals   = safeCol(T, 'quip');  % From curated file (already computed)
icbFlag    = quipVals >= THRESH.quip_icb;

T.group4 = nan(height(T), 1);
canGroup = ~isnan(safeCol(T, 'apathy_composite')) & ~isnan(quipVals);
T.group4(canGroup & ~apathyFlag & ~icbFlag) = 1;  % Neither
T.group4(canGroup &  apathyFlag & ~icbFlag) = 2;  % Apathy only
T.group4(canGroup & ~apathyFlag &  icbFlag) = 3;  % ICB only
T.group4(canGroup &  apathyFlag &  icbFlag) = 4;  % Both

for g = 1:4
    fprintf('    %15s : %5d observations\n', groupLabels{g}, sum(T.group4 == g, 'omitnan'));
end
fprintf('    %15s : %5d\n', 'Ungroupable', sum(~canGroup));

%% ============  4. BASELINE PD COHORT  ===================================
fprintf('\n[4] Selecting baseline PD cohort...\n');

isPD = safeCol(T, 'COHORT') == 1;  % COHORT=1 is PD in PPMI
isBL = safeCol(T, 'YEAR') == 0;
bl = T(isBL & isPD & ~isnan(T.group4), :);
fprintf('    Baseline PD with valid groups: N = %d\n', height(bl));
for g = 1:4
    fprintf('    %15s : %4d\n', groupLabels{g}, sum(bl.group4 == g));
end
fprintf('\n');

%% ============  5. INTERNAL CONSISTENCY  =================================

fprintf('[5] Internal consistency of apathy composite...\n');

if ~isempty(gdsRaw)
    % --- Merge individual GDS apathy items from gdsRaw into bl -----------
    gds_items = gdsRaw(:, intersect(gdsRaw.Properties.VariableNames, ...
        {'PATNO','EVENT_ID','GDSDROPD','GDSHOME','GDSENRGY'}, 'stable'));
    [~, uIdx] = unique(gds_items(:, {'PATNO','EVENT_ID'}), 'rows', 'first');
    gds_items = gds_items(uIdx, :);
    bl = outerjoin(bl, gds_items, 'Keys', {'PATNO','EVENT_ID'}, ...
                   'MergeKeys', true, 'Type', 'left');

    % --- Build item matrix (all scored in the apathy direction) ----------
    apathyItemNames = {'NP1APAT','GDSDROPD','GDSHOME','GDSENRGY'};
    apathyItemMat = nan(height(bl), 4);

    apathyItemMat(:, 1) = safeCol(bl, 'NP1APAT');            % 0-4, higher = more apathetic
    apathyItemMat(:, 2) = safeNumeric(bl.GDSDROPD);          % Yes(1) = apathetic
    apathyItemMat(:, 3) = safeNumeric(bl.GDSHOME);           % Yes(1) = apathetic
    apathyItemMat(:, 4) = 1 - safeNumeric(bl.GDSENRGY);     % Yes(1) = energetic → reverse

    validAlpha = all(~isnan(apathyItemMat), 2);
    if sum(validAlpha) > 30
        alphaVal = cronbach_alpha(apathyItemMat(validAlpha, :));
        fprintf('    Cronbach''s alpha = %.3f (n = %d)\n', alphaVal, sum(validAlpha));
        fprintf('    Item-total correlations (corrected):\n');
        for ai = 1:4
            restScore = sum(apathyItemMat(validAlpha, setdiff(1:4, ai)), 2);
            r = corr(apathyItemMat(validAlpha, ai), restScore, 'Type', 'Spearman');
            fprintf('      %-10s : r = %.3f\n', apathyItemNames{ai}, r);
        end
    else
        fprintf('    Insufficient complete cases (n = %d)\n', sum(validAlpha));
    end
end
fprintf('\n');

%% ========================================================================
%  23. LONG-FORMAT TABLES FOR THE LINEAR MIXED-EFFECTS MODELS
%  Tpd_lme (all PD visits) and lmeT (one row per visit with the baseline covariates). Apathy
%  and ICB are defined at baseline (time-invariant, z-scored); the models themselves are
%  fitted by figure2_cognitive_trajectories.m and figure5_genetic_moderation.m.
% =========================================================================
fprintf('[23] Building the long-format tables for the mixed-effects models...\n');

% --- Prepare long-format analysis table ----------------------------------
Tpd_lme = T(isPD, :);

% Baseline continuous scores (carry forward to all visits per patient)
blFlags = Tpd_lme(safeCol(Tpd_lme,'YEAR')==0, :);
blFlagSub = table(blFlags.PATNO, ...
    safeCol(blFlags,'apathy_composite'), ...       % continuous apathy
    safeCol(blFlags,'quip'), ...                    % continuous ICB
    safeCol(blFlags,'depression_factor'), ...
    safeCol(blFlags,'SEX'), ...
    safeCol(blFlags,'EDUCYRS'), ...
    safeCol(blFlags,'duration_yrs'), ...
    safeCol(blFlags,'LEDD'), ...
    safeCol(blFlags,'NHY'), ...
    safeCol(blFlags,'age'), ...
    safeCol(blFlags,'updrs3_score'), ...
    'VariableNames',{'PATNO','bl_apathy_raw','bl_icb_raw','bl_depfactor', ...
                     'bl_sex','bl_educ','bl_duration','bl_ledd','bl_hy','bl_age','bl_updrs3_score'});
% Remove duplicate PATNOs (take first)
[~,uI] = unique(blFlagSub.PATNO,'first');
blFlagSub = blFlagSub(uI,:);

% Z-score apathy and ICB for interpretability (1-unit = 1 SD)
apMu = mean(blFlagSub.bl_apathy_raw,'omitnan');
apSD = std(blFlagSub.bl_apathy_raw,0,'omitnan');
icbMu = mean(blFlagSub.bl_icb_raw,'omitnan');
icbSD = std(blFlagSub.bl_icb_raw,0,'omitnan');
blFlagSub.bl_apathy_z = (blFlagSub.bl_apathy_raw - apMu) ./ apSD;
blFlagSub.bl_icb_z    = (blFlagSub.bl_icb_raw    - icbMu) ./ icbSD;

fprintf('    Apathy composite: M=%.2f, SD=%.2f (range %d-%d)\n', ...
    apMu, apSD, min(blFlagSub.bl_apathy_raw), max(blFlagSub.bl_apathy_raw));
fprintf('    QUIP total:       M=%.2f, SD=%.2f (range %d-%d)\n', ...
    icbMu, icbSD, min(blFlagSub.bl_icb_raw), max(blFlagSub.bl_icb_raw));

% Merge baseline flags into full longitudinal table
Tpd_lme = outerjoin(Tpd_lme, blFlagSub, 'Keys','PATNO','MergeKeys',true,'Type','left');

% Time variable (continuous years)
Tpd_lme.time = safeCol(Tpd_lme, 'YEAR');

% MoCA
Tpd_lme.moca_val = safeCol(Tpd_lme, 'moca');

% Compute number of missing yearly visits per patient (out of years 0-5)
patIDs = unique(Tpd_lme.PATNO);
nMissVec = nan(numel(patIDs), 1);
for pi = 1:numel(patIDs)
    pIdx = Tpd_lme.PATNO == patIDs(pi);
    yrsPresent = unique(safeCol(Tpd_lme(pIdx,:),'YEAR'));
    yrsPresent = yrsPresent(~isnan(yrsPresent) & yrsPresent >= 0 & yrsPresent <= 5);
    nMissVec(pi) = 6 - numel(yrsPresent);
end
Tpd_lme.bl_nmiss = nan(height(Tpd_lme), 1);
for pi = 1:numel(patIDs)
    pIdx = Tpd_lme.PATNO == patIDs(pi);
    Tpd_lme.bl_nmiss(pIdx) = nMissVec(pi);
end

% --- Build fitlme-compatible table ---------------------------------------
lmeT = table( ...
    Tpd_lme.PATNO, ...
    Tpd_lme.time, ...
    Tpd_lme.moca_val, ...
    Tpd_lme.bl_apathy_z, ...
    Tpd_lme.bl_icb_z, ...
    Tpd_lme.bl_depfactor, ...
    Tpd_lme.bl_sex, ...
    Tpd_lme.bl_educ, ...
    Tpd_lme.bl_nmiss, ...
    Tpd_lme.bl_duration, ...
    Tpd_lme.bl_ledd, ...
    Tpd_lme.bl_hy, ...
    Tpd_lme.bl_age, ...
    Tpd_lme.bl_updrs3_score, ...
    'VariableNames', {'PATNO','time','moca','apathy','icb', ...
                      'depfactor','sex','educ','nmiss','duration','ledd','hy','age','updrs3_score'});

% Sex remains categorical
lmeT.sex = categorical(lmeT.sex);

% Valid rows
okBase = ~isnan(lmeT.apathy) & ~isnan(lmeT.icb) & ~isnan(lmeT.time);

% Helper: coefficient table extraction
coeffToTable = @(mdl) lme_coeff_to_table(mdl);

%% ========================================================================
%  23b. MERGE OF THE BIFACTOR g SCORES
%  Run after behavioural/cognition_bifactor_cfa.R has exported cfa_factor_scores.csv
% =========================================================================
fprintf('[23b] Importing CFA factor scores...\n');

% --- Import CFA factor scores from R output --------------------------------
cfaPath = fullfile(outputDir, 'cfa_factor_scores.csv');  % adjust path as needed
cfaScores = readtable(cfaPath);
cfaScores.PATNO = string(cfaScores.PATNO);
Tpd_lme = T(isPD, :);
fprintf('    Loaded CFA factor scores: %d rows\n', height(cfaScores));

% Standardise column names (R export may use dots instead of underscores)
cfaScores.Properties.VariableNames = strrep(cfaScores.Properties.VariableNames, '.', '_');
Tpd_lme = outerjoin(Tpd_lme, cfaScores, ...
    'LeftKeys',{'PATNO','YEAR'}, 'RightKeys',{'PATNO','YEAR'}, ...
    'MergeKeys',true, 'Type','left');

% Clean up merge column
if ismember('YEAR_merge', Tpd_lme.Properties.VariableNames)
    Tpd_lme.YEAR_merge = [];
end
fprintf('    Factor scores merged. Available: bifactor_g\n');
fprintf('    Non-missing bifactor_g: %d / %d rows\n', ...
    sum(~isnan(Tpd_lme.bifactor_g)), height(Tpd_lme));

% ---- Key join of the factor scores. lmeT rows follow the curated-cut order
% (PATNO, EVENT_ID) while Tpd_lme is sorted by the outerjoin keys (PATNO, YEAR); the two
% orders differ for unscheduled visits, so every copy below uses w219_loc, the Tpd_lme row
% with the same PATNO and YEAR as each lmeT row, instead of the row position.
w219_kL = string(lmeT.PATNO) + "_" + string(double(lmeT.time));
w219_kR = string(Tpd_lme.PATNO) + "_" + string(double(Tpd_lme.YEAR));
assert(numel(unique(w219_kR)) == numel(w219_kR), 'duplicate PATNO+YEAR keys in Tpd_lme');
[w219_tf, w219_loc] = ismember(w219_kL, w219_kR);
assert(all(w219_tf), '%d lmeT rows without a PATNO+YEAR match', sum(~w219_tf));
% --- Rebuild LME table with factor scores ----------------------------------
lmeT.bifactor_g = Tpd_lme.bifactor_g(w219_loc);


% --- Local functions -----------------------------------------------------
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

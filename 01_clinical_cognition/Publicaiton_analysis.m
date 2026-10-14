% Role: Master clinical script: cohort, apathy-composite and QUIP-S phenotyping, the four-group
%     classification and the longitudinal tables with the bifactor g scores merged (lmeT, Tpd_lme,
%     TpdSub),. Its workspace is the input of 06_figures/Figure1.m, Figure2.m and Figure3.m and of
%     05_mediation_moderation/Figure4.m. It reads cfa_factor_scores.csv (written by
%     02_cfa_bifactor/ppm_cfa_Gemeni_bothModels.R) from its own output folder outputDir
%     ($PPMI_CLINICAL_ROOT/revision/figures/regenerated): copy the file there before running.
% Manuscript item: cohort and four-group counts (Results; Figure 1a); master workspace of Figures 1,
%     2, 4 and 5 and of Table 1; apathy-composite reliability (console output).
%% ========================================================================
%  PPMI Analysis: Apathy/Anhedonia & Impulse-Control Behaviours in PD
%  Multi-modal Investigation of Co-existing and Distinctive Features
% =========================================================================
%  Author : Bahaa Atallah
%  Date   : March 2026
%  Data   : PPMI Curated Data Cut (Public, 2024-12-11) + Raw CSV Instruments
%
%  Apathy/anhedonia composite (validated 4-item measure):
%     MDS-UPDRS Part I item 1.5 (NP1APAT, range 0-4)
%     + GDS-15 item 2  (GDSDROPD, "Dropped activities?", Yes=1)
%     + GDS-15 item 9  (GDSHOME,  "Prefer to stay home?", Yes=1)
%     + GDS-15 item 13 (GDSENRGY, "Full of energy?", REVERSED: No→1)
%     Range: 0-7
%     References: O'Callaghan & Lewis (2021); Pluck & Brown (2002)
%
%  Depression factor (GDS-15 excl. apathy items 2,9,13 and memory item 10):
%     Items 1,3,4,5,6,7,8,11,14,15 with reverse scoring on items 1,5,7,11
%     Range: 0-11
%
%  ICB: QUIP-Current-Short ≥ 1 positive domain
%
%  Genetics: APOE ε4 carrier status; genetic subgroup (Sporadic/LRRK2/GBA)
% =========================================================================

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
volFile     = fullfile(dlDir,  'FS7_ASEG_VOL_16Mar2026.csv');

% --- Output directory ----------------------------------------------------
% REVISION DRIVER COPY (2026-08-29): identical to GD/Publicaiton_analysis.m
% except this one line: ALL outputs are redirected into the revision
% workspace so no existing generated file is overwritten.
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
nGrp       = 4;
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

%% ========  6. TABLE 1: BASELINE DEMOGRAPHICS  ==========================
fprintf('[6] Table 1: Baseline demographics...\n');

t1Defs = {
    'age_at_visit',       'Age (years)',             'cont';
    'SEX',                'Sex (% female)',          'binary_0f';
    'EDUCYRS',            'Education (years)',       'cont';
    'duration_yrs',       'Disease duration (yrs)',  'cont';
    'moca',               'MoCA total (0-30)',       'cont';
    'apathy_composite',   'Apathy composite (0-7)',  'cont';
    'depression_factor',  'Depression (10-item GDS)', 'cont';
    'gds',                'GDS total (0-15)',        'cont';
    'quip',               'QUIP-S ICB domains',      'cont';
    'stai_trait',         'STAI trait anxiety',       'cont';
    'updrs3_score',       'MDS-UPDRS III motor',     'cont';
    'updrs2_score',       'MDS-UPDRS II ADL',        'cont';
    'MSEADLG',            'Schwab-England ADL (%)',  'cont';
    'NHY',                'Hoehn & Yahr stage',      'cont';
    'ess',                'ESS sleepiness (0-24)',   'cont';
    'rem',                'RBDSQ total',             'cont';
    'scopa',              'SCOPA-AUT total',         'cont';
    'MIA_CAUDATE_BILAT',       'DaT caudate (mean)',      'cont';
    'MIA_PUTAMEN_BILAT',       'DaT putamen (mean)',      'cont';
    'LEDD',               'Levodopa equiv. dose',    'cont';
};

summaryRows = {};
for v = 1:size(t1Defs, 1)
    vname = t1Defs{v,1}; vlabel = t1Defs{v,2}; vtype = t1Defs{v,3};
    row = {vlabel};
    for g = 1:4
        vals = safeCol(bl, vname);
        vals = vals(bl.group4 == g);
        vals = vals(~isnan(vals));
        n = numel(vals);
        if n == 0
            row{end+1} = '---'; %#ok
        elseif strcmp(vtype, 'binary_0f')
            row{end+1} = sprintf('%.1f%% (%d)', 100*sum(vals==0)/n, n); %#ok
        else
            row{end+1} = sprintf('%.2f +/- %.2f', mean(vals), std(vals)); %#ok
        end
    end
    % Kruskal-Wallis
    allV = safeCol(bl, vname); grpV = bl.group4;
    ok = ~isnan(allV) & ~isnan(grpV);
    if sum(ok) > 10
        p = kruskalwallis(allV(ok), grpV(ok), 'off');
        row{end+1} = ternary(p < 0.001, '< .001', sprintf('%.3f', p)); %#ok
    else
        row{end+1} = '---'; %#ok
    end
    summaryRows = [summaryRows; row]; %#ok
end

header = ['Variable', groupLabelsShort, 'p (KW)'];
summaryT = cell2table(summaryRows, 'VariableNames', header);
disp(summaryT);
% writetable(summaryT, fullfile(outputDir, 'Table1_Baseline_Demographics.csv'));
% fprintf('    -> Table1_Baseline_Demographics.csv\n\n');

%% =======  7. FIG 1: PREVALENCE & CO-OCCURRENCE  ========================
fprintf('[7] Figure 1: Prevalence...\n');

nBL = height(bl);
counts = arrayfun(@(g) sum(bl.group4 == g), 1:4);
pcts   = 100 * counts / nBL;

fig1 = figure('Units','centimeters','Position',[2 2 17.4 7],'Color',COL.bg);

% (a) Donut chart
ax1 = axes('Position',[0.02 0.08 0.42 0.82]); hold on;
theta0 = pi/2; cumP = [0, cumsum(pcts)]/100;
for g = 1:4
    t = linspace(theta0-cumP(g)*2*pi, theta0-cumP(g+1)*2*pi, 100);
    fill([cos(t), 0.55*cos(fliplr(t))], [sin(t), 0.55*sin(fliplr(t))], ...
        COL.groups(g,:), 'EdgeColor','w','LineWidth',1.5);
end
text(0,0, {sprintf('N = %d', nBL),'Baseline PD'}, 'HorizontalAlignment','center', ...
    'FontSize',8,'FontWeight','bold','Color',COL.accent);
axis equal off;
text(-1.15,1.15,'(a)','FontSize',10,'FontWeight','bold','Color',COL.accent);

% (b) Horizontal bar
ax2 = axes('Position',[0.52 0.15 0.44 0.72]);
bh = barh(1:4, pcts, 0.55, 'FaceColor','flat','EdgeColor','none');
bh.CData = COL.groups;
set(ax2,'YTick',1:4,'YTickLabel',groupLabels,'YDir','reverse','FontSize',8, ...
    'Box','off','TickDir','out','XGrid','on','GridColor',COL.grid,'GridAlpha',0.5);
xlabel('Prevalence (%)','FontSize',9);
for g = 1:4
    text(pcts(g)+1, g, sprintf('%.1f%%  (n=%d)',pcts(g),counts(g)), ...
        'FontSize',7.5,'Color',COL.txt,'VerticalAlignment','middle');
end
xlim([0 max(pcts)*1.4]);
title('Prevalence at Baseline','FontSize',10,'FontWeight','bold','Color',COL.accent);
text(-0.15,-0.05,'(b)','Units','normalized','FontSize',10,'FontWeight','bold','Color',COL.accent);

exportFigure(fig1, fullfile(outputDir, 'Fig1_Prevalence'));
fprintf('    -> Fig1_Prevalence\n\n');

%%
%% ======= MAIN SCRIPT: CLINICAL PROFILE RADAR =================
fprintf('[Process] Preparing data for publication radar plot...\n');

% 1. Define dimensions and logic
radarDefs = {
    'updrs3_score',     'Motor',           false;
    'moca',             'Cognition',       true;  % Invert (Higher = better)
    'apathy_composite', 'Apathy',          false;
    'depression_factor','Depression',      false;
    'stai_trait',       'Anxiety',         false;
    'quip',             'Impulsivity',     false;
    'ess',              'Sleepiness',      false;
    'rem',              'REM Sleep',       false;
    'scopa',            'Autonomic',       false;
    'MSEADLG',          'ADL Function',    true;  % Invert (Higher = better)
};

nAxes = size(radarDefs, 1);
axisLabels = radarDefs(:, 2)';
isInverted = [radarDefs{:, 3}];
varNames = radarDefs(:, 1);

% 2. Extract and Normalize Data (0 to 1 scale)
allData = bl{:, varNames};
lo = prctile(allData, 1, 1);
hi = prctile(allData, 99, 1);
ranges = hi - lo;
ranges(ranges == 0) = 1; % Prevent divide-by-zero

nGroups = 4;
plotData = nan(nGroups, nAxes);

for g = 1:nGroups
    idx = (bl.group4 == g);
    if ~any(idx), continue; end
    
    % Get raw mean for the group
    m = nanmean(allData(idx, :), 1);
    
    % Normalize to 0-1
    val = (m - lo) ./ ranges;
    
    % Invert specific dimensions so 0=Best, 1=Worst across the board
    val(isInverted) = 1 - val(isInverted);
    
    % Clamp to strictly 0 to 1
    plotData(g, :) = max(0, min(1, val));
end

% 3. Plotting Setup
% Define colors if not already in your workspace
colors = COL.groups(1:4, :); 
groupNames = groupLabels(1:4);

% Create figure and call the plotting function
figR = figure('Units', 'centimeters', 'Position', [2 2 16 16], 'Color', COL.bg);
clf(figR);

% Call the modular function
radarPlotPub(plotData, axisLabels, groupNames, colors, COL);

% Add Custom Titles
text(gca, 0, 1.35, 'Clinical Profile by Group', ...
    'FontSize', 12, 'FontWeight', 'bold', 'HorizontalAlignment', 'center', 'Color', COL.accent);
text(gca, 0, 1.25, 'Normalised means (0 = Best, 1 = Worst)', ...
    'FontSize', 8, 'HorizontalAlignment', 'center', 'Color', [0.5 0.5 0.5]);

% exportFigure(figR, fullfile(outputDir, 'Fig_Radar_ClinicalProfile'));

%% ========================================================================
%  24. SANKEY DIAGRAM & MOTIVATIONAL INSTABILITY PHENOTYPING
% =========================================================================
fprintf('[24] Sankey diagram: group transitions over time...\n');
Tpd = T(isPD & ~isnan(T.group4), :);

% --- 1. BUILD TRANSITION MATRIX ------------------------------------------
sankeyYears = [0, 2, 4];  
sankeyPatients = unique(Tpd.PATNO);
nSP = numel(sankeyPatients);
patGroupMat = nan(nSP, numel(sankeyYears));
yrCol = safeCol(Tpd, 'YEAR');

for ti = 1:numel(sankeyYears)
    for pi = 1:nSP
        idx = Tpd.PATNO == sankeyPatients(pi) & yrCol == sankeyYears(ti);
        if any(idx)
            grps = Tpd.group4(idx);
            grps = grps(~isnan(grps));
            if ~isempty(grps)
                patGroupMat(pi, ti) = grps(1);
            end
        end
    end
end

% Keep only patients with data at all time points
complete = all(~isnan(patGroupMat), 2);
patGroupComplete = patGroupMat(complete, :);
completePatIDs = sankeyPatients(complete); % Need this to map back to baseline later

fprintf('    Patients with group data at years %s: %d / %d\n', ...
    mat2str(sankeyYears), sum(complete), nSP);

%% --- 2. PUBLICATION-QUALITY SANKEY DIAGRAM ------------------------------
if sum(complete) > 30
    fig17 = figure('Units','centimeters','Position',[2 2 18 11],'Color',COL.bg);
    ax = axes('Color',COL.bg); hold on;
    nTrans = numel(sankeyYears) - 1;  
    nG = 4;
    xPositions = 1:numel(sankeyYears);
    
    % Sankey layout parameters
    gap = 0.04; % Gap between group blocks
    scale = 1 - (gap * (nG-1)); % Remaining space for data
    bw = 0.12; % Bar width
    
    % Draw blocks and flows
    nodeY = cell(numel(sankeyYears), 1);
    
    for ti = 1:numel(sankeyYears)
        grpCounts = arrayfun(@(g) sum(patGroupComplete(:,ti)==g), 1:nG);
        totalN = sum(grpCounts);
        
        % Calculate Y positions with gaps
        h = (grpCounts / totalN) * scale;
        starts = zeros(1, nG);
        currentY = 0;
        for g = 1:nG
            starts(g) = currentY;
            currentY = currentY + h(g) + gap;
        end
        nodeY{ti} = [starts; starts + h]; % [bottoms; tops]
        
        % Draw stacked bars
        for g = 1:nG
            yBot = nodeY{ti}(1, g);
            yTop = nodeY{ti}(2, g);
            if yTop <= yBot, continue; end
            
            rectangle('Position', [xPositions(ti)-bw, yBot, 2*bw, yTop-yBot], ...
                'FaceColor', COL.groups(g,:), 'EdgeColor', 'w', 'LineWidth', 0.8, ...
                'Curvature', 0.05);
            
            if (yTop-yBot) > 0.04
                text(xPositions(ti), (yBot+yTop)/2, sprintf('%s\nn=%d', groupLabelsShort{g}, grpCounts(g)), ...
                    'HorizontalAlignment','center','FontSize',6,'FontWeight','bold', 'Color','w');
            end
        end
    end
    
    % Draw flow bands
    for ti = 1:nTrans
        grpCounts1 = arrayfun(@(g) sum(patGroupComplete(:,ti)==g), 1:nG);
        offL = zeros(nG,1); offR = zeros(nG,1);
        
        for g1 = 1:nG
            for g2 = 1:nG
                nTr = sum(patGroupComplete(:,ti)==g1 & patGroupComplete(:,ti+1)==g2);
                if nTr == 0, continue; end
                
                hFlow = (nTr / sum(grpCounts1)) * scale;
                yL_bot = nodeY{ti}(1, g1) + offL(g1);
                yL_top = yL_bot + hFlow;
                offL(g1) = offL(g1) + hFlow;
                
                yR_bot = nodeY{ti+1}(1, g2) + offR(g2);
                yR_top = yR_bot + hFlow;
                offR(g2) = offR(g2) + hFlow;
                
                % Smooth S-curve
                xL = xPositions(ti) + bw;
                xR = xPositions(ti+1) - bw;
                tt = linspace(0, 1, 50);
                xCurve = xL + (xR-xL) * (3*tt.^2 - 2*tt.^3);
                yBotCurve = yL_bot + (yR_bot-yL_bot) * (3*tt.^2 - 2*tt.^3);
                yTopCurve = yL_top + (yR_top-yL_top) * (3*tt.^2 - 2*tt.^3);
                
                % Color blending
                if g1 == g2
                    bandCol = COL.groups(g1,:); alpha = 0.25;
                else
                    bandCol = (COL.groups(g1,:) + COL.groups(g2,:)) / 2; alpha = 0.15;
                end
                
                fill([xCurve, fliplr(xCurve)], [yBotCurve, fliplr(yTopCurve)], ...
                    bandCol, 'FaceAlpha', alpha, 'EdgeColor', 'none');
            end
        end
    end
    
    set(ax,'XTick',xPositions,'XTickLabel',arrayfun(@(y)sprintf('Year %d',y),sankeyYears,'UniformOutput',false), ...
        'FontSize',9,'Box','off','YTick',[],'TickDir','out');
    ylim([-0.05 1.05]); xlim([0.5 numel(sankeyYears)+0.5]);
    title('Group Transitions Over Time', 'FontSize',11,'FontWeight','bold','Color',COL.accent);
    % exportFigure(fig17, fullfile(outputDir, 'Fig_Sankey_Transitions'));
end

%% --- 3. IDENTIFY MOTIVATIONAL STABILITY & DETAILED SUBGROUPS ------------
fprintf('\n[25] Profiling Motivational Instability Subgroups...\n');

% Initialize tracking arrays
isStable   = false(size(patGroupComplete, 1), 1);
isUnstable = false(size(patGroupComplete, 1), 1);
detGroup   = nan(size(patGroupComplete, 1), 1);

% Logic Definitions for Detailed Subgroups:
% 1: Persistent Neither (Always 1)
% 2: Persistent Apathy  (Always 2)
% 3: Persistent ICB     (Always 3)
% 4: Apathy -> New ICB  (Starts 2, moves to 3 or 4)
% 5: ICB -> New Apathy  (Starts 3, moves to 2 or 4)

for i = 1:size(patGroupComplete,1)
    path = patGroupComplete(i, :);
    
    % --- Main 2-Group Logic ---
    if all(path == path(1))
        isStable(i) = true;
    end
    
    for t = 1:(size(patGroupComplete,2)-1)
        if ismember(path(t), [2,3,4]) && ismember(path(t+1), [2,3,4]) && path(t) ~= path(t+1)
            isUnstable(i) = true;
        end
    end
    
    % --- Granular 5-Group Logic ---
    if all(path == 1), detGroup(i) = 1; end % Persist Neither
    if all(path == 2), detGroup(i) = 2; end % Persist Apathy
    if all(path == 3), detGroup(i) = 3; end % Persist ICB
    
    % Apathy developing new ICB
    if path(1) == 2 && any(ismember(path(2:end), [3, 4]))
        detGroup(i) = 4;
    end
    
    % ICB developing new Apathy
    if path(1) == 3 && any(ismember(path(2:end), [2, 4]))
        detGroup(i) = 5;
    end
end

% Map back to baseline table
bl.StabilityGroup = nan(height(bl), 1);
bl.StabilityGroup(ismember(bl.PATNO, completePatIDs(isStable))) = 1;
bl.StabilityGroup(ismember(bl.PATNO, completePatIDs(isUnstable))) = 2;

bl.DetGroup = nan(height(bl), 1);
for i = 1:length(completePatIDs)
    bl.DetGroup(bl.PATNO == completePatIDs(i)) = detGroup(i);
end
%%


%%
%% --- 4. RADAR & STRIATAL PLOTS (BASELINE & LONGITUDINAL TRAJECTORIES) ---
% Define valid longitudinal subgroups
idx_valid = (bl.StabilityGroup == 1 | bl.StabilityGroup == 2);

% -------------------------------------------------------------------------
% (A) GLOBAL COVARIATE ADJUSTMENT (Calculated on the full baseline cohort)
% -------------------------------------------------------------------------
c_dat = safeCol(bl, 'MIA_CAUDATE_BILAT');
p_dat = safeCol(bl, 'MIA_PUTAMEN_BILAT');
age_v = safeCol(bl, 'age'); 

sex_raw = safeCol(bl, 'SEX');
if iscell(sex_raw) || iscategorical(sex_raw) || isstring(sex_raw)
    sex_v = double(categorical(sex_raw));
else
    sex_v = sex_raw; 
end

validScan = ~isnan(c_dat) & ~isnan(p_dat);
bl.PLOT_STRIATAL = nan(height(bl), 1);
if any(validScan)
    bl.PLOT_STRIATAL(validScan) = mean([c_dat(validScan), p_dat(validScan)], 2); 
end

validAll = validScan & ~isnan(age_v) & ~isnan(sex_v);
isAdjusted = false;
if sum(validAll) > 5 
    X_cov = [age_v(validAll), sex_v(validAll)];
    if size(unique(X_cov, 'rows'), 1) > 2
        mdl = fitlm(X_cov, bl.PLOT_STRIATAL(validAll));
        bl.PLOT_STRIATAL(validAll) = mdl.Residuals.Raw + nanmean(bl.PLOT_STRIATAL(validAll));
        isAdjusted = true;
    end
end

% Now extract the longitudinal subset containing the adjusted data
blSub = bl(idx_valid, :);

% -------------------------------------------------------------------------
% (B) FIGURE GENERATION (1x5 Layout)
% -------------------------------------------------------------------------
if height(blSub) > 5
    figSub = figure('Units','centimeters','Position',[2 2 38 8.5],'Color',COL.bg);
    
    subColors = [0.4 0.6 0.4; 0.8 0.3 0.3]; % Green for Stable, Red for Unstable
    pathColors = [0.4 0.4 0.7; 0.8 0.3 0.3]; % Purple for Persistent Path, Red for Unstable
    detColors = [0.7 0.7 0.7; 0.3 0.5 0.7; 0.8 0.5 0.3; 0.6 0.4 0.7; 0.7 0.3 0.4]; 
    
    % --- PLOT 1: RADAR PLOT (Main 2 Groups) ------------------------------
    axR = subplot(1,5,1); hold on;
    radarDefs = {
        'updrs3_score',      'Motor',           false;
        'moca',              'Cognition',       true;  
        'apathy_composite',  'Apathy',          false;
        'quip',              'Impulsivity',     false;
        'MIA_CAUDATE_BILAT', 'Caudate DaT',     true;  
        'MIA_PUTAMEN_BILAT', 'Putamen DaT',     true;  
    };
    
    varNames = radarDefs(:,1); axisLabels = radarDefs(:,2)'; isInverted = [radarDefs{:,3}];
    nVars = length(varNames); lo = nan(1, nVars); hi = nan(1, nVars);
    
    for v = 1:nVars
        colD = safeCol(bl, varNames{v}); valD = colD(~isnan(colD));
        if ~isempty(valD), lo(v) = prctile(valD, 1); hi(v) = prctile(valD, 99); end
    end
    ranges = hi - lo; ranges(ranges == 0) = 1; 
    
    radarMeans = nan(2, nVars);
    for g = 1:2
        for v = 1:nVars
            colD = safeCol(blSub(blSub.StabilityGroup == g, :), varNames{v});
            valD = colD(~isnan(colD));
            if ~isempty(valD)
                val = (mean(valD) - lo(v)) / ranges(v);
                if isInverted(v), val = 1 - val; end
                radarMeans(g, v) = max(0, min(1, val));
            end
        end
    end
    
    axes(axR); 
    radarPlotPub(radarMeans, axisLabels, {'Stable', 'Unstable'}, subColors, COL);
    title('(a) Baseline Profiles (0=Best, 1=Worst)', 'FontSize', 9, 'Color', COL.accent);

    % --- PLOT 2: CROSS-SECTIONAL BASELINE (4 Groups) ---------------------
    axS1 = subplot(1,5,2); hold on;
    m_base = nan(1,4); se_base = nan(1,4); vals_base = cell(1,4);
    
    for g = 1:4
        v = bl.PLOT_STRIATAL(bl.group4 == g); 
        vals_base{g} = v(~isnan(v));
        if ~isempty(vals_base{g})
            m_base(g) = mean(vals_base{g});
            se_base(g) = std(vals_base{g})/sqrt(length(vals_base{g}));
            bar(g, m_base(g), 0.6, 'FaceColor', COL.groups(g,:), 'EdgeColor', 'none', 'FaceAlpha', 0.85);
        end
    end
    errorbar(1:4, m_base, se_base, 'k.', 'LineWidth', 1, 'CapSize', 3);
    
    % Stats generator for multiple bars
    y_bracket = max(m_base + se_base) + 0.1;
    for i=1:3
        for j=i+1:4
            if length(vals_base{i})>1 && length(vals_base{j})>1
                [~, p_val] = ttest2(vals_base{i}, vals_base{j});
                if p_val < 0.05
                    plot([i, i, j, j], [y_bracket, y_bracket+0.03, y_bracket+0.03, y_bracket], 'Color', [0.3 0.3 0.3]);
                    if p_val < 0.001, p_str = '***'; elseif p_val < 0.01, p_str = '**'; else, p_str = '*'; end
                    text((i+j)/2, y_bracket+0.06, p_str, 'HorizontalAlignment','center','FontSize',9, 'Color', COL.txt);
                    y_bracket = y_bracket + 0.12;
                end
            end
        end
    end
    ylim([min(0, min(m_base-se_base)-0.2), y_bracket]);
    set(axS1, 'XTick', 1:4, 'XTickLabel', groupLabelsShort, 'Box', 'off'); xtickangle(25); 
    if isAdjusted, ylabel('Adjusted DaT Binding'); else, ylabel('Raw DaT Binding'); end
    title('(b) Baseline Cross-Sectional', 'FontSize', 9, 'Color', COL.accent);

    % --- PLOT 3: MAIN LONGITUDINAL (Stable vs Unstable) ------------------
    axS2 = subplot(1,5,3); hold on;
    val_g1 = blSub.PLOT_STRIATAL(blSub.StabilityGroup == 1); val_g1 = val_g1(~isnan(val_g1));
    val_g2 = blSub.PLOT_STRIATAL(blSub.StabilityGroup == 2); val_g2 = val_g2(~isnan(val_g2));
    
    m_main = [mean(val_g1), mean(val_g2)]; se_main = [std(val_g1)/sqrt(length(val_g1)), std(val_g2)/sqrt(length(val_g2))];
    for g=1:2
        if ~isnan(m_main(g)), bar(g, m_main(g), 0.6, 'FaceColor', subColors(g,:), 'EdgeColor', 'none', 'FaceAlpha', 0.85); end
    end
    errorbar(1:2, m_main, se_main, 'k.', 'LineWidth', 1, 'CapSize', 4);
    
    y_bracket = max(m_main + se_main) + 0.1;
    if length(val_g1) > 1 && length(val_g2) > 1
        [~, p_val] = ttest2(val_g1, val_g2);
        if p_val < 0.05
            plot([1, 1, 2, 2], [y_bracket, y_bracket+0.03, y_bracket+0.03, y_bracket], 'Color', [0.3 0.3 0.3]);
            if p_val < 0.001, p_str = '***'; elseif p_val < 0.01, p_str = '**'; else, p_str = '*'; end
            text(1.5, y_bracket+0.06, p_str, 'HorizontalAlignment','center','FontSize',10, 'Color', COL.txt);
            y_bracket = y_bracket + 0.12;
        end
    end
    ylim([min(0, min(m_main-se_main)-0.2), y_bracket]);
    set(axS2, 'XTick', 1:2, 'XTickLabel', {'All Stable', 'All Unstable'}, 'Box', 'off'); xtickangle(15); 
    title('(c) Longitudinal Stability', 'FontSize', 9, 'Color', COL.accent);

    % --- PLOT 4: PERSISTENT PATHOLOGY vs UNSTABLE ------------------------
    axS3 = subplot(1,5,4); hold on;
    val_pers = blSub.PLOT_STRIATAL(blSub.DetGroup == 2 | blSub.DetGroup == 3); val_pers = val_pers(~isnan(val_pers));
    m_path = [mean(val_pers), mean(val_g2)]; se_path = [std(val_pers)/sqrt(length(val_pers)), std(val_g2)/sqrt(length(val_g2))];
    
    for g = 1:2
        if ~isnan(m_path(g)), bar(g, m_path(g), 0.6, 'FaceColor', pathColors(g,:), 'EdgeColor', 'none', 'FaceAlpha', 0.85); end
    end
    errorbar(1:2, m_path, se_path, 'k.', 'LineWidth', 1, 'CapSize', 4);
    
    y_bracket = max(m_path + se_path) + 0.1;
    if length(val_pers) > 1 && length(val_g2) > 1
        [~, p_val] = ttest2(val_pers, val_g2);
        if p_val < 0.05
            plot([1, 1, 2, 2], [y_bracket, y_bracket+0.03, y_bracket+0.03, y_bracket], 'Color', [0.3 0.3 0.3]);
            if p_val < 0.001, p_str = '***'; elseif p_val < 0.01, p_str = '**'; else, p_str = '*'; end
            text(1.5, y_bracket+0.06, p_str, 'HorizontalAlignment','center','FontSize',10, 'Color', COL.txt);
            y_bracket = y_bracket + 0.12;
        end
    end
    ylim([min(0, min(m_path-se_path)-0.2), y_bracket]);
    set(axS3, 'XTick', 1:2, 'XTickLabel', {'Persist Path.', 'Unstable'}, 'Box', 'off'); xtickangle(15);
    title('(d) Static vs Shifting Path.', 'FontSize', 9, 'Color', COL.accent);

    % --- PLOT 5: GRANULAR TRAJECTORIES (5 Groups) ------------------------
    axS4 = subplot(1,5,5); hold on;
    detNames = {'Persist Neither', 'Persist Apathy', 'Persist ICB', 'Apathy \rightarrow ICB', 'ICB \rightarrow Apathy'};
    
    m_det = nan(1, 5); se_det = nan(1, 5); vals_det = cell(1,5);
    for g = 1:5
        v = blSub.PLOT_STRIATAL(blSub.DetGroup == g); 
        vals_det{g} = v(~isnan(v));
        if ~isempty(vals_det{g})
            m_det(g) = mean(vals_det{g});
            se_det(g) = std(vals_det{g}) / sqrt(length(vals_det{g}));
            bar(g, m_det(g), 0.6, 'FaceColor', detColors(g,:), 'EdgeColor', 'none', 'FaceAlpha', 0.85);
            text(g, min(0, m_det(g)-se_det(g)-0.05), sprintf('n=%d', length(vals_det{g})), 'HorizontalAlignment','center','FontSize',6,'Color',[0.4 0.4 0.4]);
        end
    end
    errorbar(1:5, m_det, se_det, 'k.', 'LineWidth', 1, 'CapSize', 3);
    
    y_bracket = max(m_det + se_det) + 0.1;
    for i=1:4
        for j=i+1:5
            if length(vals_det{i})>1 && length(vals_det{j})>1
                [~, p_val] = ttest2(vals_det{i}, vals_det{j});
                if p_val < 0.05
                    plot([i, i, j, j], [y_bracket, y_bracket+0.03, y_bracket+0.03, y_bracket], 'Color', [0.3 0.3 0.3]);
                    if p_val < 0.001, p_str = '***'; elseif p_val < 0.01, p_str = '**'; else, p_str = '*'; end
                    text((i+j)/2, y_bracket+0.06, p_str, 'HorizontalAlignment','center','FontSize',9, 'Color', COL.txt);
                    y_bracket = y_bracket + 0.12;
                end
            end
        end
    end
    ylim([min(0, min(m_det-se_det)-0.2), y_bracket]);
    set(axS4, 'XTick', 1:5, 'XTickLabel', detNames, 'Box', 'off'); xtickangle(30);
    title('(e) Granular Trajectories', 'FontSize', 9, 'Color', COL.accent);
    
    sgtitle('Phenotyping Motivational Instability: Cross-Sectional to Longitudinal Trajectories', 'FontSize', 12, 'FontWeight', 'bold', 'Color', COL.accent);
    % exportFigure(figSub, fullfile(outputDir, 'Fig_Motivational_Trajectories_Full'));
end

%%
%% ===  10. FIG 4: LONGITUDINAL TRAJECTORIES (ALL SUBGROUP DIVISIONS) ======
fprintf('[10] Figure 4: Longitudinal trajectories across all 4 subgroup divisions...\n');

years = 0:5;
longVars  = {'moca','updrs3_score','apathy_composite','depression_factor','quip'};
longNames = {'MoCA (Cognition)','UPDRS-III (Motor)','Apathy Composite', ...
             'Depression (GDS)','QUIP (ICB)'};
longYlab  = {'MoCA Score (0-30)','UPDRS-III Score','Score (0-7)', ...
             'Score (0-10)','ICB Domains'};

% 1. Map both StabilityGroup and DetGroup from 'bl' to the longitudinal 'Tpd' table
Tpd.StabilityGroup = nan(height(Tpd), 1);
Tpd.DetGroup = nan(height(Tpd), 1);
for i = 1:height(bl)
    idx = (Tpd.PATNO == bl.PATNO(i));
    Tpd.StabilityGroup(idx) = bl.StabilityGroup(i);
    Tpd.DetGroup(idx) = bl.DetGroup(i);
end

% 2. Filter safely
TpdSub = Tpd(~isnan(Tpd.StabilityGroup), :);

% Colors from previous plots
subColors = [0.4 0.6 0.4; 0.8 0.3 0.3]; 
pathColors = [0.4 0.4 0.7; 0.8 0.3 0.3]; 
detColors = [0.7 0.7 0.7; 0.3 0.5 0.7; 0.8 0.5 0.3; 0.6 0.4 0.7; 0.7 0.3 0.4]; 

% Create a massive 4-row, 5-column figure
fig4 = figure('Units','centimeters','Position',[2 2 32 24],'Color',COL.bg);
nRows = 4;

for row = 1:nRows
    
    % Setup logic for the current row
    if row == 1
        nGrps = 4; cMap = COL.groups(1:4,:); 
        rowTitle = 'Original Baseline Groups'; lgNames = groupLabelsShort;
    elseif row == 2
        nGrps = 2; cMap = subColors; 
        rowTitle = 'Stable vs Unstable'; lgNames = {'Stable', 'Unstable'};
    elseif row == 3
        nGrps = 2; cMap = pathColors; 
        rowTitle = 'Persistent Path. vs Unstable'; lgNames = {'Persist Path.', 'Unstable'};
    elseif row == 4
        nGrps = 5; cMap = detColors; 
        rowTitle = 'Granular Trajectories'; lgNames = {'Persist Neither', 'Persist Apathy', 'Persist ICB', 'Apathy \rightarrow ICB', 'ICB \rightarrow Apathy'};
    end
    
    for lv = 1:numel(longVars)
        ax = subplot(nRows, numel(longVars), (row-1)*numel(longVars) + lv);
        set(ax,'Color',COL.bg); hold on;
        
        ph4 = gobjects(nGrps,1); 
        all_mY = nan(nGrps, numel(years)); 
        all_seY = nan(nGrps, numel(years)); 
        
        % --- Calculate and Plot Trajectories ---
        for g = 1:nGrps
            mY = nan(1, numel(years)); semY = nan(1, numel(years)); nY = nan(1, numel(years));
            
            for yi = 1:numel(years)
                % Define group indexing based on the row logic
                if row == 1
                    idx = (TpdSub.group4 == g);
                elseif row == 2
                    idx = (TpdSub.StabilityGroup == g);
                elseif row == 3
                    if g == 1, idx = (TpdSub.DetGroup == 2 | TpdSub.DetGroup == 3); end
                    if g == 2, idx = (TpdSub.StabilityGroup == 2); end
                elseif row == 4
                    idx = (TpdSub.DetGroup == g);
                end
                
                yr_idx = idx & (TpdSub.YEAR == years(yi));
                v = safeCol(TpdSub(yr_idx, :), longVars{lv});
                v = v(~isnan(v));
                
                nY(yi) = numel(v);
                if nY(yi) > 0, mY(yi) = mean(v); semY(yi) = std(v)/sqrt(nY(yi)); end
            end
            
            vYrs = nY > 0;
            all_mY(g,:) = mY; all_seY(g,:) = semY;
            ph4(g) = plotTrajectory(ax, years, mY, semY, nY, vYrs, cMap(g,:), 'endpoints');
        end
        
        % --- Add Year-by-Year Statistical Significance ---
        for yi = 1:numel(years)
            % Re-extract data for testing
            grpData = cell(1, nGrps);
            for g = 1:nGrps
                if row == 1, idx = (TpdSub.group4 == g);
                elseif row == 2, idx = (TpdSub.StabilityGroup == g);
                elseif row == 3
                    if g==1, idx = (TpdSub.DetGroup == 2 | TpdSub.DetGroup == 3); else, idx = (TpdSub.StabilityGroup == 2); end
                elseif row == 4, idx = (TpdSub.DetGroup == g);
                end
                v = safeCol(TpdSub(idx & (TpdSub.YEAR == years(yi)), :), longVars{lv});
                grpData{g} = v(~isnan(v));
            end
            
            y_loc = nanmax(all_mY(:,yi) + all_seY(:,yi)) + 0.05 * abs(nanmax(all_mY(:,yi)));
            
            % 2-Group t-tests (Rows 2 and 3)
            if nGrps == 2 && length(grpData{1}) > 2 && length(grpData{2}) > 2
                [~, p_val] = ttest2(grpData{1}, grpData{2});
                if p_val < 0.05
                    if p_val < 0.001, p_str = '***'; elseif p_val < 0.01, p_str = '**'; else, p_str = '*'; end
                    text(years(yi), y_loc, p_str, 'HorizontalAlignment', 'center', 'FontSize', 12, 'Color', COL.txt, 'FontWeight', 'bold');
                end
            end
            
            % Multi-Group ANOVA (Rows 1 and 4)
            if nGrps > 2
                % Combine for anova1
                allVals = []; allGrps = [];
                for g = 1:nGrps
                    if length(grpData{g}) > 1
                        allVals = [allVals; grpData{g}]; %#ok
                        allGrps = [allGrps; repmat(g, length(grpData{g}), 1)]; %#ok
                    end
                end
                if length(unique(allGrps)) > 1 % Need at least 2 groups with data to run ANOVA
                    p_val = anova1(allVals, allGrps, 'off');
                    if p_val < 0.05
                        text(years(yi), y_loc, sprintf('p=%.3f', p_val), 'HorizontalAlignment', 'center', 'FontSize', 7, 'Color', COL.accent, 'FontWeight', 'bold');
                    end
                end
            end
        end
        
        % --- Formatting ---
        xlabel('Year','FontSize',8); 
        if lv == 1, ylabel(rowTitle, 'FontSize', 10, 'FontWeight', 'bold', 'Color', COL.accent); end
        if row == 1, title(longNames{lv}, 'FontSize', 10, 'FontWeight', 'bold', 'Color', COL.accent); end
        
        set(ax,'FontSize',7,'Box','off','XGrid','on','YGrid','on', 'GridColor',COL.grid,'GridAlpha',0.4,'TickDir','out');
        xlim([-0.3 max(years)+0.3]);
        
        % Expand Y-limits slightly to fit asterisks/p-values
        yl = ylim; ylim([yl(1), yl(2) + diff(yl)*0.18]);
        
        % Add Legend to the first column of each row
        if lv == 1
            valid_ph = ph4(isgraphics(ph4));
            valid_lg = lgNames(isgraphics(ph4));
            lg = legend(valid_ph, valid_lg, 'Location','best','FontSize',6,'Box','off');
            lg.ItemTokenSize = [12 12];
        end
        hold off;
    end
end

sgtitle('Longitudinal Trajectories Across All Clinical Groupings', 'FontSize', 14, 'FontWeight', 'bold', 'Color', COL.accent);
exportFigure(fig4, fullfile(outputDir, 'Fig4_Longitudinal_MegaGrid_Complete'));
fprintf('    -> Fig4_Longitudinal_MegaGrid_Complete\n\n');

%%
%% ===  11. FIG 5: LONGITUDINAL COGNITIVE DOMAIN TRAJECTORIES  ===========
fprintf('[11] Figure 5: Longitudinal cognitive domain composites...\n');

% -------------------------------------------------------------------------
% 1. CREATE COGNITIVE COMPOSITES (Z-SCORED)
% -------------------------------------------------------------------------
% Define raw variables and their clinical direction
% Direction: 1 = Higher is better, -1 = Lower is better (e.g., time in seconds)
rawVars    = {'DVT_TOTAL_RECALL', 'DVT_RECOG_DISC_INDEX', 'lns', 'VLTANIM', 'SDMTOTAL', 'TMT_B', 'bjlot'};
directions = [ 1,                  1,                      1,     1,         1,         -1,      1];

% Pre-allocate Z-score columns in Tpd
for i = 1:numel(rawVars)
    v = rawVars{i};
    % Calculate Baseline Mean and SD to anchor the Z-scores
    bl_vals = safeCol(Tpd(Tpd.YEAR == 0, :), v);
    bl_vals = bl_vals(~isnan(bl_vals));
    mu = mean(bl_vals);
    sig = std(bl_vals);
    
    % Z-score the entire longitudinal column and apply direction inversion
    raw_col = safeCol(Tpd, v);
    z_col = directions(i) * ((raw_col - mu) / sig);
    Tpd.(['Z_', v]) = z_col;
end

% Group Z-scores into Domains using row-wise nanmean
Tpd.comp_memory = nanmean([Tpd.Z_DVT_TOTAL_RECALL, Tpd.Z_DVT_RECOG_DISC_INDEX], 2);
Tpd.comp_exec   = nanmean([Tpd.Z_lns, Tpd.Z_VLTANIM, Tpd.Z_SDMTOTAL, Tpd.Z_TMT_B], 2);
Tpd.comp_visuo  = Tpd.Z_bjlot;

% -------------------------------------------------------------------------
% 2. FIGURE PREPARATION
% -------------------------------------------------------------------------
years = 0:5;
cogVars  = {'moca', 'comp_memory', 'comp_exec', 'comp_visuo'};
cogNames = {'Global Cognition (MoCA)', 'Memory Composite', 'Executive & Attention', 'Visuospatial Function'};
cogYLab  = {'Raw Score (0-30)', 'Z-Score', 'Z-Score', 'Z-Score'};

% Filter safely
TpdSub = Tpd(~isnan(Tpd.StabilityGroup), :);

% Colors
subColors = [0.4 0.6 0.4; 0.8 0.3 0.3]; 
pathColors = [0.4 0.4 0.7; 0.8 0.3 0.3]; 
detColors = [0.7 0.7 0.7; 0.3 0.5 0.7; 0.8 0.5 0.3; 0.6 0.4 0.7; 0.7 0.3 0.4]; 

% Create a massive 4-row, 4-column figure
fig5 = figure('Units','centimeters','Position',[2 2 28 24],'Color',COL.bg);
nRows = 4;

% -------------------------------------------------------------------------
% 3. MEGA-GRID PLOTTING LOOP
% -------------------------------------------------------------------------
for row = 1:nRows
    
    % Setup logic for the current row
    if row == 1
        nGrps = 4; cMap = COL.groups(1:4,:); 
        rowTitle = 'Original Baseline Groups'; lgNames = groupLabelsShort;
    elseif row == 2
        nGrps = 2; cMap = subColors; 
        rowTitle = 'Stable vs Unstable'; lgNames = {'Stable', 'Unstable'};
    elseif row == 3
        nGrps = 2; cMap = pathColors; 
        rowTitle = 'Persistent Path. vs Unstable'; lgNames = {'Persist Path.', 'Unstable'};
    elseif row == 4
        nGrps = 5; cMap = detColors; 
        rowTitle = 'Granular Trajectories'; lgNames = {'Persist Neither', 'Persist Apathy', 'Persist ICB', 'Apathy \rightarrow ICB', 'ICB \rightarrow Apathy'};
    end
    
    for cv = 1:numel(cogVars)
        ax = subplot(nRows, numel(cogVars), (row-1)*numel(cogVars) + cv);
        set(ax,'Color',COL.bg); hold on;
        
        ph5 = gobjects(nGrps,1); 
        all_mY = nan(nGrps, numel(years)); 
        all_seY = nan(nGrps, numel(years)); 
        
        % --- Calculate and Plot Trajectories ---
        for g = 1:nGrps
            mY = nan(1, numel(years)); semY = nan(1, numel(years)); nY = nan(1, numel(years));
            
            for yi = 1:numel(years)
                % Group indexing logic
                if row == 1, idx = (TpdSub.group4 == g);
                elseif row == 2, idx = (TpdSub.StabilityGroup == g);
                elseif row == 3
                    if g == 1, idx = (TpdSub.DetGroup == 2 | TpdSub.DetGroup == 3); end
                    if g == 2, idx = (TpdSub.StabilityGroup == 2); end
                elseif row == 4, idx = (TpdSub.DetGroup == g);
                end
                
                yr_idx = idx & (TpdSub.YEAR == years(yi));
                v = safeCol(TpdSub(yr_idx, :), cogVars{cv});
                v = v(~isnan(v));
                
                nY(yi) = numel(v);
                if nY(yi) > 0, mY(yi) = mean(v); semY(yi) = std(v)/sqrt(nY(yi)); end
            end
            
            vYrs = nY > 0;
            all_mY(g,:) = mY; all_seY(g,:) = semY;
            ph5(g) = plotTrajectory(ax, years, mY, semY, nY, vYrs, cMap(g,:), 'endpoints');
        end
        
        % --- Add Year-by-Year Statistical Significance ---
        for yi = 1:numel(years)
            grpData = cell(1, nGrps);
            for g = 1:nGrps
                if row == 1, idx = (TpdSub.group4 == g);
                elseif row == 2, idx = (TpdSub.StabilityGroup == g);
                elseif row == 3
                    if g==1, idx = (TpdSub.DetGroup == 2 | TpdSub.DetGroup == 3); else, idx = (TpdSub.StabilityGroup == 2); end
                elseif row == 4, idx = (TpdSub.DetGroup == g);
                end
                v = safeCol(TpdSub(idx & (TpdSub.YEAR == years(yi)), :), cogVars{cv});
                grpData{g} = v(~isnan(v));
            end
            
            y_loc = nanmax(all_mY(:,yi) + all_seY(:,yi)) + 0.05 * abs(nanmax(all_mY(:,yi)));
            
            % 2-Group t-tests
            if nGrps == 2 && length(grpData{1}) > 2 && length(grpData{2}) > 2
                [~, p_val] = ttest2(grpData{1}, grpData{2});
                if p_val < 0.05
                    if p_val < 0.001, p_str = '***'; elseif p_val < 0.01, p_str = '**'; else, p_str = '*'; end
                    text(years(yi), y_loc, p_str, 'HorizontalAlignment', 'center', 'FontSize', 12, 'Color', COL.txt, 'FontWeight', 'bold');
                end
            end
            
            % Multi-Group ANOVA
            if nGrps > 2
                allVals = []; allGrps = [];
                for g = 1:nGrps
                    if length(grpData{g}) > 1
                        allVals = [allVals; grpData{g}]; %#ok
                        allGrps = [allGrps; repmat(g, length(grpData{g}), 1)]; %#ok
                    end
                end
                if length(unique(allGrps)) > 1 
                    p_val = anova1(allVals, allGrps, 'off');
                    if p_val < 0.05
                        text(years(yi), y_loc, sprintf('p=%.3f', p_val), 'HorizontalAlignment', 'center', 'FontSize', 7, 'Color', COL.accent, 'FontWeight', 'bold');
                    end
                end
            end
        end
        
        % --- Formatting ---
        xlabel('Year','FontSize',8); 
        if cv == 1, ylabel(rowTitle, 'FontSize', 10, 'FontWeight', 'bold', 'Color', COL.accent); 
        else, ylabel(cogYLab{cv}, 'FontSize', 7); end
        
        if row == 1, title(cogNames{cv}, 'FontSize', 10, 'FontWeight', 'bold', 'Color', COL.accent); end
        
        set(ax,'FontSize',7,'Box','off','XGrid','on','YGrid','on', 'GridColor',COL.grid,'GridAlpha',0.4,'TickDir','out');
        xlim([-0.3 max(years)+0.3]);
        
        % Expand Y-limits to fit stats and plot 0-line for Z-scores
        yl = ylim; 
        if cv > 1, plot([-1, max(years)+1], [0 0], 'k--', 'Color', [0.7 0.7 0.7], 'LineWidth', 0.5); end
        ylim([yl(1), yl(2) + diff(yl)*0.18]);
        
        % Add Legend
        if cv == 1
            valid_ph = ph5(isgraphics(ph5));
            valid_lg = lgNames(isgraphics(ph5));
            lg = legend(valid_ph, valid_lg, 'Location','best','FontSize',6,'Box','off');
            lg.ItemTokenSize = [12 12];
        end
        hold off;
    end
end

sgtitle('Longitudinal Cognitive Domain Decline Across Subgroups', 'FontSize', 14, 'FontWeight', 'bold', 'Color', COL.accent);
exportFigure(fig5, fullfile(outputDir, 'Fig5_Cognitive_Composites_MegaGrid'));
fprintf('    -> Fig5_Cognitive_Composites_MegaGrid\n\n');

%%
%% ===  14. FIG 8: COVARIATE-ADJUSTED WITHIN-SUBJECT CHANGE (Yr 0 -> 4) ==
fprintf('[14] Figure 8: Within-subject cognitive change (Adjusted)...\n');

bl_t = T(safeCol(T,'YEAR')==0 & isPD & ~isnan(T.group4), :);
y4_t = T(safeCol(T,'YEAR')==4 & isPD, :);

[~, ia, ib] = intersect(bl_t.PATNO, y4_t.PATNO);
fprintf('    Matched Year 0->4: %d participants\n', numel(ia));

% 1. Create Baseline and Year 4 Z-Scored Composites
rawVars = {'DVT_TOTAL_RECALL', 'DVT_RECOG_DISC_INDEX', 'lns', 'VLTANIM', 'SDMTOTAL', 'TMT_B', 'bjlot'};
dirs    = [1, 1, 1, 1, 1, -1, 1]; 

for i = 1:numel(rawVars)
    v = rawVars{i};
    valid_bl = safeCol(bl_t(ia, :), v); 
    valid_bl = valid_bl(~isnan(valid_bl));
    mu = mean(valid_bl); sig = std(valid_bl);
    
    bl_t.(['Z_', v]) = dirs(i) * ((safeCol(bl_t, v) - mu) / sig);
    y4_t.(['Z_', v]) = dirs(i) * ((safeCol(y4_t, v) - mu) / sig);
end

bl_t.comp_memory = nanmean([bl_t.Z_DVT_TOTAL_RECALL, bl_t.Z_DVT_RECOG_DISC_INDEX], 2);
bl_t.comp_exec   = nanmean([bl_t.Z_lns, bl_t.Z_VLTANIM, bl_t.Z_SDMTOTAL, bl_t.Z_TMT_B], 2);
bl_t.comp_visuo  = bl_t.Z_bjlot;

y4_t.comp_memory = nanmean([y4_t.Z_DVT_TOTAL_RECALL, y4_t.Z_DVT_RECOG_DISC_INDEX], 2);
y4_t.comp_exec   = nanmean([y4_t.Z_lns, y4_t.Z_VLTANIM, y4_t.Z_SDMTOTAL, y4_t.Z_TMT_B], 2);
y4_t.comp_visuo  = y4_t.Z_bjlot;

% 2. Plotting Setup
chgVars  = {'moca', 'comp_memory', 'comp_exec', 'comp_visuo'};
chgNames = {'Adj. \Delta MoCA', 'Adj. \Delta Memory', 'Adj. \Delta Executive', 'Adj. \Delta Visuospatial'};

if numel(ia) > 20
    fig8 = figure('Units','centimeters','Position',[2 2 22 6],'Color',COL.bg);
    
    for cv = 1:4
        ax = subplot(1,4,cv); set(ax,'Color',COL.bg); hold on;
        
        % Extract raw data
        v0 = safeCol(bl_t(ia,:), chgVars{cv});
        v4 = safeCol(y4_t(ib,:), chgVars{cv});
        delta_raw = v4 - v0; 
        
        grpBL = bl_t.group4(ia);
        
        % Extract Covariates for these specific patients
        age_c = safeCol(bl_t(ia,:), 'age');
        edu_c = safeCol(bl_t(ia,:), 'EDUCYRS'); % Check name
        sex_c = safeCol(bl_t(ia,:), 'SEX');
        if iscell(sex_c) || iscategorical(sex_c) || isstring(sex_c)
            sex_c = double(categorical(sex_c));
        end
        
        % Filter for complete cases to run the adjustment model
        valid = ~isnan(delta_raw) & ~isnan(age_c) & ~isnan(sex_c) & ~isnan(edu_c) & ~isnan(v0) & ~isnan(grpBL);
        
        % 3. Covariate Adjustment (Age, Sex, Educ, + BASELINE SCORE)
        delta_adj = nan(size(delta_raw));
        if sum(valid) > 10
            X_cov = [age_c(valid), sex_c(valid), edu_c(valid), v0(valid)];
            % Check for rank deficiency 
            if size(unique(X_cov, 'rows'), 1) > 4
                mdl = fitlm(X_cov, delta_raw(valid));
                delta_adj(valid) = mdl.Residuals.Raw + nanmean(delta_raw(valid));
            else
                delta_adj(valid) = delta_raw(valid); % Fallback if data is too sparse
            end
        end
        
        % 4. Plotting the Adjusted Deltas
        for g = 1:4
            gd = delta_adj(grpBL==g & valid);
            if isempty(gd), continue; end
            
            jit = 0.15 * (rand(numel(gd),1) - 0.5);
            scatter(g+jit, gd, 12, COL.groups(g,:), 'filled', 'MarkerFaceAlpha', 0.4);
            
            q = quantile(gd, [0.25 0.5 0.75]);
            if q(3) > q(1)
                rectangle('Position', [g-0.25, q(1), 0.50, q(3)-q(1)], ...
                    'EdgeColor', COL.groups(g,:), 'LineWidth', 1.2, 'Curvature', 0.1);
            end
            plot([g-0.25, g+0.25], [q(2) q(2)], '-', 'Color', COL.groups(g,:), 'LineWidth', 2.5);
        end
        
        yline(0, '--', 'Color', [0.6 0.6 0.6], 'LineWidth', 1);
        
        % Parametric ANOVA on the adjusted continuous residuals
        if sum(valid) > 10
            p = anova1(delta_adj(valid), grpBL(valid), 'off');
            text(0.5, 0.97, sprintf('p = %.3f', p), 'Units', 'normalized', ...
                'HorizontalAlignment', 'center', 'FontSize', 8, 'FontWeight', 'bold', 'Color', COL.txt);
        end
        
        % Formatting
        set(ax,'XTick',1:4,'XTickLabel',groupLabelsShort,'FontSize',7,'Box','off', ...
            'YGrid','on','GridColor',COL.grid,'TickDir','out');
        xtickangle(25); 
        ylabel(chgNames{cv}, 'FontSize', 9, 'FontWeight', 'bold');
        title(chgNames{cv}, 'FontSize', 9, 'FontWeight', 'bold', 'Color', COL.accent);
        
        yl = ylim;
        ylim([yl(1) - diff(yl)*0.05, yl(2) + diff(yl)*0.1]);
        hold off;
    end
    
    sgtitle('Covariate-Adjusted Cognitive Change (Year 0 \rightarrow Year 4)', 'FontSize', 12, 'FontWeight', 'bold', 'Color', COL.accent);
    exportFigure(fig8, fullfile(outputDir, 'Fig8_CogChange_Adjusted'));
end
%%
%% ========================================================================
%  23. LINEAR MIXED-EFFECTS MODELS: Apathy × ICB × Time on Cognition
%      (Continuous apathy & ICB scores)
%  Outcome 1: MoCA total
%  Outcome 2: Cognitive composite (PC1 of neuropsych battery)
%  Model A (unadjusted): Outcome ~ Apathy*ICB*Time + (1+Time|PATNO)
%  Model B (adjusted):   + Sex + Education + nMissing + Duration + LEDD
%                         + HY + Depression
%  Apathy and ICB defined at BASELINE (time-invariant, z-scored) to avoid
%  reverse causation with time-varying cognition.
% =========================================================================
fprintf('[23] Linear mixed-effects models: Apathy x ICB x Time (continuous)...\n');

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

% --- Cognitive composite via PCA -----------------------------------------
fprintf('    Computing cognitive composite (PCA of neuropsych battery)...\n');
pcaVars = {'lns','SDMTOTAL','VLTANIM','DVT_TOTAL_RECALL','bjlot'};
pcaNames = {'LNS','SDMT','Sem.Fluency','HVLT Recall','Benton JLO'};
pcaMat = nan(height(Tpd_lme), numel(pcaVars));
for vi = 1:numel(pcaVars)
    pcaMat(:,vi) = safeCol(Tpd_lme, pcaVars{vi});
end

% Z-score each test using the full PD sample (pooled mean/SD)
pcaMu = mean(pcaMat, 1, 'omitnan');
pcaSD = std(pcaMat, 0, 1, 'omitnan');
pcaZ = (pcaMat - pcaMu) ./ pcaSD;

% PCA on complete rows only, then project all rows
completeRows = all(~isnan(pcaZ), 2);
nComplete = sum(completeRows);
fprintf('    Complete cases for PCA: %d / %d\n', nComplete, height(Tpd_lme));

if nComplete > 100
    [coeff, ~, ~, ~, explained] = pca(pcaZ(completeRows,:));
    fprintf('    PC1 explains %.1f%% of variance\n', explained(1));
    fprintf('    PC1 loadings:\n');
    for vi = 1:numel(pcaVars)
        fprintf('      %-15s : %.3f\n', pcaNames{vi}, coeff(vi,1));
    end
    Tpd_lme.cog_pc1 = nan(height(Tpd_lme), 1);
    Tpd_lme.cog_pc1(completeRows) = pcaZ(completeRows,:) * coeff(:,1);
    % Flip sign if higher PC1 = worse cognition
    if mean(coeff(:,1)) < 0
        Tpd_lme.cog_pc1 = -Tpd_lme.cog_pc1;
        coeff(:,1) = -coeff(:,1);
        fprintf('    (PC1 sign flipped so higher = better cognition)\n');
    end
else
    Tpd_lme.cog_pc1 = nan(height(Tpd_lme), 1);
    fprintf('    WARNING: Insufficient complete cases for PCA\n');
end

% --- Build fitlme-compatible table ---------------------------------------
lmeT = table( ...
    Tpd_lme.PATNO, ...
    Tpd_lme.time, ...
    Tpd_lme.moca_val, ...
    Tpd_lme.cog_pc1, ...
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
    'VariableNames', {'PATNO','time','moca','cog_pc1','apathy','icb', ...
                      'depfactor','sex','educ','nmiss','duration','ledd','hy','age','updrs3_score'});

% Sex remains categorical
lmeT.sex = categorical(lmeT.sex);

% Valid rows
okBase = ~isnan(lmeT.apathy) & ~isnan(lmeT.icb) & ~isnan(lmeT.time);

% Helper: coefficient table extraction
coeffToTable = @(mdl) lme_coeff_to_table(mdl);

% ---- Model 1: MoCA -------------------------------------------------------
fprintf('\n    ============================================================\n');
fprintf('    MODEL 1: MoCA ~ Apathy(z) * ICB(z) * Time\n');
fprintf('    ============================================================\n');

okM1 = okBase & ~isnan(lmeT.moca);
lmeT1 = lmeT(okM1, :);
fprintf('    Observations: %d  |  Patients: %d\n', height(lmeT1), numel(unique(lmeT1.PATNO)));

try
    % Model 1A: Unadjusted
    fprintf('\n    --- Model 1A (Unadjusted) ---\n');
    mdl1a = fitlme(lmeT1, ...
        'moca ~ apathy*icb*time + (1 + time | PATNO)', ...
        'FitMethod','REML');
    disp(mdl1a);
    fe1a = coeffToTable(mdl1a);
    fprintf('\n    Key effects (unadjusted):\n');
    for ri = 1:height(fe1a)
        nm = fe1a.Name{ri};
        if contains(nm,'time') || contains(nm,'apathy') || contains(nm,'icb')
            fprintf('      %-35s  b=%+7.3f  SE=%.3f  t=%.2f  p=%.4f %s\n', ...
                nm, fe1a.Estimate(ri), fe1a.SE(ri), fe1a.tStat(ri), fe1a.pValue(ri), ...
                pToStars(fe1a.pValue(ri)));
        end
    end

    % Model 1B: Adjusted
    fprintf('\n    --- Model 1B (Adjusted, incl. depression) ---\n');
    okAdj = okM1 & ~isnan(lmeT.educ) & ~isnan(lmeT.nmiss) & ...
            ~isnan(lmeT.duration) & ~isnan(lmeT.ledd) & ~isnan(lmeT.updrs3_score) & ...
            ~isundefined(lmeT.sex) & ~isnan(lmeT.depfactor);
    lmeT1b = lmeT(okAdj, :);
    fprintf('    Observations: %d  |  Patients: %d\n', height(lmeT1b), numel(unique(lmeT1b.PATNO)));

    mdl1b = fitlme(lmeT1b, ...
        'moca ~ apathy*icb*time + age + sex + educ + nmiss + duration + ledd + updrs3_score + depfactor + (1 + time | PATNO)', ...
        'FitMethod','REML');
    disp(mdl1b);
    fe1b = coeffToTable(mdl1b);
    fprintf('\n    Key effects (adjusted):\n');
    for ri = 1:height(fe1b)
        nm = fe1b.Name{ri};
        if contains(nm,'time') || contains(nm,'apathy') || contains(nm,'icb')
            fprintf('      %-35s  b=%+7.3f  SE=%.3f  t=%.2f  p=%.4f %s\n', ...
                nm, fe1b.Estimate(ri), fe1b.SE(ri), fe1b.tStat(ri), fe1b.pValue(ri), ...
                pToStars(fe1b.pValue(ri)));
        end
    end

    % LRT: 3-way vs 2-way
    fprintf('\n    --- Likelihood ratio test: 3-way vs 2-way model ---\n');
    mdl1b_2way = fitlme(lmeT1b, ...
        'moca ~ apathy*icb + apathy*time + icb*time + age + sex + educ + nmiss + duration + ledd + updrs3_score + depfactor + (1 + time | PATNO)', ...
        'FitMethod','ML');
    mdl1b_3way = fitlme(lmeT1b, ...
        'moca ~ apathy*icb*time + age + sex + educ + nmiss + duration + ledd + updrs3_score + depfactor + (1 + time | PATNO)', ...
        'FitMethod','ML');
    comp1 = compare(mdl1b_2way, mdl1b_3way);
    fprintf('    Chi2 = %.2f, df = %d, p = %.4f %s\n', ...
        comp1.LRStat(2), comp1.deltaDF(2), comp1.pValue(2), pToStars(comp1.pValue(2)));

    writetable(fe1a, fullfile(outputDir, 'Table6_LME_MoCA_Unadjusted.csv'));
    writetable(fe1b, fullfile(outputDir, 'Table7_LME_MoCA_Adjusted.csv'));
    fprintf('    -> Table6_LME_MoCA_Unadjusted.csv\n');
    fprintf('    -> Table7_LME_MoCA_Adjusted.csv\n');

catch ME
    fprintf('    ERROR in MoCA LME: %s\n', ME.message);
    fprintf('    (requires Statistics & Machine Learning Toolbox)\n');
end

% ---- Model 2: Cognitive Composite (PC1) ----------------------------------
fprintf('\n    ============================================================\n');
fprintf('    MODEL 2: CogPC1 ~ Apathy(z) * ICB(z) * Time\n');
fprintf('    ============================================================\n');

okM2 = okBase & ~isnan(lmeT.cog_pc1);
lmeT2 = lmeT(okM2, :);
fprintf('    Observations: %d  |  Patients: %d\n', height(lmeT2), numel(unique(lmeT2.PATNO)));

try
    % Model 2A: Unadjusted
    fprintf('\n    --- Model 2A (Unadjusted) ---\n');
    mdl2a = fitlme(lmeT2, ...
        'cog_pc1 ~ apathy*icb*time + (1 + time | PATNO)', ...
        'FitMethod','REML');
    disp(mdl2a);
    fe2a = coeffToTable(mdl2a);
    fprintf('\n    Key effects (unadjusted):\n');
    for ri = 1:height(fe2a)
        nm = fe2a.Name{ri};
        if contains(nm,'time') || contains(nm,'apathy') || contains(nm,'icb')
            fprintf('      %-35s  b=%+7.3f  SE=%.3f  t=%.2f  p=%.4f %s\n', ...
                nm, fe2a.Estimate(ri), fe2a.SE(ri), fe2a.tStat(ri), fe2a.pValue(ri), ...
                pToStars(fe2a.pValue(ri)));
        end
    end

    % Model 2B: Adjusted
    fprintf('\n    --- Model 2B (Adjusted, incl. depression) ---\n');
    okAdj2 = okM2 & ~isnan(lmeT.educ) & ~isnan(lmeT.nmiss) & ...
             ~isnan(lmeT.duration) & ~isnan(lmeT.ledd) & ~isnan(lmeT.updrs3_score) & ...
             ~isundefined(lmeT.sex) & ~isnan(lmeT.depfactor) & ~isnan(lmeT.age);
    lmeT2b = lmeT(okAdj2, :);
    fprintf('    Observations: %d  |  Patients: %d\n', height(lmeT2b), numel(unique(lmeT2b.PATNO)));

    mdl2b = fitlme(lmeT2b, ...
        'cog_pc1 ~ apathy*icb*time + age + sex + educ + nmiss + duration + ledd + updrs3_score + depfactor + (1 + time | PATNO)', ...
        'FitMethod','REML');
    disp(mdl2b);
    fe2b = coeffToTable(mdl2b);
    fprintf('\n    Key effects (adjusted):\n');
    for ri = 1:height(fe2b)
        nm = fe2b.Name{ri};
        if contains(nm,'time') || contains(nm,'apathy') || contains(nm,'icb')
            fprintf('      %-35s  b=%+7.3f  SE=%.3f  t=%.2f  p=%.4f %s\n', ...
                nm, fe2b.Estimate(ri), fe2b.SE(ri), fe2b.tStat(ri), fe2b.pValue(ri), ...
                pToStars(fe2b.pValue(ri)));
        end
    end

    % LRT: 3-way vs 2-way
    fprintf('\n    --- Likelihood ratio test: 3-way vs 2-way model ---\n');
    mdl2b_2way = fitlme(lmeT2b, ...
        'cog_pc1 ~ apathy*icb + apathy*time + icb*time + age + sex + educ + nmiss + duration + ledd + updrs3_score + depfactor + (1 + time | PATNO)', ...
        'FitMethod','ML');
    mdl2b_3way = fitlme(lmeT2b, ...
        'cog_pc1 ~ apathy*icb*time + age + sex + educ + nmiss + duration + ledd + updrs3_score + depfactor + (1 + time | PATNO)', ...
        'FitMethod','ML');
    comp2 = compare(mdl2b_2way, mdl2b_3way);
    fprintf('    Chi2 = %.2f, df = %d, p = %.4f %s\n', ...
        comp2.LRStat(2), comp2.deltaDF(2), comp2.pValue(2), pToStars(comp2.pValue(2)));

    writetable(fe2a, fullfile(outputDir, 'Table8_LME_CogPC1_Unadjusted.csv'));
    writetable(fe2b, fullfile(outputDir, 'Table9_LME_CogPC1_Adjusted.csv'));
    fprintf('    -> Table8_LME_CogPC1_Unadjusted.csv\n');
    fprintf('    -> Table9_LME_CogPC1_Adjusted.csv\n');

catch ME
    fprintf('    ERROR in CogPC1 LME: %s\n', ME.message);
    fprintf('    (requires Statistics & Machine Learning Toolbox)\n');
end

% ---- Fig 21: Model-predicted trajectories with CI ribbons ----------------
fprintf('\n    Generating model-predicted trajectory figure...\n');
try
    fig21 = figure('Units','centimeters','Position',[2 2 17.4 8],'Color',COL.bg);
    set(fig21,'DefaultAxesFontName','Helvetica','DefaultTextFontName','Helvetica');

    predTime = linspace(0, 5, 100)';
    % Four prototypical profiles: Low/High apathy x Low/High ICB at +/-1 SD
    predLabels = {'Low Ap / Low ICB','High Ap / Low ICB', ...
                  'Low Ap / High ICB','High Ap / High ICB'};
    apZvals  = [-1  1 -1  1];   % z-scored: -1 SD / +1 SD
    icbZvals = [-1 -1  1  1];
    predCols = COL.groups;       % Okabe-Ito 4-group palette

    % Collect adjusted models and their data
    mdls = {}; mdlData = {}; yLabels = {}; panelTitles = {};
    if exist('mdl1b','var')
        mdls{end+1} = mdl1b; mdlData{end+1} = lmeT1b;
        yLabels{end+1} = 'Predicted MoCA';
        panelTitles{end+1} = 'A   MoCA';
    end
    if exist('mdl2b','var')
        mdls{end+1} = mdl2b; mdlData{end+1} = lmeT2b;
        yLabels{end+1} = 'Predicted Cog PC1 (z)';
        panelTitles{end+1} = 'B   Cognitive Composite';
    end
    nPanels = numel(mdls);

    for mi = 1:nPanels
        ax = subplot(1, nPanels, mi);
        set(ax, 'Color', COL.bg, 'FontSize', 7); hold on;
        dt = mdlData{mi};

        % Covariate medians for prediction
        medEduc = median(dt.educ,'omitnan');
        medNmiss = median(dt.nmiss,'omitnan');
        medDur  = median(dt.duration,'omitnan');
        medLedd = median(dt.ledd,'omitnan');
        medHy   = median(dt.hy,'omitnan');
        medDep  = median(dt.depfactor,'omitnan');
        sexCats = categories(dt.sex); sexCnts = countcats(dt.sex);
        [~,mxI] = max(sexCnts); modSex = sexCats{mxI};

        nPT = numel(predTime);
        phP = gobjects(4,1);

        for gi = 1:4
            predTbl = table( ...
                repmat(apZvals(gi), nPT, 1), ...
                repmat(icbZvals(gi), nPT, 1), ...
                predTime, ...
                repmat(categorical(cellstr(modSex)), nPT, 1), ...
                repmat(medEduc, nPT, 1), ...
                repmat(medNmiss, nPT, 1), ...
                repmat(medDur, nPT, 1), ...
                repmat(medLedd, nPT, 1), ...
                repmat(medupdrs3_score, nPT, 1), ...
                repmat(medDep, nPT, 1), ...
                repmat(dt.PATNO(1), nPT, 1), ...
                'VariableNames', {'apathy','icb','time','sex','educ', ...
                                  'nmiss','duration','ledd','updrs3_score','depfactor','PATNO'});

            [yPred, yCI] = predict(mdls{mi}, predTbl, 'Conditional', false, ...
                                   'Simultaneous', false);

            % Shaded 95% CI ribbon
            fill(ax, [predTime; flipud(predTime)], ...
                 [yCI(:,1); flipud(yCI(:,2))], ...
                 predCols(gi,:), 'FaceAlpha', 0.10, 'EdgeColor', 'none', ...
                 'HandleVisibility', 'off');

            % Trajectory line
            phP(gi) = plot(ax, predTime, yPred, '-', ...
                'Color', predCols(gi,:), 'LineWidth', 1.8);
        end

        % Axis formatting
        xlabel('Years from Baseline', 'FontSize', 8);
        ylabel(yLabels{mi}, 'FontSize', 8);
        title(panelTitles{mi}, 'FontSize', 9, 'FontWeight', 'bold', ...
              'Color', COL.accent, 'HorizontalAlignment', 'left', ...
              'Units', 'normalized', 'Position', [0 1.04 0]);
        xlim([-0.2 5.2]);
        formatAxis(ax, COL);

        % Legend (bottom-left, compact)
        if mi == nPanels
            lg = legend(phP, predLabels, 'Location', 'southwest', ...
                'FontSize', 5.5, 'Box', 'off', 'NumColumns', 1);
            lg.ItemTokenSize = [12 8];
        end
    end

    % Annotation: what +/- 1SD means in original units
    annotation(fig21, 'textbox', [0.01 0.005 0.98 0.04], ...
        'String', sprintf('Profiles at \\pm1 SD: Apathy composite (M=%.1f, SD=%.1f), QUIP total (M=%.1f, SD=%.1f). Covariates held at sample median.', ...
                          apMu, apSD, icbMu, icbSD), ...
        'FontSize', 5.5, 'FontName', 'Helvetica', 'Color', [0.5 0.5 0.5], ...
        'EdgeColor', 'none', 'HorizontalAlignment', 'center', ...
        'FitBoxToText', 'off', 'VerticalAlignment', 'bottom');

    exportFigure(fig21, fullfile(outputDir, 'Fig21_LME_Predicted_Trajectories'));
    fprintf('    -> Fig21_LME_Predicted_Trajectories\n');
catch ME
    fprintf('    Prediction figure skipped: %s\n', ME.message);
end

% ---- Fig 22: Forest plot of interaction terms ----------------------------
fprintf('    Generating forest plot of model coefficients...\n');
try
    if exist('fe1b','var') && exist('fe2b','var')
        fig22 = figure('Units','centimeters','Position',[2 2 17.4 10],'Color',COL.bg);
        set(fig22,'DefaultAxesFontName','Helvetica','DefaultTextFontName','Helvetica');

        panelNames = {'A   MoCA (adjusted)', 'B   Cognitive Composite (adjusted)'};
        feList = {fe1b, fe2b};

        for mi = 1:2
            fe = feList{mi};
            ax = subplot(1, 2, mi);
            set(ax, 'Color', COL.bg, 'FontSize', 7); hold on;

            % Select terms of interest (apathy, icb, time, and their interactions)
            keepIdx = [];
            for ri = 1:height(fe)
                nm = fe.Name{ri};
                if contains(nm,'apathy') || contains(nm,'icb') || strcmp(nm,'time')
                    keepIdx(end+1) = ri; %#ok
                end
            end
            nK = numel(keepIdx);

            % --- Colour palette for significance bands ---
            colSig3 = COL.overlap;         % p < .001
            colSig2 = COL.apathy;          % p < .01
            colSig1 = COL.icb;             % p < .05
            colNS   = [0.65 0.65 0.65];    % p >= .05

            maxAbsCI = 0;
            for ki = 1:nK
                ri = keepIdx(ki);
                ci_lo = fe.Estimate(ri) - 1.96*fe.SE(ri);
                ci_hi = fe.Estimate(ri) + 1.96*fe.SE(ri);
                maxAbsCI = max(maxAbsCI, max(abs([ci_lo ci_hi])));
            end

            for ki = 1:nK
                ri = keepIdx(ki);
                est   = fe.Estimate(ri);
                se    = fe.SE(ri);
                ci_lo = est - 1.96*se;
                ci_hi = est + 1.96*se;
                p     = fe.pValue(ri);
                yp    = nK - ki + 1;

                % Colour by significance tier
                if p < 0.001,     mc = colSig3;
                elseif p < 0.01,  mc = colSig2;
                elseif p < 0.05,  mc = colSig1;
                else,             mc = colNS;
                end

                % Horizontal CI bar
                plot(ax, [ci_lo, ci_hi], [yp yp], '-', ...
                    'Color', mc, 'LineWidth', 2);
                % Point estimate (filled circle)
                plot(ax, est, yp, 'o', 'Color', mc, ...
                    'MarkerFaceColor', mc, 'MarkerSize', 6, 'LineWidth', 0.5);

                % Annotation: b, 95% CI, stars
                annStr = sprintf('  b = %+.3f [%.3f, %.3f] %s', est, ci_lo, ci_hi, pToStars(p));
                text(ax, maxAbsCI*1.08, yp, annStr, ...
                    'FontSize', 5.5, 'Color', COL.txt, ...
                    'VerticalAlignment', 'middle', 'FontName', 'Helvetica');
            end

            % Reference line at zero
            xline(ax, 0, '-', 'Color', [0.75 0.75 0.75], 'LineWidth', 0.6);

            % Clean tick labels: prettify coefficient names
            tickLabels = fe.Name(keepIdx);
            tickLabels = strrep(tickLabels, 'apathy', 'Apathy');
            tickLabels = strrep(tickLabels, 'icb', 'ICB');
            tickLabels = strrep(tickLabels, 'time', 'Time');
            tickLabels = strrep(tickLabels, ':', ' \times ');

            set(ax, 'YTick', 1:nK, 'YTickLabel', flip(tickLabels), ...
                'TickDir', 'out', 'Box', 'off', ...
                'XGrid', 'on', 'GridColor', COL.grid, 'GridAlpha', 0.3, ...
                'YLim', [0.3 nK+0.7]);
            xlabel('Estimate (95% CI)', 'FontSize', 8);
            title(panelNames{mi}, 'FontSize', 9, 'FontWeight', 'bold', ...
                  'Color', COL.accent, 'HorizontalAlignment', 'left', ...
                  'Units', 'normalized', 'Position', [0 1.04 0]);

            formatAxis(ax, COL);
            hold off;
        end

        % Significance legend (bottom-right annotation)
        annotation(fig22, 'textbox', [0.60 0.005 0.39 0.04], ...
            'String', '*** p<.001   ** p<.01   * p<.05   ns p\geq.05', ...
            'FontSize', 5.5, 'FontName', 'Helvetica', 'Color', [0.5 0.5 0.5], ...
            'EdgeColor', 'none', 'HorizontalAlignment', 'right', ...
            'FitBoxToText', 'off');

        exportFigure(fig22, fullfile(outputDir, 'Fig22_LME_Forest_Plot'));
        fprintf('    -> Fig22_LME_Forest_Plot\n');
    end
catch ME
    fprintf('    Forest plot skipped: %s\n', ME.message);
end

% ---- Summary table of all models ----------------------------------------
try
    if exist('mdl1a','var') && exist('mdl1b','var') && exist('mdl2a','var') && exist('mdl2b','var')
        sumT = table( ...
            {'MoCA unadj.';'MoCA adj.';'CogPC1 unadj.';'CogPC1 adj.'}, ...
            [mdl1a.NumObservations; mdl1b.NumObservations; mdl2a.NumObservations; mdl2b.NumObservations], ...
            [numel(unique(lmeT1.PATNO)); numel(unique(lmeT1b.PATNO)); ...
             numel(unique(lmeT2.PATNO)); numel(unique(lmeT2b.PATNO))], ...
            [mdl1a.LogLikelihood; mdl1b.LogLikelihood; mdl2a.LogLikelihood; mdl2b.LogLikelihood], ...
            [mdl1a.ModelCriterion.AIC; mdl1b.ModelCriterion.AIC; mdl2a.ModelCriterion.AIC; mdl2b.ModelCriterion.AIC], ...
            [mdl1a.ModelCriterion.BIC; mdl1b.ModelCriterion.BIC; mdl2a.ModelCriterion.BIC; mdl2b.ModelCriterion.BIC], ...
            'VariableNames', {'Model','N_obs','N_patients','LogLik','AIC','BIC'});
        disp(sumT);
        % writetable(sumT, fullfile(outputDir, 'Table10_LME_Model_Summary.csv'));
        fprintf('    -> Table10_LME_Model_Summary.csv\n');
    end
catch ME
    fprintf('    Model summary table skipped: %s\n', ME.message);
end
fprintf('\n');

%% Using G Factors
%% ========================================================================
%% ========================================================================
%  23b. LME MODELS USING CFA FACTOR SCORES (g + domain factors)
%  Drop-in replacement for the PCA composite section.
%  Run AFTER the R CFA script has exported CFA_factor_scores.csv
% =========================================================================
fprintf('[23b] Importing CFA factor scores and running domain-specific LMEs...\n');

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

fprintf('    Factor scores merged. Available: bifactor_g, bifactor_att_exec, bifactor_memory, bifactor_visuo\n');
fprintf('    Non-missing bifactor_g: %d / %d rows\n', ...
    sum(~isnan(Tpd_lme.bifactor_g)), height(Tpd_lme));

% ---- Key join of the factor scores (2026-10-08). lmeT rows follow the curated-cut order
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
lmeT.bifactor_att_exec = Tpd_lme.bifactor_att_exec(w219_loc);
lmeT.bifactor_memory = Tpd_lme.bifactor_memory(w219_loc);
lmeT.bifactor_visuo = Tpd_lme.bifactor_visuo(w219_loc);

try
lmeT.obs_g = Tpd_lme.obs_g(w219_loc);
lmeT.obs_att_exec_resid = Tpd_lme.obs_att_exec_resid(w219_loc);
lmeT.obs_memory_resid = Tpd_lme.obs_memory_resid(w219_loc);
lmeT.obs_visuo_resid = Tpd_lme.obs_visuo_resid(w219_loc);

lmeT.obs_att_exec = Tpd_lme.obs_att_exec(w219_loc);
lmeT.obs_memory = Tpd_lme.obs_memory(w219_loc);
lmeT.obs_visuo = Tpd_lme.obs_visuo(w219_loc);

lmeT.cfa_g = Tpd_lme.cfa_g(w219_loc);
lmeT.cfa_att_exec = Tpd_lme.cfa_att_exec(w219_loc);
lmeT.cfa_memory = Tpd_lme.cfa_memory(w219_loc);
lmeT.cfa_visuo = Tpd_lme.cfa_visuo(w219_loc);


end

% --- Add APOE e4 to LME table ---
apoe4_raw = safeCol(Tpd_lme, 'APOE_e4');
apoe4_raw = apoe4_raw(w219_loc);
lmeT.apoe4 = double(apoe4_raw >= 1);  % Binarize: 0=non-carrier, 1=carrier
lmeT.apoe4(isnan(apoe4_raw)) = NaN;
% --- Run LMEs for each cognitive outcome -----------------------------------
% cogOutcomes = {'bifactor_g', 'bifactor_att_exec', 'bifactor_memory', 'bifactor_visuo'};
cogOutcomes = {'bifactor_g', 'obs_att_exec_resid', 'obs_memory_resid', 'obs_visuo_resid'};

cogLabels   = {'General Cognition (g)', 'Attention/Executive', ...
               'Memory', 'Visuospatial'};
tableNums   = {'Table8','Table9','Table10','Table11'};

allFE = {};  % collect for combined forest plot

for ci = 1:numel(cogOutcomes)
    outVar = cogOutcomes{ci};
    outLabel = cogLabels{ci};


    

    fprintf('\n    ============================================================\n');
    fprintf('    MODEL: %s ~ Apathy(z) * ICB(z) * Time\n', outLabel);
    fprintf('    ============================================================\n');

    okM = okBase & ~isnan(lmeT.(outVar));
    lmeD = lmeT(okM, :);
    fprintf('    Observations: %d  |  Patients: %d\n', ...
        height(lmeD), numel(unique(lmeD.PATNO)));

    try
        % Unadjusted
        formula_unadj = sprintf('%s ~ apathy*icb*time + (1 + time | PATNO)', outVar);
        mdl_unadj = fitlme(lmeD, formula_unadj, 'FitMethod','REML');
        fe_unadj = coeffToTable(mdl_unadj);

        % Adjusted
        okAdj = okM & ~isnan(lmeT.educ) & ~isnan(lmeT.nmiss) & ...
                ~isnan(lmeT.duration) & ~isnan(lmeT.ledd) & ~isnan(lmeT.hy) & ...
                ~isundefined(lmeT.sex) & ~isnan(lmeT.depfactor) & ~isnan(lmeT.age);
        lmeD_adj = lmeT(okAdj, :);
        fprintf('    Adjusted n: %d obs, %d patients\n', ...
            height(lmeD_adj), numel(unique(lmeD_adj.PATNO)));

        formula_adj = sprintf('%s ~ apathy*icb*time + age + sex + educ + nmiss + duration + ledd + updrs3_score + depfactor + (1 + time | PATNO)', outVar);
        mdl_adj = fitlme(lmeD_adj, formula_adj, 'FitMethod','REML');
        fe_adj = coeffToTable(mdl_adj);

        fprintf('\n    Key effects (adjusted):\n');
        for ri = 1:height(fe_adj)
            nm = fe_adj.Name{ri};
            if contains(nm,'time') || contains(nm,'apathy') || contains(nm,'icb')
                fprintf('      %-35s  b=%+7.3f  SE=%.3f  t=%.2f  p=%.4f %s\n', ...
                    nm, fe_adj.Estimate(ri), fe_adj.SE(ri), fe_adj.tStat(ri), ...
                    fe_adj.pValue(ri), pToStars(fe_adj.pValue(ri)));
            end
        end

        % LRT: 3-way vs 2-way
        formula_2way = sprintf('%s ~ apathy*icb + apathy*time + icb*time + sex + educ + nmiss + duration + ledd + updrs3_score + depfactor + (1 + time | PATNO)', outVar);
        formula_3way = sprintf('%s ~ apathy*icb*time + age + sex + educ + nmiss + duration + ledd + updrs3_score + depfactor + (1 + time | PATNO)', outVar);
        mdl_2w = fitlme(lmeD_adj, formula_2way, 'FitMethod','ML');
        mdl_3w = fitlme(lmeD_adj, formula_3way, 'FitMethod','ML');
        comp = compare(mdl_2w, mdl_3w);
        fprintf('    LRT 3-way: Chi2=%.2f, df=%d, p=%.4f %s\n', ...
            comp.LRStat(2), comp.deltaDF(2), comp.pValue(2), pToStars(comp.pValue(2)));


        

        % Save
        writetable(fe_adj, fullfile(outputDir, sprintf('%s_LME_%s_Adjusted.csv', tableNums{ci}, outVar)));
        fprintf('    -> %s_LME_%s_Adjusted.csv\n', tableNums{ci}, outVar);

        % Store for combined forest plot
        fe_adj.Domain = repmat({outLabel}, height(fe_adj), 1);
        allFE{ci} = fe_adj;

    catch ME
        fprintf('    ERROR: %s\n', ME.message);
    end
end

% --- Fig 23: Combined forest plot (all four domains) -----------------------
fprintf('\n    Generating combined domain forest plot...\n');
try
    allFEtab = vertcat(allFE{:});

    fig23 = figure('Units','centimeters','Position',[2 2 17.4 14],'Color',COL.bg);
    set(fig23,'DefaultAxesFontName','Helvetica','DefaultTextFontName','Helvetica');

    domainCols = [COL.groups(1,:); COL.apathy; COL.icb; COL.overlap];

    for ci = 1:4
        ax = subplot(2, 2, ci);
        set(ax, 'Color', COL.bg, 'FontSize', 6.5); hold on;

        fe = allFE{ci};
        keepIdx = [];
        for ri = 1:height(fe)
            nm = fe.Name{ri};
            if contains(nm,'apathy') || contains(nm,'icb') || strcmp(nm,'time')
                keepIdx(end+1) = ri; %#ok
            end
        end
        nK = numel(keepIdx);

        colSig  = domainCols(ci,:);
        colNS   = [0.65 0.65 0.65];

        for ki = 1:nK
            ri = keepIdx(ki);
            est   = fe.Estimate(ri);
            ci_lo = est - 1.96*fe.SE(ri);
            ci_hi = est + 1.96*fe.SE(ri);
            p     = fe.pValue(ri);
            yp    = nK - ki + 1;

            mc = colNS;
            if p < 0.05, mc = colSig; end
            if p < 0.001, mc = colSig * 0.8; end  % darker for strongest effects

            plot(ax, [ci_lo ci_hi], [yp yp], '-', 'Color', mc, 'LineWidth', 2);
            plot(ax, est, yp, 'o', 'Color', mc, 'MarkerFaceColor', mc, ...
                'MarkerSize', 5.5, 'LineWidth', 0.5);
        end

        xline(ax, 0, '-', 'Color', [0.75 0.75 0.75], 'LineWidth', 0.5);

        tickLabels = fe.Name(keepIdx);
        tickLabels = strrep(tickLabels, 'apathy', 'Apathy');
        tickLabels = strrep(tickLabels, 'icb', 'ICB');
        tickLabels = strrep(tickLabels, 'time', 'Time');
        tickLabels = strrep(tickLabels, ':', ' \times ');

        panelLetter = char('A' + ci - 1);
        set(ax, 'YTick', 1:nK, 'YTickLabel', flip(tickLabels), ...
            'TickDir', 'out', 'Box', 'off', 'YLim', [0.3 nK+0.7], ...
            'XGrid', 'on', 'GridColor', COL.grid, 'GridAlpha', 0.25);
        xlabel('Estimate (95% CI)', 'FontSize', 7);
        title(sprintf('%s   %s', panelLetter, cogLabels{ci}), ...
            'FontSize', 8.5, 'FontWeight', 'bold', 'Color', domainCols(ci,:), ...
            'HorizontalAlignment', 'left', 'Units', 'normalized', ...
            'Position', [0 1.06 0]);
        formatAxis(ax, COL);
    end

    sgtitle('LME Fixed Effects by Cognitive Domain (CFA factor scores, adjusted)', ...
        'FontSize', 10, 'FontWeight', 'bold', 'Color', COL.accent);
    exportFigure(fig23, fullfile(outputDir, 'Fig23_LME_CFA_Domain_Forest'));
    fprintf('    -> Fig23_LME_CFA_Domain_Forest\n');

catch ME
    fprintf('    Combined forest plot skipped: %s\n', ME.message);
end

fprintf('\n');

%%
%% ========================================================================
%  23b. LME MODELS USING CFA FACTOR SCORES (g + domain factors)
%  Drop-in replacement for the PCA composite section.
%  Run AFTER the R CFA script has exported CFA_factor_scores.csv
% =========================================================================
fprintf('[23b] Importing CFA factor scores and running domain-specific LMEs...\n');

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
fprintf('    Factor scores merged. Available: bifactor_g, bifactor_att_exec, bifactor_memory, bifactor_visuo\n');
fprintf('    Non-missing bifactor_g: %d / %d rows\n', ...
    sum(~isnan(Tpd_lme.bifactor_g)), height(Tpd_lme));

% ---- Key join of the factor scores (2026-10-08). lmeT rows follow the curated-cut order
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
lmeT.bifactor_att_exec = Tpd_lme.bifactor_att_exec(w219_loc);
lmeT.bifactor_memory = Tpd_lme.bifactor_memory(w219_loc);
lmeT.bifactor_visuo = Tpd_lme.bifactor_visuo(w219_loc);

% --- Add APOE e4 to LME table ----------------------------------------------
apoe4_raw = safeCol(Tpd_lme, 'APOE_e4');
apoe4_raw = apoe4_raw(w219_loc);
lmeT.apoe4 = double(apoe4_raw >= 1);  % Binarize: 0=non-carrier, 1=carrier
lmeT.apoe4(isnan(apoe4_raw)) = NaN;

% --- Run LMEs for each cognitive outcome -----------------------------------
cogOutcomes = {'bifactor_g', 'bifactor_att_exec', 'bifactor_memory', 'bifactor_visuo'};
cogLabels   = {'General Cognition (g)', 'Attention/Executive', ...
               'Memory', 'Visuospatial'};
tableNums   = {'Table8','Table9','Table10','Table11'};
allFE = {};  % collect for combined forest plot

for ci = 1:numel(cogOutcomes)
    outVar = cogOutcomes{ci};
    outLabel = cogLabels{ci};
    
    fprintf('\n    ============================================================\n');
    fprintf('    MODEL: %s ~ (APOE4, Apathy, ICB, Time) up to 3-way\n', outLabel);
    fprintf('    ============================================================\n');
    
    % Require APOE4 to be non-missing for the models
    okM = okBase & ~isnan(lmeT.(outVar)) & ~isnan(lmeT.apoe4);
    lmeD = lmeT(okM, :);
    fprintf('    Observations: %d  |  Patients: %d\n', ...
        height(lmeD), numel(unique(lmeD.PATNO)));
        
    try
        % Unadjusted Model (All interactions up to 3-way)
        formula_unadj = sprintf('%s ~ apathy*icb*time + apoe4*apathy*time + apoe4*icb*time + apoe4*apathy*icb + (1 + time | PATNO)', outVar);
        mdl_unadj = fitlme(lmeD, formula_unadj, 'FitMethod','REML');
        fe_unadj = coeffToTable(mdl_unadj);
        
        % Adjusted Model (All interactions up to 3-way + Covariates)
        okAdj = okM & ~isnan(lmeT.educ) & ~isnan(lmeT.nmiss) & ...
                ~isnan(lmeT.duration) & ~isnan(lmeT.ledd) & ~isnan(lmeT.hy) & ...
                ~isundefined(lmeT.sex) & ~isnan(lmeT.depfactor);
        lmeD_adj = lmeT(okAdj, :);
        fprintf('    Adjusted n: %d obs, %d patients\n', ...
            height(lmeD_adj), numel(unique(lmeD_adj.PATNO)));
            
        formula_adj = sprintf('%s ~ apoe4*icb*time + apoe4*apathy*time + sex + educ + nmiss + duration + ledd + hy + depfactor + (1 + time | PATNO)', outVar);
        mdl_adj = fitlme(lmeD_adj, formula_adj, 'FitMethod','REML', 'DummyVarCoding', 'effects');
        fe_adj = coeffToTable(mdl_adj);
        
        fprintf('\n    Key effects (adjusted):\n');
        for ri = 1:height(fe_adj)
            nm = fe_adj.Name{ri};
            % Now includes APOE4 interactions in the console printout
            if contains(nm,'time') || contains(nm,'apathy') || contains(nm,'icb') || contains(nm,'apoe4')
                fprintf('      %-40s  b=%+7.3f  SE=%.3f  t=%.2f  p=%.4f %s\n', ...
                    nm, fe_adj.Estimate(ri), fe_adj.SE(ri), fe_adj.tStat(ri), ...
                    fe_adj.pValue(ri), pToStars(fe_adj.pValue(ri)));
            end
        end
        
        % LRT: 3-way Apathy*ICB*Time vs model without it (controlling for genetic interactions)
        formula_2way = sprintf('%s ~ apathy*icb + apathy*time + icb*time + apoe4*apathy*time + apoe4*icb*time + apoe4*apathy*icb + sex + educ + nmiss + duration + ledd + hy + depfactor + (1 + time | PATNO)', outVar);
        formula_3way = formula_adj;
        mdl_2w = fitlme(lmeD_adj, formula_2way, 'FitMethod','ML');
        mdl_3w = fitlme(lmeD_adj, formula_3way, 'FitMethod','ML');
        comp = compare(mdl_2w, mdl_3w);
        fprintf('    LRT (Apathy:ICB:Time vs 2-way base): Chi2=%.2f, df=%d, p=%.4f %s\n', ...
            comp.LRStat(2), comp.deltaDF(2), comp.pValue(2), pToStars(comp.pValue(2)));
        
        % Save to CSV
        writetable(fe_adj, fullfile(outputDir, sprintf('%s_LME_%s_Adjusted.csv', tableNums{ci}, outVar)));
        fprintf('    -> %s_LME_%s_Adjusted.csv\n', tableNums{ci}, outVar);
        
        % Store for combined forest plot
        fe_adj.Domain = repmat({outLabel}, height(fe_adj), 1);
        allFE{ci} = fe_adj;
        
    catch ME
        fprintf('    ERROR: %s\n', ME.message);
    end
end

% --- Fig 23: Combined forest plot (all four domains) -----------------------
fprintf('\n    Generating combined domain forest plot...\n');
try
    allFEtab = vertcat(allFE{:});
    fig23 = figure('Units','centimeters','Position',[2 2 17.4 16],'Color',COL.bg); % Slightly taller to fit more vars
    set(fig23,'DefaultAxesFontName','Helvetica','DefaultTextFontName','Helvetica');
    domainCols = [COL.groups(1,:); COL.apathy; COL.icb; COL.overlap];
    
    for ci = 1:4
        ax = subplot(2, 2, ci);
        set(ax, 'Color', COL.bg, 'FontSize', 6.5); hold on;
        fe = allFE{ci};
        keepIdx = [];
        
        % Plot all clinical and genetic interaction terms
        for ri = 1:height(fe)
            nm = fe.Name{ri};
            if contains(nm,'apathy') || contains(nm,'icb') || strcmp(nm,'time') || contains(nm,'apoe4')
                keepIdx(end+1) = ri; %#ok
            end
        end
        
        nK = numel(keepIdx);
        colSig  = domainCols(ci,:);
        colNS   = [0.65 0.65 0.65];
        
        for ki = 1:nK
            ri = keepIdx(ki);
            est   = fe.Estimate(ri);
            ci_lo = est - 1.96*fe.SE(ri);
            ci_hi = est + 1.96*fe.SE(ri);
            p     = fe.pValue(ri);
            yp    = nK - ki + 1;
            
            mc = colNS;
            if p < 0.05, mc = colSig; end
            if p < 0.001, mc = colSig * 0.8; end  % darker for strongest effects
            
            plot(ax, [ci_lo ci_hi], [yp yp], '-', 'Color', mc, 'LineWidth', 2);
            plot(ax, est, yp, 'o', 'Color', mc, 'MarkerFaceColor', mc, ...
                'MarkerSize', 5.5, 'LineWidth', 0.5);
                
            % Add significance stars next to the CI
            if p < 0.05
                text(ax, ci_hi + 0.01 * abs(ci_hi - ci_lo + 0.1), yp, ...
                    pToStars(p), 'FontSize', 7, 'Color', mc, 'VerticalAlignment', 'middle', 'FontWeight', 'bold');
            end
        end
        
        xline(ax, 0, '-', 'Color', [0.75 0.75 0.75], 'LineWidth', 0.5);
        
        % Clean up tick labels for readability
        tickLabels = fe.Name(keepIdx);
        tickLabels = strrep(tickLabels, 'apoe4_1', 'APOE4');
        tickLabels = strrep(tickLabels, 'apoe4', 'APOE4');
        tickLabels = strrep(tickLabels, 'apathy_1', 'Apathy');
        tickLabels = strrep(tickLabels, 'apathy', 'Apathy');
        tickLabels = strrep(tickLabels, 'icb_1', 'ICB');
        tickLabels = strrep(tickLabels, 'icb', 'ICB');
        tickLabels = strrep(tickLabels, 'time', 'Time');
        tickLabels = strrep(tickLabels, ':', ' \times ');
        
        panelLetter = char('A' + ci - 1);
        set(ax, 'YTick', 1:nK, 'YTickLabel', flip(tickLabels), ...
            'TickDir', 'out', 'Box', 'off', 'YLim', [0.3 nK+0.7], ...
            'XGrid', 'on', 'GridColor', COL.grid, 'GridAlpha', 0.25);
            
        % Expand X limit slightly to fit asterisks
        xl = xlim(ax);
        xlim(ax, [xl(1), xl(2) + diff(xl)*0.15]);
        
        xlabel('Estimate (95% CI)', 'FontSize', 7);
        title(sprintf('%s   %s', panelLetter, cogLabels{ci}), ...
            'FontSize', 8.5, 'FontWeight', 'bold', 'Color', domainCols(ci,:), ...
            'HorizontalAlignment', 'left', 'Units', 'normalized', ...
            'Position', [0 1.06 0]);
        formatAxis(ax, COL);
    end
    
    sgtitle('LME Fixed Effects by Cognitive Domain (CFA Bifactor scores, adjusted for APOE)', ...
        'FontSize', 10, 'FontWeight', 'bold', 'Color', COL.accent);
    exportFigure(fig23, fullfile(outputDir, 'Fig23_LME_CFA_Domain_Forest'));
    fprintf('    -> Fig23_LME_CFA_Domain_Forest\n');
catch ME
    fprintf('    Combined forest plot skipped: %s\n', ME.message);
end
fprintf('\n');
%% ===  11. FIG 5: LONGITUDINAL COGNITIVE DOMAIN TRAJECTORIES (CFA) =======
fprintf('[11] Figure 5: Longitudinal cognitive domain composites (CFA)...\n');

% -------------------------------------------------------------------------
% 1. FIGURE PREPARATION
% -------------------------------------------------------------------------
% Define the exact CFA column names we extracted from R
cogVars  = {'bifactor_g', 'bifactor_att_exec', 'bifactor_memory', 'bifactor_visuo'};
cogNames = {'Global Cognition (g)', 'Executive & Attention', 'Memory', 'Visuospatial Function'};
cogYLab  = {'CFA Score (z)', 'CFA Score (z)', 'CFA Score (z)', 'CFA Score (z)'};

years = 0:5;

% Ensure we are using the dataset that has the CFA scores and Stability Groups
% Note: If your StabilityGroups are in 'T' or 'Tpd', make sure the CFA scores 
% are merged into that same table before this step.
TpdSub = Tpd_lme;

% Colors
subColors  = [0.4 0.6 0.4; 0.8 0.3 0.3]; 
pathColors = [0.4 0.4 0.7; 0.8 0.3 0.3]; 
detColors  = [0.7 0.7 0.7; 0.3 0.5 0.7; 0.8 0.5 0.3; 0.6 0.4 0.7; 0.7 0.3 0.4]; 

% Create a massive 4-row, 4-column figure
fig5 = figure('Units','centimeters','Position',[2 2 28 24],'Color',COL.bg);
nRows = 4;

% -------------------------------------------------------------------------
% 2. MEGA-GRID PLOTTING LOOP
% -------------------------------------------------------------------------
for row = 1:1
    
    % Setup logic for the current row
    if row == 1
        nGrps = 4; cMap = COL.groups(1:4,:); 
        rowTitle = 'Original Baseline Groups'; lgNames = groupLabelsShort;
    elseif row == 2
        nGrps = 2; cMap = subColors; 
        rowTitle = 'Stable vs Unstable'; lgNames = {'Stable', 'Unstable'};
    elseif row == 3
        nGrps = 2; cMap = pathColors; 
        rowTitle = 'Persistent Path. vs Unstable'; lgNames = {'Persist Path.', 'Unstable'};
    elseif row == 4
        nGrps = 5; cMap = detColors; 
        rowTitle = 'Granular Trajectories'; lgNames = {'Persist Neither', 'Persist Apathy', 'Persist ICB', 'Apathy \rightarrow ICB', 'ICB \rightarrow Apathy'};
    end
    
    for cv = 1:numel(cogVars)
        ax = subplot(nRows, numel(cogVars), (row-1)*numel(cogVars) + cv);
        set(ax,'Color',COL.bg); hold on;
        
        ph5 = gobjects(nGrps,1); 
        all_mY = nan(nGrps, numel(years)); 
        all_seY = nan(nGrps, numel(years)); 
        
        % --- Calculate and Plot Trajectories ---
        for g = 1:nGrps
            mY = nan(1, numel(years)); semY = nan(1, numel(years)); nY = nan(1, numel(years));
            
            for yi = 1:numel(years)
                % Group indexing logic
                if row == 1, idx = (TpdSub.group4 == g);
                elseif row == 2, idx = (TpdSub.StabilityGroup == g);
                elseif row == 3
                    if g == 1, idx = (TpdSub.DetGroup == 2 | TpdSub.DetGroup == 3); 
                    elseif g == 2, idx = (TpdSub.StabilityGroup == 2); 
                    end
                elseif row == 4, idx = (TpdSub.DetGroup == g);
                end
                
                yr_idx = idx & (TpdSub.YEAR == years(yi));
                v = safeCol(TpdSub(yr_idx, :), cogVars{cv});
                v = v(~isnan(v));
                
                nY(yi) = numel(v);
                if nY(yi) > 0
                    mY(yi) = mean(v); 
                    semY(yi) = std(v)/sqrt(nY(yi)); 
                end
            end
            
            vYrs = nY > 0;
            all_mY(g,:) = mY; all_seY(g,:) = semY;
            ph5(g) = plotTrajectory(ax, years, mY, semY, nY, vYrs, cMap(g,:), 'endpoints');
        end
        
        % --- Add Year-by-Year Statistical Significance ---
        for yi = 1:numel(years)
            grpData = cell(1, nGrps);
            for g = 1:nGrps
                if row == 1, idx = (TpdSub.group4 == g);
                elseif row == 2, idx = (TpdSub.StabilityGroup == g);
                elseif row == 3
                    if g==1, idx = (TpdSub.DetGroup == 2 | TpdSub.DetGroup == 3); 
                    else, idx = (TpdSub.StabilityGroup == 2); 
                    end
                elseif row == 4, idx = (TpdSub.DetGroup == g);
                end
                
                v = safeCol(TpdSub(idx & (TpdSub.YEAR == years(yi)), :), cogVars{cv});
                grpData{g} = v(~isnan(v));
            end
            
            y_loc = nanmax(all_mY(:,yi) + all_seY(:,yi)) + 0.10 * abs(nanmax(all_mY(:,yi)));
            if isnan(y_loc) || y_loc == 0, y_loc = 0.2; end % Fallback for pure 0
            
            % 2-Group t-tests
            if nGrps == 2 && length(grpData{1}) > 2 && length(grpData{2}) > 2
                [~, p_val] = ttest2(grpData{1}, grpData{2});
                if p_val < 0.05
                    if p_val < 0.001, p_str = '***'; elseif p_val < 0.01, p_str = '**'; else, p_str = '*'; end
                    text(years(yi), y_loc, p_str, 'HorizontalAlignment', 'center', 'FontSize', 12, 'Color', COL.txt, 'FontWeight', 'bold');
                end
            end
            
            % Multi-Group ANOVA
            if nGrps > 2
                allVals = []; allGrps = [];
                for g = 1:nGrps
                    if length(grpData{g}) > 1
                        allVals = [allVals; grpData{g}]; %#ok
                        allGrps = [allGrps; repmat(g, length(grpData{g}), 1)]; %#ok
                    end
                end
                if length(unique(allGrps)) > 1 
                    p_val = anova1(allVals, allGrps, 'off');
                    if p_val < 0.05
                        text(years(yi), y_loc, sprintf('p=%.3f', p_val), 'HorizontalAlignment', 'center', 'FontSize', 7, 'Color', COL.accent, 'FontWeight', 'bold');
                    end
                end
            end
        end
        
        % --- Formatting ---
        xlabel('Year','FontSize',8); 
        
        % Y-Axis Labels (Only on the far left column)
        if cv == 1
            ylabel({rowTitle, cogYLab{cv}}, 'FontSize', 9, 'FontWeight', 'bold', 'Color', COL.accent); 
        end
        
        % Titles (Only on the top row)
        if row == 1
            title(cogNames{cv}, 'FontSize', 11, 'FontWeight', 'bold', 'Color', COL.accent); 
        end
        
        set(ax,'FontSize',7,'Box','off','XGrid','on','YGrid','on', 'GridColor',COL.grid,'GridAlpha',0.4,'TickDir','out');
        xlim([-0.3 max(years)+0.3]);
        
        % Plot 0-line for Z-scores
        plot([-1, max(years)+1], [0 0], '--', 'Color', [0.6 0.6 0.6], 'LineWidth', 0.8); 
        
        % Expand Y-limits to fit stats and plot
        yl = ylim; 
        ylim([yl(1), yl(2) + diff(yl)*0.18]);
        
        % Add Legend
        if cv == 1
            valid_ph = ph5(isgraphics(ph5));
            valid_lg = lgNames(isgraphics(ph5));
            lg = legend(valid_ph, valid_lg, 'Location','southwest','FontSize',6,'Box','off');
            lg.ItemTokenSize = [10 10];
        end
        hold off;
    end
end

sgtitle('Longitudinal Cognitive Decline (CFA Factors) Across Motivational Phenotypes', ...
    'FontSize', 14, 'FontWeight', 'bold', 'Color', COL.accent);
    
exportFigure(fig5, fullfile(outputDir, 'Fig5_CFA_Cognitive_Trajectories_MegaGrid'));
fprintf('    -> Fig5_CFA_Cognitive_Trajectories_MegaGrid\n\n');

%%
%% ===  18. FIG 13: APOE e4 ANALYSIS  ====================================
fprintf('[18] APOE e4 analysis (using CFA Global Cognition)...\n');
% --- QUICK FIX: Map Baseline CFA scores into the 'bl' table ---
if ~ismember('cfa_g', bl.Properties.VariableNames) && exist('Tpd_lme', 'var')
    fprintf('    Mapping Year 0 CFA scores into baseline table...\n');
    
    % Initialize empty target column in bl
    bl.cfa_g = nan(height(bl), 1);
    
    % Force pure double IDs to bypass string/cell/scientific notation formatting
    if iscell(bl.PATNO) || iscategorical(bl.PATNO) || isstring(bl.PATNO)
        bl_patno = str2double(string(bl.PATNO));
    else
        bl_patno = double(bl.PATNO);
    end
    
    if iscell(Tpd_lme.PATNO) || iscategorical(Tpd_lme.PATNO) || isstring(Tpd_lme.PATNO)
        lme_patno = str2double(string(Tpd_lme.PATNO));
    else
        lme_patno = double(Tpd_lme.PATNO);
    end
    
    lme_year = double(safeCol(Tpd_lme, 'YEAR'));
    
    % Direct injection loop: Match Patient ID at Year == 0
    matched_bl = 0;
    for i = 1:height(bl)
        idx = (lme_patno == bl_patno(i)) & (lme_year == 0);
        if any(idx)
            match_idx = find(idx, 1); % Take the first valid match
            bl.cfa_g(i) = Tpd_lme.cfa_g(match_idx);
            matched_bl = matched_bl + 1;
        end
    end
    fprintf('    Successfully mapped CFA baseline scores for %d patients.\n', matched_bl);
end
% --------------------------------------------------------------
% Ensure we are pulling from the table that has the CFA scores
if ~ismember('cfa_g', bl.Properties.VariableNames)
    warning('cfa_g not found in baseline table. Ensure CFA merge happened before this step.');
end

apoe4_raw = safeCol(bl, 'APOE_e4');
apoe4 = double(apoe4_raw >= 1);  % Binarize: 0=non-carrier, 1=carrier (1 or 2 copies)
apoe4(isnan(apoe4_raw)) = NaN;
hasAPOE = sum(~isnan(apoe4)) > 20;

if hasAPOE
    n_carrier = sum(apoe4==1); n_noncarrier = sum(apoe4==0);
    n_homo = sum(apoe4_raw==2);
    fprintf('    APOE e4 carriers: %d (het: %d, homo: %d)  |  Non-carriers: %d\n', ...
        n_carrier, n_carrier-n_homo, n_homo, n_noncarrier);
        
    fig13 = figure('Units','centimeters','Position',[2 2 17.4 7],'Color',COL.bg);
    
    % (a) Group prevalence by APOE e4
    ax1 = subplot(1,3,1); set(ax1,'Color',COL.bg); hold on;
    for ae = 0:1
        aeIdx = apoe4 == ae & ~isnan(bl.group4);
        nAE = sum(aeIdx);
        pct = arrayfun(@(g) 100*sum(bl.group4==g & aeIdx)/nAE, 1:4);
        xpos = (1:4) + (ae-0.5)*0.35;
        bar(xpos, pct, 0.30, 'FaceColor', ternary(ae==0,[0.6 0.6 0.6],COL.apathy), ...
            'EdgeColor','none','FaceAlpha',0.8);
    end
    set(ax1,'XTick',1:4,'XTickLabel',groupLabelsShort,'FontSize',7,'Box','off', ...
        'YGrid','on','GridColor',COL.grid,'TickDir','out');
    xtickangle(30); ylabel('Prevalence (%)','FontSize',8);
    lg=legend({sprintf('e4- (n=%d)',n_noncarrier), sprintf('e4+ (n=%d)',n_carrier)}, ...
        'Location','northeast','Box','off','FontSize',6);
    lg.ItemTokenSize=[8 8];
    title('(a) Group by APOE \epsilon4','FontSize',8,'FontWeight','bold','Color',COL.accent);
    
    % (b) Apathy composite by APOE e4 × group
    ax2 = subplot(1,3,2); set(ax2,'Color',COL.bg); hold on;
    for ae = 0:1
        for g = 1:4
            idx = apoe4==ae & bl.group4==g;
            vals = safeCol(bl,'apathy_composite'); vals = vals(idx); vals = vals(~isnan(vals));
            if isempty(vals), continue; end
            xp = g + (ae-0.5)*0.3;
            jit = 0.06*(rand(numel(vals),1)-0.5);
            mcol = ternary(ae==0, COL.groups(g,:)*0.6+0.4, COL.groups(g,:));
            scatter(xp+jit, vals, 10, mcol, 'filled', 'MarkerFaceAlpha',0.3);
            plot([xp-0.08 xp+0.08], [median(vals) median(vals)], '-', 'Color',mcol,'LineWidth',2);
        end
    end
    set(ax2,'XTick',1:4,'XTickLabel',groupLabelsShort,'FontSize',7,'Box','off', ...
        'YGrid','on','GridColor',COL.grid,'TickDir','out');
    xtickangle(30); ylabel('Apathy Composite','FontSize',8);
    title('(b) Apathy × APOE \epsilon4','FontSize',8,'FontWeight','bold','Color',COL.accent);
    
    % (c) CFA Global Cognition by APOE e4 × group
    ax3 = subplot(1,3,3); set(ax3,'Color',COL.bg); hold on;
    for ae = 0:1
        for g = 1:4
            idx = apoe4==ae & bl.group4==g;
            vals = safeCol(bl,'cfa_g'); vals = vals(idx); vals = vals(~isnan(vals));
            if isempty(vals), continue; end
            xp = g + (ae-0.5)*0.3;
            jit = 0.06*(rand(numel(vals),1)-0.5);
            mcol = ternary(ae==0, COL.groups(g,:)*0.6+0.4, COL.groups(g,:));
            scatter(xp+jit, vals, 10, mcol, 'filled', 'MarkerFaceAlpha',0.3);
            plot([xp-0.08 xp+0.08], [median(vals) median(vals)], '-', 'Color',mcol,'LineWidth',2);
        end
    end
    yline(ax3, 0, '--', 'Color', [0.7 0.7 0.7], 'LineWidth', 0.5);
    set(ax3,'XTick',1:4,'XTickLabel',groupLabelsShort,'FontSize',7,'Box','off', ...
        'YGrid','on','GridColor',COL.grid,'TickDir','out');
    xtickangle(30); ylabel('Global Cognition (cfa_g)','FontSize',8,'Interpreter','none');
    title('(c) CFA Cognition × APOE \epsilon4','FontSize',8,'FontWeight','bold','Color',COL.accent);
    
    sgtitle('APOE \epsilon4 Carrier Status Analysis','FontSize',10,'FontWeight','bold','Color',COL.accent);
    exportFigure(fig13, fullfile(outputDir, 'Fig13_APOE_Analysis'));
    fprintf('    -> Fig13_APOE_Analysis\n');
    
    % Statistical tests
    fprintf('\n    APOE e4 effects (Mann-Whitney U):\n');
    apoeTestVars = {'apathy_composite','quip','cfa_g','depression_factor','updrs3_score','stai_trait'};
    apoeTestNames = {'Apathy','QUIP','CFA Global Cog','Depression','Motor','Anxiety'};
    for vi = 1:numel(apoeTestVars)
        v1 = safeCol(bl, apoeTestVars{vi});
        g1 = v1(apoe4==0 & ~isnan(v1)); g2 = v1(apoe4==1 & ~isnan(v1));
        if numel(g1)>5 && numel(g2)>5
            p = ranksum(g1, g2);
            fprintf('      %-15s  e4−: %5.2f±%.2f  e4+: %5.2f±%.2f  p=%.4f %s\n', ...
                apoeTestNames{vi}, mean(g1),std(g1), mean(g2),std(g2), p, pToStars(p));
        end
    end
    
    % Odds ratio: APOE e4 -> apathy
    a11 = sum(apoe4==1 & safeCol(bl,'apathy_composite')>=THRESH.apathy_comp);
    a10 = sum(apoe4==1 & safeCol(bl,'apathy_composite')<THRESH.apathy_comp);
    a01 = sum(apoe4==0 & safeCol(bl,'apathy_composite')>=THRESH.apathy_comp);
    a00 = sum(apoe4==0 & safeCol(bl,'apathy_composite')<THRESH.apathy_comp);
    if all([a11 a10 a01 a00] > 0)
        OR = (a11*a00)/(a10*a01);
        logOR = log(OR); SE_logOR = sqrt(1/a11+1/a10+1/a01+1/a00);
        CI = exp([logOR-1.96*SE_logOR, logOR+1.96*SE_logOR]);
        fprintf('\n    APOE e4 -> Apathy: OR = %.2f (95%% CI: %.2f-%.2f)\n', OR, CI(1), CI(2));
    end
end
fprintf('\n');

%% ===  19. FIG 14: GENETIC SUBGROUP ANALYSIS  ============================
fprintf('[19] Genetic subgroup analysis...\n');
if ismember('subgroup', bl.Properties.VariableNames)
    sgCol = string(bl.subgroup);
    
    sgClean = sgCol;
    sgClean(contains(sgClean,'LRRK2') & ~contains(sgClean,'GBA')) = "LRRK2";
    sgClean(contains(sgClean,'GBA')   & ~contains(sgClean,'LRRK2')) = "GBA";
    sgClean(contains(sgClean,'LRRK2') &  contains(sgClean,'GBA'))  = "LRRK2+GBA";
    sgClean(sgClean == "SNCA")  = "SNCA";
    sgClean(sgClean == "PRKN" | sgClean == "PINK1") = "PRKN/PINK1";
    
    sgUniAll = unique(sgClean);
    sgCounts = arrayfun(@(s) sum(sgClean==s), sgUniAll);
    validSG = sgUniAll(sgCounts >= 20);
    
    fprintf('    PD genetic subgroups (n>=20): %s\n', strjoin(validSG, ', '));
    if numel(validSG) >= 2
        isPDsub = ismember(sgClean, validSG);
        blSG = bl(isPDsub, :);
        sgLabels = sgClean(isPDsub);
        fig14 = figure('Units','centimeters','Position',[2 2 17.4 10],'Color',COL.bg);
        
        % (a) Apathy/ICB prevalence by genetic subgroup
        ax1 = subplot(2,3,1); set(ax1,'Color',COL.bg); hold on;
        sgUni = validSG; nSG = numel(sgUni);
        prevApSG = nan(nSG,1); prevIcbSG = nan(nSG,1); prevBothSG = nan(nSG,1); nSGcounts = nan(nSG,1);
        for si = 1:nSG
            idx = sgLabels == sgUni(si) & ~isnan(blSG.group4);
            n = sum(idx); nSGcounts(si) = n;
            if n > 0
                prevApSG(si)   = 100*sum((blSG.group4==2|blSG.group4==4) & idx)/n;
                prevIcbSG(si)  = 100*sum((blSG.group4==3|blSG.group4==4) & idx)/n;
                prevBothSG(si) = 100*sum(blSG.group4==4 & idx)/n;
            end
        end
        x = 1:nSG; bw = 0.25;
        bar(x-bw, prevApSG, bw, 'FaceColor',COL.apathy,'EdgeColor','none','FaceAlpha',0.8);
        bar(x,    prevIcbSG, bw, 'FaceColor',COL.icb,   'EdgeColor','none','FaceAlpha',0.8);
        bar(x+bw, prevBothSG, bw, 'FaceColor',COL.overlap,'EdgeColor','none','FaceAlpha',0.8);
        
        set(ax1,'XTick',x,'XTickLabel',sgUni,'FontSize',6.5,'Box','off', ...
            'YGrid','on','GridColor',COL.grid,'TickDir','out');
        xtickangle(30); ylabel('Prevalence (%)','FontSize',7);
        lg=legend({'Apathy','ICB','Both'},'Location','northwest','Box','off','FontSize',5.5);
        lg.ItemTokenSize=[6 6];
        title('(a) Prevalence by Subgroup','FontSize',8,'FontWeight','bold','Color',COL.accent);
        for si = 1:nSG
            text(x(si), -4, sprintf('n=%d',nSGcounts(si)), 'HorizontalAlignment','center','FontSize',5,'Color',[0.5 0.5 0.5]);
        end
        
        % (b-f) Key variables by genetic subgroup
        sgTestVars  = {'apathy_composite','cfa_g','updrs3_score','stai_trait','depression_factor'};
        sgTestNames = {'(b) Apathy','(c) CFA Global Cog','(d) Motor','(e) Anxiety','(f) Depression'};
        
        for vi = 1:numel(sgTestVars)
            ax = subplot(2,3,vi+1); set(ax,'Color',COL.bg); hold on;
            vals = safeCol(blSG, sgTestVars{vi});
            for si = 1:nSG
                gv = vals(sgLabels==sgUni(si)); gv = gv(~isnan(gv));
                if isempty(gv), continue; end
                jit = 0.15*(rand(numel(gv),1)-0.5);
                sgCols = [COL.neither; COL.apathy; COL.icb; COL.overlap; 0.5 0.5 0.5];
                scatter(si+jit, gv, 10, sgCols(min(si,size(sgCols,1)),:), 'filled','MarkerFaceAlpha',0.25);
                q = quantile(gv, [0.25 0.5 0.75]);
                plot([si-0.12 si+0.12], [q(2) q(2)], '-k', 'LineWidth',1.5);
                if q(3) > q(1)
                    rectangle('Position',[si-0.15, q(1), 0.30, q(3)-q(1)], ...
                        'EdgeColor','k','LineWidth',0.8,'Curvature',0.1);
                end
            end
            if vi == 2, yline(ax, 0, '--', 'Color', [0.7 0.7 0.7], 'LineWidth', 0.5); end % 0-line for CFA
            
            ok = ~isnan(vals) & sgLabels ~= "";
            if sum(ok) > 10
                p = kruskalwallis(vals(ok), sgLabels(ok), 'off');
                text(0.95,0.95, sprintf('p=%.3f%s',p,pToStars(p)), 'Units','normalized', ...
                    'HorizontalAlignment','right','VerticalAlignment','top','FontSize',6,'Color',COL.txt);
            end
            set(ax,'XTick',1:nSG,'XTickLabel',sgUni,'FontSize',6,'Box','off', ...
                'YGrid','on','GridColor',COL.grid,'TickDir','out');
            xtickangle(30);
            title(sgTestNames{vi},'FontSize',7,'FontWeight','bold','Color',COL.accent);
            hold off;
        end
        sgtitle('Genetic Subgroup Analysis (PPMI)','FontSize',10,'FontWeight','bold','Color',COL.accent);
        exportFigure(fig14, fullfile(outputDir, 'Fig14_Genetic_Subgroups'));
        fprintf('    -> Fig14_Genetic_Subgroups\n');
        
        fprintf('\n    Genetic subgroup comparisons (KW):\n');
        for vi = 1:numel(sgTestVars)
            vals = safeCol(blSG, sgTestVars{vi});
            ok = ~isnan(vals);
            if sum(ok) > 10
                p = kruskalwallis(vals(ok), sgLabels(ok), 'off');
                fprintf('      %-18s  p = %.4f %s\n', sgTestNames{vi}(5:end), p, pToStars(p));
            end
        end
    end
end
fprintf('\n');

%% ===  20. FIG 15: LONGITUDINAL COGNITION BY GENETICS  ===================
fprintf('[20] Figure 15: Longitudinal cognition stratified by genetics...\n');

% Ensure we pull from Tpd_lme because it contains the CFA factor scores
if exist('Tpd_lme', 'var')
    T_genetics = Tpd_lme;
else
    warning('Tpd_lme not found. Defaulting to T. cfa_g may be missing.');
    T_genetics = T(isPD, :);
end

apoe4_all_raw = safeCol(T_genetics, 'APOE_e4');
apoe4_all = double(apoe4_all_raw >= 1); 
apoe4_all(isnan(apoe4_all_raw)) = NaN;
hasAPOElong = sum(~isnan(apoe4_all)) > 50;

if hasAPOElong
    fig15 = figure('Units','centimeters','Position',[2 2 17.4 6.5],'Color',COL.bg);
    years = 0:5; % Ensure years is defined
    
    % (a) CFA Cognition by APOE e4 x Apathy group
    ax1 = subplot(1,3,1); set(ax1,'Color',COL.bg); hold on;
    apComp_all = safeCol(T_genetics, 'apathy_composite');
    strataColors = [COL.neither; COL.apathy; COL.icb; COL.overlap];
    strataNames = {'e4-/NoAp','e4-/Ap','e4+/NoAp','e4+/Ap'};
    strataDef = {apoe4_all==0 & apComp_all<THRESH.apathy_comp, apoe4_all==0 & apComp_all>=THRESH.apathy_comp, ...
                 apoe4_all==1 & apComp_all<THRESH.apathy_comp, apoe4_all==1 & apComp_all>=THRESH.apathy_comp};
    phS = gobjects(4,1);
    for si = 1:4
        sIdx = strataDef{si};
        tSub = T_genetics(sIdx, :);
        if height(tSub) < 10, continue; end
        nYrs = numel(years);
        mY = nan(nYrs,1); sY = nan(nYrs,1); nY = nan(nYrs,1);
        for yi = 1:nYrs
            yIdx = safeCol(tSub,'YEAR') == years(yi);
            vals = safeCol(tSub,'cfa_g'); vals = vals(yIdx); vals = vals(~isnan(vals));
            nY(yi) = numel(vals);
            if numel(vals) >= 5
                mY(yi) = mean(vals); sY(yi) = std(vals)/sqrt(numel(vals));
            end
        end
        vYrs = ~isnan(mY);
        phS(si) = plotTrajectory(ax1, years, mY, sY, nY, vYrs, strataColors(si,:), 'none');
    end
    yline(ax1, 0, '--', 'Color', [0.7 0.7 0.7], 'LineWidth', 0.5);
    xlabel('Year','FontSize',7); ylabel('Global Cognition (cfa_g)','FontSize',7,'Interpreter','none');
    title('(a) CFA Cog: APOE \epsilon4 x Apathy','FontSize',8,'FontWeight','bold','Color',COL.accent);
    vph = phS(isgraphics(phS)); vn = strataNames(isgraphics(phS));
    if ~isempty(vph), lg=legend(vph,vn,'Location','southwest','FontSize',5,'Box','off'); lg.ItemTokenSize=[6 6]; end
    formatAxis(ax1, COL); xlim([-0.3 max(years)+0.3]);
    
    % (b) CFA Cognition by genetic subgroup
    ax2 = subplot(1,3,2); set(ax2,'Color',COL.bg); hold on;
    if ismember('subgroup', T_genetics.Properties.VariableNames)
        sgAll = string(T_genetics.subgroup);
        sgAllClean = sgAll;
        sgAllClean(contains(sgAllClean,'LRRK2') & ~contains(sgAllClean,'GBA')) = "LRRK2";
        sgAllClean(contains(sgAllClean,'GBA')   & ~contains(sgAllClean,'LRRK2')) = "GBA";
        sgAllClean(contains(sgAllClean,'LRRK2') &  contains(sgAllClean,'GBA'))  = "LRRK2+GBA";
        sgAllClean(sgAllClean == "SNCA")  = "SNCA";
        sgAllClean(sgAllClean == "PRKN" | sgAllClean == "PINK1") = "PRKN/PINK1";
        
        sgAllUni = unique(sgAllClean);
        sgAllCnts = arrayfun(@(s) sum(sgAllClean==s), sgAllUni);
        pdSubs = sgAllUni(sgAllCnts >= 50);
        sgCols = lines(numel(pdSubs));
        phSG = gobjects(numel(pdSubs),1);
        for si = 1:numel(pdSubs)
            tSub = T_genetics(sgAllClean == pdSubs(si), :);
            nYrs = numel(years);
            mY = nan(nYrs,1); sY = nan(nYrs,1); nY = nan(nYrs,1);
            for yi = 1:nYrs
                yIdx = safeCol(tSub,'YEAR') == years(yi);
                vals = safeCol(tSub,'cfa_g'); vals = vals(yIdx); vals = vals(~isnan(vals));
                nY(yi) = numel(vals);
                if numel(vals) >= 5
                    mY(yi) = mean(vals); sY(yi) = std(vals)/sqrt(numel(vals));
                end
            end
            vYrs = ~isnan(mY);
            phSG(si) = plotTrajectory(ax2, years, mY, sY, nY, vYrs, sgCols(si,:), 'none');
        end
        yline(ax2, 0, '--', 'Color', [0.7 0.7 0.7], 'LineWidth', 0.5);
        vph = phSG(isgraphics(phSG)); vn = pdSubs(isgraphics(phSG));
        if ~isempty(vph), lg=legend(vph,vn,'Location','southwest','FontSize',5,'Box','off'); lg.ItemTokenSize=[6 6]; end
    end
    xlabel('Year','FontSize',7); ylabel('Global Cognition (cfa_g)','FontSize',7,'Interpreter','none');
    title('(b) CFA Cog by Genetic Subgroup','FontSize',8,'FontWeight','bold','Color',COL.accent);
    formatAxis(ax2, COL); xlim([-0.3 max(years)+0.3]);
    
    % (c) Motor by genetic subgroup
    ax3 = subplot(1,3,3); set(ax3,'Color',COL.bg); hold on;
    if ismember('subgroup', T_genetics.Properties.VariableNames)
        phSG3 = gobjects(numel(pdSubs),1);
        for si = 1:numel(pdSubs)
            tSub = T_genetics(sgAllClean == pdSubs(si), :);
            nYrs = numel(years);
            mY = nan(nYrs,1); sY = nan(nYrs,1); nY = nan(nYrs,1);
            for yi = 1:nYrs
                yIdx = safeCol(tSub,'YEAR') == years(yi);
                vals = safeCol(tSub,'updrs3_score'); vals = vals(yIdx); vals = vals(~isnan(vals));
                nY(yi) = numel(vals);
                if numel(vals) >= 5
                    mY(yi) = mean(vals); sY(yi) = std(vals)/sqrt(numel(vals));
                end
            end
            vYrs = ~isnan(mY);
            phSG3(si) = plotTrajectory(ax3, years, mY, sY, nY, vYrs, sgCols(si,:), 'none');
        end
        vph = phSG3(isgraphics(phSG3)); vn = pdSubs(isgraphics(phSG3));
        if ~isempty(vph), lg=legend(vph,vn,'Location','northwest','FontSize',5,'Box','off'); lg.ItemTokenSize=[6 6]; end
    end
    xlabel('Year','FontSize',7); ylabel('UPDRS-III Score','FontSize',7);
    title('(c) Motor by Genetic Subgroup','FontSize',8,'FontWeight','bold','Color',COL.accent);
    formatAxis(ax3, COL); xlim([-0.3 max(years)+0.3]);
    
    sgtitle('Longitudinal Outcomes by Genetics','FontSize',10,'FontWeight','bold','Color',COL.accent);
    exportFigure(fig15, fullfile(outputDir, 'Fig15_Genetics_Longitudinal'));
    fprintf('    -> Fig15_Genetics_Longitudinal\n\n');
end

%%
%% ===  20. FIG 15: LONGITUDINAL COGNITION & MOTOR BY GENETICS  ===================
fprintf('[20] Figure 15: Longitudinal cognition and motor outcomes stratified by genetics...\n');

% Ensure we pull from Tpd_lme because it contains the CFA factor scores
if exist('Tpd_lme', 'var')
    T_genetics = Tpd_lme;
else
    warning('Tpd_lme not found. Defaulting to T. cfa_g may be missing.');
    T_genetics = T(isPD, :);
end

% --- Prep Variables ---
years = 0:5; 
nYrs = numel(years);
apoe4_all_raw = safeCol(T_genetics, 'APOE_e4');
apoe4_all = double(apoe4_all_raw >= 1); 
apoe4_all(isnan(apoe4_all_raw)) = NaN;

apComp_all = safeCol(T_genetics, 'apathy_composite');
strataColors = [COL.neither; COL.apathy; COL.icb; COL.overlap];

hasAPOElong = sum(~isnan(apoe4_all)) > 50;
if hasAPOElong
    % Expanded figure size to accommodate a 3x3 grid
    fig15 = figure('Units','centimeters','Position',[2 2 22 18],'Color',COL.bg);
    
    % --- Prep Genetic Subgroups & Filter 'Normosomic' ---
    sgAllClean = strings(height(T_genetics), 1);
    pdSubs = [];
    if ismember('subgroup', T_genetics.Properties.VariableNames)
        sgAll = string(T_genetics.subgroup);
        sgAllClean = sgAll;
        sgAllClean(contains(sgAllClean,'LRRK2') & ~contains(sgAllClean,'GBA')) = "LRRK2";
        sgAllClean(contains(sgAllClean,'GBA')   & ~contains(sgAllClean,'LRRK2')) = "GBA";
        sgAllClean(contains(sgAllClean,'LRRK2') &  contains(sgAllClean,'GBA'))  = "LRRK2+GBA";
        sgAllClean(sgAllClean == "SNCA")  = "SNCA";
        sgAllClean(sgAllClean == "PRKN" | sgAllClean == "PINK1") = "PRKN/PINK1";
        
        sgAllUni = unique(sgAllClean);
        sgAllCnts = arrayfun(@(s) sum(sgAllClean==s), sgAllUni);
        pdSubs = sgAllUni(sgAllCnts >= 50);
        
        % Remove normosomic / normosmic group from visualization list
        removeIdx = contains(lower(pdSubs), 'normosomic') | contains(lower(pdSubs), 'normosmic');
        pdSubs(removeIdx) = [];
    end
    
    % Define specific gene positives/negatives for apathy split panels
    apoe_pos = apoe4_all == 1;
    apoe_neg = apoe4_all == 0;
    
    % INCLUDE GBA REGARDLESS OF LRRK2 STATUS
    gba_pos = contains(sgAllClean, "GBA"); 
    gba_neg = ~gba_pos & sgAllClean ~= ""; 
    
    snca_pos = contains(sgAllClean, "SNCA");
    snca_neg = ~snca_pos & sgAllClean ~= "";
    
    % Configuration for automated 2x3 panel generation (Gene x Apathy splits)
    % {ax_idx, data_col, gene_pos, gene_neg, gene_name, y_label, title_str}
    configs = {
        {1, 'cfa_g', apoe_pos, apoe_neg, 'e4', 'Global Cognition (cfa_g)', '(a) CFA Cog: APOE \epsilon4 x Apathy'},
        {2, 'updrs3_score', apoe_pos, apoe_neg, 'e4', 'UPDRS-III Score', '(b) Motor: APOE \epsilon4 x Apathy'},
        {4, 'cfa_g', gba_pos, gba_neg, 'GBA', 'Global Cognition (cfa_g)', '(d) CFA Cog: GBA x Apathy'},
        {5, 'updrs3_score', gba_pos, gba_neg, 'GBA', 'UPDRS-III Score', '(e) Motor: GBA x Apathy'},
        {7, 'cfa_g', snca_pos, snca_neg, 'SNCA', 'Global Cognition (cfa_g)', '(g) CFA Cog: SNCA x Apathy'},
        {8, 'updrs3_score', snca_pos, snca_neg, 'SNCA', 'UPDRS-III Score', '(h) Motor: SNCA x Apathy'}
    };

    % Loop through configurations to build the Gene x Apathy panels
    for ci = 1:length(configs)
        c = configs{ci};
        ax_idx = c{1}; valCol = c{2}; gPos = c{3}; gNeg = c{4}; gName = c{5}; yLab = c{6}; tStr = c{7};
        
        ax = subplot(3,3,ax_idx); set(ax,'Color',COL.bg); hold on;
        
        sNames = {[gName '-/NoAp'], [gName '-/Ap'], [gName '+/NoAp'], [gName '+/Ap']};
        sDef = {gNeg & apComp_all<THRESH.apathy_comp, gNeg & apComp_all>=THRESH.apathy_comp, ...
                gPos & apComp_all<THRESH.apathy_comp, gPos & apComp_all>=THRESH.apathy_comp};
            
        phS = gobjects(4,1);
        for si = 1:4
            sIdx = sDef{si};
            tSub = T_genetics(sIdx, :);
            if height(tSub) < 10, continue; end
            
            mY = nan(nYrs,1); sY = nan(nYrs,1); nY = nan(nYrs,1);
            for yi = 1:nYrs
                yIdx = safeCol(tSub,'YEAR') == years(yi);
                vals = safeCol(tSub, valCol); vals = vals(yIdx); vals = vals(~isnan(vals));
                nY(yi) = numel(vals);
                if numel(vals) >= 5
                    mY(yi) = mean(vals); sY(yi) = std(vals)/sqrt(numel(vals));
                end
            end
            vYrs = ~isnan(mY);
            if any(vYrs)
                phS(si) = plotTrajectory(ax, years, mY, sY, nY, vYrs, strataColors(si,:), 'none');
            end
        end
        
        yline(ax, 0, '--', 'Color', [0.7 0.7 0.7], 'LineWidth', 0.5);
        xlabel(ax, 'Year','FontSize',7); ylabel(ax, yLab,'FontSize',7,'Interpreter','none');
        title(ax, tStr,'FontSize',8,'FontWeight','bold','Color',COL.accent);
        
        vph = phS(isgraphics(phS)); vn = sNames(isgraphics(phS));
        if ~isempty(vph)
            lg=legend(ax, vph, vn, 'Location','southwest','FontSize',5,'Box','off'); lg.ItemTokenSize=[6 6]; 
        end
        formatAxis(ax, COL); xlim(ax, [-0.3 max(years)+0.3]);
    end

    % --- General Genetic Subgroup Panels (Column 3) ---
    if ~isempty(pdSubs)
        sgCols = lines(numel(pdSubs));
        
        % (c) CFA Cognition by genetic subgroup
        ax3 = subplot(3,3,3); set(ax3,'Color',COL.bg); hold on;
        phSG1 = gobjects(numel(pdSubs),1);
        for si = 1:numel(pdSubs)
            tSub = T_genetics(sgAllClean == pdSubs(si), :);
            mY = nan(nYrs,1); sY = nan(nYrs,1); nY = nan(nYrs,1);
            for yi = 1:nYrs
                yIdx = safeCol(tSub,'YEAR') == years(yi);
                vals = safeCol(tSub,'cfa_g'); vals = vals(yIdx); vals = vals(~isnan(vals));
                nY(yi) = numel(vals);
                if numel(vals) >= 5
                    mY(yi) = mean(vals); sY(yi) = std(vals)/sqrt(numel(vals));
                end
            end
            vYrs = ~isnan(mY);
            if any(vYrs)
                phSG1(si) = plotTrajectory(ax3, years, mY, sY, nY, vYrs, sgCols(si,:), 'none');
            end
        end
        yline(ax3, 0, '--', 'Color', [0.7 0.7 0.7], 'LineWidth', 0.5);
        xlabel(ax3, 'Year','FontSize',7); ylabel(ax3, 'Global Cognition (cfa_g)','FontSize',7,'Interpreter','none');
        title(ax3, '(c) CFA Cog by Subgroup','FontSize',8,'FontWeight','bold','Color',COL.accent);
        vph = phSG1(isgraphics(phSG1)); vn = pdSubs(isgraphics(phSG1));
        if ~isempty(vph), lg=legend(ax3, vph,vn,'Location','southwest','FontSize',5,'Box','off'); lg.ItemTokenSize=[6 6]; end
        formatAxis(ax3, COL); xlim(ax3, [-0.3 max(years)+0.3]);

        % (f) Motor by genetic subgroup
        ax6 = subplot(3,3,6); set(ax6,'Color',COL.bg); hold on;
        phSG2 = gobjects(numel(pdSubs),1);
        for si = 1:numel(pdSubs)
            tSub = T_genetics(sgAllClean == pdSubs(si), :);
            mY = nan(nYrs,1); sY = nan(nYrs,1); nY = nan(nYrs,1);
            for yi = 1:nYrs
                yIdx = safeCol(tSub,'YEAR') == years(yi);
                vals = safeCol(tSub,'updrs3_score'); vals = vals(yIdx); vals = vals(~isnan(vals));
                nY(yi) = numel(vals);
                if numel(vals) >= 5
                    mY(yi) = mean(vals); sY(yi) = std(vals)/sqrt(numel(vals));
                end
            end
            vYrs = ~isnan(mY);
            if any(vYrs)
                phSG2(si) = plotTrajectory(ax6, years, mY, sY, nY, vYrs, sgCols(si,:), 'none');
            end
        end
        xlabel(ax6, 'Year','FontSize',7); ylabel(ax6, 'UPDRS-III Score','FontSize',7);
        title(ax6, '(f) Motor by Subgroup','FontSize',8,'FontWeight','bold','Color',COL.accent);
        vph = phSG2(isgraphics(phSG2)); vn = pdSubs(isgraphics(phSG2));
        if ~isempty(vph), lg=legend(ax6, vph,vn,'Location','northwest','FontSize',5,'Box','off'); lg.ItemTokenSize=[6 6]; end
        formatAxis(ax6, COL); xlim(ax6, [-0.3 max(years)+0.3]);
        
        % (i) APOE Cross-Prevalence by Subgroup
        ax9 = subplot(3,3,9); set(ax9,'Color',COL.bg); hold on;
        prev_apoe = nan(numel(pdSubs), 1);
        for si = 1:numel(pdSubs)
            idx = (sgAllClean == pdSubs(si));
            a_stat = apoe4_all(idx);
            a_stat = a_stat(~isnan(a_stat));
            if ~isempty(a_stat)
                prev_apoe(si) = mean(a_stat) * 100;
            end
        end
        b = bar(ax9, categorical(pdSubs), prev_apoe, 'FaceColor', COL.accent, 'EdgeColor', 'none', 'BarWidth', 0.6);
        ylabel(ax9, '% APOE \epsilon4 Positive', 'FontSize', 7);
        title(ax9, '(i) APOE \epsilon4 Prevalence', 'FontSize', 8, 'FontWeight', 'bold', 'Color', COL.accent);
        xtickangle(ax9, 30);
        formatAxis(ax9, COL);
    end
    
    sgtitle('Longitudinal Outcomes by Genetics','FontSize',10,'FontWeight','bold','Color',COL.accent);
    exportFigure(fig15, fullfile(outputDir, 'Fig15_Genetics_Longitudinal_Expanded'));
    fprintf('    -> Fig15_Genetics_Longitudinal_Expanded\n\n');
end
%% ═══════════════════════════════════════════════════════════════════════
%% [8] SUBCORTICAL VOLUMES — Group difference (HC-derived ICV correction)
%% ═══════════════════════════════════════════════════════════════════════
%% [8] SUBCORTICAL VOLUMES — Group difference (HC-derived ICV correction)
%% ═══════════════════════════════════════════════════════════════════════
fprintf('\n[8] Subcortical volumes by group (ANCOVA & HC-Residuals)...\n');

% --- Safety: Initialize counters and check group variable ---
if ~exist('figN', 'var'), figN = 0; end
if ~exist('tblN', 'var'), tblN = 0; end
bl.grp = bl.group4;
nGrp = 4; % explicitly define number of groups for plotting loops
hasFS = false;

if isfile(volFile)
    fsV = readtable(volFile,'VariableNamingRule','preserve');
    fsV.PATNO = string(fsV.PATNO);
    
    % Safely filter for Baseline using string comparison
    fsV_bl = fsV(strcmp(string(fsV.EVENT_ID), "BL"), :);
    [~,uI] = unique(fsV_bl.PATNO,'first'); 
    fsV_bl = fsV_bl(uI,:);
    
    % ROI Definitions
    roiDef = {
        'Hippocampus',     'Left_Hippocampus',      'Right_Hippocampus';
        'Amygdala',        'Left_Amygdala',         'Right_Amygdala';
        'Caudate',         'Left_Caudate',          'Right_Caudate';
        'Putamen',         'Left_Putamen',          'Right_Putamen';
        'Thalamus',        'Left_Thalamus',         'Right_Thalamus';
        'Accumbens',       'Left_Accumbens_area',   'Right_Accumbens_area';
        'Pallidum',        'Left_Pallidum',         'Right_Pallidum'};
        
    allCols = unique([roiDef(:,2); roiDef(:,3); {'EstimatedTotalIntraCranialVol'}]);
    allCols = allCols(allCols ~= "");
    fsCols = intersect(allCols, fsV_bl.Properties.VariableNames);
    fsCols = reshape(fsCols, 1, []);  
    fsKey  = fsV_bl(:, [{'PATNO'}, fsCols]);
    
    % Merge into PD baseline table
    bl = outerjoin(bl, fsKey, 'Keys','PATNO', 'MergeKeys',true, 'Type','left');
    
    % Extract Healthy Controls (Cohort 2) at Baseline safely
    hc_bl = T(strcmp(string(T.COHORT), "2") & strcmp(string(T.EVENT_ID), "BL"), :);
    hc_bl = outerjoin(hc_bl, fsKey, 'Keys','PATNO', 'MergeKeys',true, 'Type','inner');
    
    % Compute bilateral averages for PD and HC
    for ri = 1:size(roiDef,1)
        lbl = roiDef{ri,1}; lCol = roiDef{ri,2}; rCol = roiDef{ri,3};
        if isempty(rCol) || strcmp(rCol,'')
            bl.(['fs_' lbl]) = sCol(bl, lCol);
            hc_bl.(['fs_' lbl]) = sCol(hc_bl, lCol);
        else
            bl.(['fs_' lbl]) = (sCol(bl,lCol) + sCol(bl,rCol)) / 2;
            hc_bl.(['fs_' lbl]) = (sCol(hc_bl,lCol) + sCol(hc_bl,rCol)) / 2;
        end
    end
    
    bl.eTIV = sCol(bl,'EstimatedTotalIntraCranialVol');
    hc_bl.eTIV = sCol(hc_bl,'EstimatedTotalIntraCranialVol');
    
    % --- DERIVE NORMATIVE ICV BETAS FROM HEALTHY CONTROLS ---
    hc_mean_icv = nanmean(hc_bl.eTIV);
    hc_beta = nan(size(roiDef,1), 1);
    
    for ri = 1:size(roiDef,1)
        lbl = roiDef{ri,1};
        v = hc_bl.(['fs_' lbl]);
        ok = ~isnan(v) & ~isnan(hc_bl.eTIV);
        if sum(ok) > 15
            b = polyfit(hc_bl.eTIV(ok), v(ok), 1);
            hc_beta(ri) = b(1);
        else
            hc_beta(ri) = 0; % Fallback if HC data is highly missing
        end
    end
    
    hasFS = sum(~isnan(bl.eTIV) & ~isnan(bl.grp)) > 30;
    fprintf('  PD with FreeSurfer + group: %d\n', sum(~isnan(bl.eTIV) & ~isnan(bl.grp)));
    fprintf('  Healthy Controls used for ICV normative modeling: %d\n', sum(~isnan(hc_bl.eTIV)));
end

if hasFS
    nROI = size(roiDef,1);
    blAgeScan = sCol(bl,'age_datscan')  ;
    blduration = sCol(bl,'duration_yrs')  ;
    blAge = sCol(bl,'age_at_visit') - blAgeScan;
    blSex = sCol(bl,'SEX');
    eTIV  = bl.eTIV;
    
    % Setup Table
    tblN = tblN+1;
    resCell = cell(nROI+1, 7);
    resCell(1,:) = {'Structure','F_group','p_group','eta2_group','Cohen_d_1v4','mean_diff_1v4','p_FDR'};
    pRaw = nan(nROI,1);
    dVals = nan(nROI,1);
    
    for ri = 1:nROI
        lbl = roiDef{ri,1};
        rawV = bl.(['fs_' lbl]);
        okAll = ~isnan(rawV) & ~isnan(eTIV) & ~isnan(blAge) & ~isnan(blSex) & ~isnan(bl.grp)  & ~isnan(blduration)  & ~isnan(blAgeScan);
        resCell{ri+1,1} = lbl;
        
        if sum(okAll) < 50
            resCell(ri+1,2:7) = {NaN};
            continue;
        end
        
        % --- 1. APPLY NORMATIVE HC RESIDUAL CORRECTION FIRST ---
        % Preserves absolute scale (mm^3) while removing natural head-size scaling
        corrV = rawV - hc_beta(ri) .* (eTIV - hc_mean_icv);
    
        % --- 2. SETUP TABLE WITH CORRECTED VOLUME ---
        grpCat = categorical(bl.grp(okAll));
        
        % ICV is dropped as a covariate because it is already mathematically removed from vol_adj
        tFit = table(corrV(okAll), grpCat, blAge(okAll), blSex(okAll), blduration(okAll),blAgeScan(okAll), ...
            'VariableNames',{'vol_adj','grp','age','sex','blduration','blAgeScan'});
        
        try
            % Primary ANCOVA on adjusted volumes
            mdl = fitlm(tFit,'vol_adj ~ grp + age+ blAgeScan  + blduration + sex');
            aovT = anova(mdl,'component',3);
            
            rowIdx = find(contains(aovT.Properties.RowNames, 'grp'));
            if ~isempty(rowIdx)
                resCell{ri+1,2} = aovT.F(rowIdx);
                resCell{ri+1,3} = aovT.pValue(rowIdx);
                pRaw(ri) = aovT.pValue(rowIdx);
                resCell{ri+1,4} = aovT.SumSq(rowIdx) / (aovT.SumSq(rowIdx) + aovT.SumSq(end));
            end
        catch ME
            fprintf('  Warning (%s): %s\n', lbl, ME.message);
        end
        
        % Cohen's d: Grp 4 vs Grp 1 (Using the cleanly corrected volumes)
        v1 = corrV(bl.grp==1 & okAll); 
        v4 = corrV(bl.grp==4 & okAll); 
        if numel(v1)>5 && numel(v4)>5
            poolSD = sqrt(((numel(v1)-1)*var(v1)+(numel(v4)-1)*var(v4))/(numel(v1)+numel(v4)-2));
            dVals(ri) = (mean(v4)-mean(v1)) / poolSD;
            resCell{ri+1,5} = dVals(ri);
            resCell{ri+1,6} = mean(v4)-mean(v1);
        end
    end
    
    % FDR correction
    pAdj = fdr_bh(pRaw(~isnan(pRaw)));
    validIdx = find(~isnan(pRaw));
    for i = 1:numel(validIdx)
        resCell{validIdx(i)+1,7} = pAdj(i);
    end
    
    resT = cell2table(resCell(2:end,:),'VariableNames',resCell(1,:));
    writetable(resT, fullfile(outputDir, sprintf('Table%d_Subcortical_ANCOVA.csv',tblN)));
    disp(resT);
    
    % ── Figure: Subcortical Volumes (Normative Corrected) ──
    figN = figN+1;
    fig = figure('Units','centimeters','Position',[2 2 17.4 12],'Color',COL.bg);
    for ri = 1:nROI
        ax = subplot(2, ceil(nROI/2), ri); hold(ax,'on');
        rawV = bl.(['fs_' roiDef{ri,1}]);
        okAll = ~isnan(rawV) & ~isnan(eTIV) & ~isnan(bl.grp);
        if sum(okAll)<30, continue; end
        
        % Apply HC correction for plotting
        corrV = rawV - hc_beta(ri) .* (eTIV - hc_mean_icv);
        
        for g = 1:nGrp
            gv = corrV(bl.grp==g & okAll);
            scatter(ax, g + 0.1*(rand(size(gv))-0.5), gv, 4, COL.groups(g,:), 'filled', 'MarkerFaceAlpha', 0.3);
            errorbar(ax, g, mean(gv), std(gv)/sqrt(numel(gv)), 'o', 'MarkerFaceColor', COL.groups(g,:), 'Color', 'k');
        end
        title(ax, roiDef{ri,1}, 'FontSize', 8, 'FontWeight','bold');
        ylabel(ax, 'Adj. Vol (mm^3)', 'FontSize', 6);
        set(ax, 'XTick', 1:4, 'XTickLabel', groupLabelsShort, 'FontSize', 6, 'TickDir', 'out', 'Box', 'off');
        xtickangle(ax, 30);
        fmtAx(ax, COL);
    end
    % saveFig(fig, fullfile(outputDir, sprintf('Fig%d_SubcortVol_Main', figN)));
    
    %% ═══════════════════════════════════════════════════════════════════
    %% [9] SLIDING MEAN — top 3 structures
    %% ═══════════════════════════════════════════════════════════════════
    fprintf('\n[9] Sliding-mean window visualization...\n');
    [~, rankI] = sort(abs(dVals),'descend');
    top3 = rankI(1:min(3, sum(~isnan(dVals))));
    
    blDur = sCol(bl,'duration_yrs');
    winHalf = 5; step = 0.5;
    figN = figN+1;
    fig = figure('Units','centimeters','Position',[2 2 17.4 6],'Color',COL.bg);
    
    for ti = 1:numel(top3)
        ri = top3(ti);
        ax = subplot(1,3,ti); hold(ax,'on');
        rawV = bl.(['fs_' roiDef{ri,1}]);
        okAll = ~isnan(rawV) & ~isnan(eTIV) & ~isnan(blAge) & ~isnan(blDur) & ~isnan(bl.grp);
        
        % Apply HC correction for the sliding window
        corrV = rawV - hc_beta(ri) .* (eTIV - hc_mean_icv);
        
        ageGrid = (floor(min(blAge(okAll))):step:ceil(max(blAge(okAll))))';
        for g = 1:nGrp
            idx = bl.grp==g & okAll;
            if sum(idx) < 5, continue; end
            age_g = blAge(idx); vol_g = corrV(idx);
            rM = nan(size(ageGrid));
            for ai = 1:numel(ageGrid)
                inW = age_g >= (ageGrid(ai)-winHalf) & age_g <= (ageGrid(ai)+winHalf);
                if sum(inW) >= 5, rM(ai) = mean(vol_g(inW)); end
            end
            plot(ax, ageGrid, rM, '-', 'Color', COL.groups(g,:), 'LineWidth', 1.5);
        end
        title(ax, roiDef{ri,1}, 'FontSize', 8, 'FontWeight','bold');
        xlabel(ax, 'Age (years)', 'FontSize', 7);
        ylabel(ax, 'Adj. Vol (mm^3)', 'FontSize', 7);
        fmtAx(ax, COL);
    end
    % saveFig(fig, fullfile(outputDir, sprintf('Fig%d_SlidingMean', figN)));
end


%%
%% ═══════════════════════════════════════════════════════════════════════
%% [10] DATSCAN SBR ANALYSIS (BASELINE & LONGITUDINAL)
%% ═══════════════════════════════════════════════════════════════════════
fprintf('\n[10] Analyzing Striatal DaTscan SBR...\n');

% --- 1. PREPARE DATSCAN VARIABLES ---
if ismember('MIA_CAUDATE_BILAT', bl.Properties.VariableNames)
    bl.dat_caudate  = safeCol(bl, 'MIA_CAUDATE_BILAT');
    bl.dat_putamen  = safeCol(bl, 'MIA_PUTAMEN_BILAT');
    bl.dat_striatum = safeCol(bl, 'MIA_STRIATUM_BILAT');
else
    warning('DaTscan variables (MIA_CAUDATE_BILAT, etc.) not found in baseline table.');
end

if ismember('MIA_CAUDATE_BILAT', Tpd.Properties.VariableNames)
    Tpd.dat_caudate  = safeCol(Tpd, 'MIA_CAUDATE_BILAT');
    Tpd.dat_putamen  = safeCol(Tpd, 'MIA_PUTAMEN_BILAT');
    Tpd.dat_striatum = safeCol(Tpd, 'MIA_STRIATUM_BILAT');
else
    warning('DaTscan variables not found in longitudinal table.');
end

blAgeScan = safeCol(bl, 'age_datscan');
if all(isnan(blAgeScan)), blAgeScan = safeCol(bl, 'age_at_visit'); end
blSex = safeCol(bl, 'SEX');
blGrp = bl.group4;

% --- 2. BASELINE ANCOVA ---
fprintf('    Baseline DaTscan ANCOVA (adjusted for age and sex):\n');
datRegions = {'dat_caudate', 'dat_putamen'};
datLabels  = {'Caudate SBR', 'Putamen SBR'};

for ri = 1:numel(datRegions)
    rawV = bl.(datRegions{ri});
    okAll = ~isnan(rawV) & ~isnan(blAgeScan) & ~isnan(blSex) & ~isnan(blGrp);
    
    if sum(okAll) > 30
        tFit = table(rawV(okAll), categorical(blGrp(okAll)), blAgeScan(okAll), blSex(okAll), ...
            'VariableNames', {'sbr', 'grp', 'age', 'sex'});
        
        mdl = fitlm(tFit, 'sbr ~ grp + age + sex');
        aovT = anova(mdl, 'component', 3);
        
        rowIdx = find(contains(aovT.Properties.RowNames, 'grp'));
        if ~isempty(rowIdx)
            pval = aovT.pValue(rowIdx);
            Fval = aovT.F(rowIdx);
            fprintf('      %-12s: F = %5.2f, p = %.4f %s\n', datLabels{ri}, Fval, pval, pToStars(pval));
        end
    end
end

% --- 3. LONGITUDINAL LME ---
fprintf('\n    Longitudinal DaTscan LME (adjusted for age and sex):\n');
for ri = 1:numel(datRegions)
    vName = datRegions{ri};
    okLME = ~isnan(Tpd.(vName)) & ~isnan(Tpd.group4) & ~isnan(Tpd.YEAR) & ~isnan(safeCol(Tpd,'SEX'));
    lmeD = Tpd(okLME, :);
    
    if height(lmeD) > 100
        lmeD.grpCat = categorical(lmeD.group4);
        formula = sprintf('%s ~ grpCat*YEAR + SEX + age_at_visit + (1 + YEAR | PATNO)', vName);
        
        try
            mdl_long = fitlme(lmeD, formula, 'FitMethod', 'REML');
            aov_long = anova(mdl_long);
            
            rowIdx = find(contains(aov_long.Term, 'grpCat:YEAR'));
            if ~isempty(rowIdx)
                pval = aov_long.pValue(rowIdx);
                Fval = aov_long.FStat(rowIdx);
                fprintf('      %-12s (Group x Time): F = %5.2f, p = %.4f %s\n', datLabels{ri}, Fval, pval, pToStars(pval));
            end
        catch ME
            fprintf('      %-12s: LME failed (%s)\n', datLabels{ri}, ME.message);
        end
    end
end

% --- 4. VISUALIZATION (4-PANEL FIGURE) ---
fprintf('\n    Generating DaTscan Visualizations...\n');
figDat = figure('Units','centimeters','Position',[2 2 20 14],'Color',COL.bg); % Taller to fit brackets
nGrp = 4;
years = 0:5;

for ri = 1:numel(datRegions)
    % =========================================================
    % TOP ROW: Baseline Bar + Scatter + Error + Stats Brackets
    % =========================================================
    axBase = subplot(2, 2, ri); hold(axBase, 'on');
    rawV = bl.(datRegions{ri});
    okAll = ~isnan(rawV) & ~isnan(blGrp);
    
    % Collect group data to find bounds for plotting stats
    grpData = cell(1, nGrp);
    mY = zeros(1, nGrp);
    seY = zeros(1, nGrp);
    for g = 1:nGrp
        gv = rawV(blGrp == g & okAll);
        grpData{g} = gv;
        if ~isempty(gv)
            mY(g) = mean(gv);
            seY(g) = std(gv)/sqrt(numel(gv));
        end
    end
    
    % Draw Bars, Jittered Scatter, and Error Bars
    for g = 1:nGrp
        if isempty(grpData{g}), continue; end
        
        % 1. Bar Plot
        b = bar(axBase, g, mY(g), 0.65);
        b.FaceColor = COL.groups(g,:);
        b.EdgeColor = 'none';
        b.FaceAlpha = 0.4; % Semi-transparent bar
        
        % 2. Jittered Scatter
        jit = 0.15 * (rand(size(grpData{g})) - 0.5);
        scatter(axBase, g + jit, grpData{g}, 8, COL.groups(g,:)*0.8, 'filled', 'MarkerFaceAlpha', 0.4);
        
        % 3. Error Bar (SEM)
        errorbar(axBase, g, mY(g), seY(g), 'k', 'LineWidth', 1.2, 'CapSize', 3);
    end
    
    % Statistical Significance Lines (Pairwise comparisons)
    pairs = [1 2; 1 3; 1 4; 2 3; 2 4; 3 4]; % All possible post-hoc pairs
    y_max = max(rawV(okAll));
    y_range = y_max - min(rawV(okAll));
    y_step = y_range * 0.08; % Step height for overlapping brackets
    curr_y = y_max + (y_range * 0.05); % Starting height for brackets
    
    sig_found = false;
    for p_idx = 1:size(pairs,1)
        g1 = pairs(p_idx, 1);
        g2 = pairs(p_idx, 2);
        
        if numel(grpData{g1}) > 5 && numel(grpData{g2}) > 5
            % Unpaired t-test for visual post-hoc representation
            [~, p_val] = ttest2(grpData{g1}, grpData{g2});
            
            if p_val < 0.05
                sig_found = true;
                % Draw Bracket
                plot(axBase, [g1, g2], [curr_y, curr_y], '-k', 'LineWidth', 1);
                plot(axBase, [g1, g1], [curr_y - y_step*0.2, curr_y], '-k', 'LineWidth', 1);
                plot(axBase, [g2, g2], [curr_y - y_step*0.2, curr_y], '-k', 'LineWidth', 1);
                
                % Determine Stars
                if p_val < 0.001, p_str = '***';
                elseif p_val < 0.01, p_str = '**';
                else, p_str = '*';
                end
                
                % Add Text
                text(axBase, mean([g1, g2]), curr_y + y_step*0.1, p_str, ...
                    'HorizontalAlignment', 'center', 'FontSize', 10, 'FontWeight', 'bold');
                
                % Move up for the next potential bracket
                curr_y = curr_y + y_step * 1.3; 
            end
        end
    end
    
    % Adjust Y-limits to fit the brackets dynamically
    if sig_found
        ylim(axBase, [0, curr_y]);
    else
        ylim(axBase, [0, y_max + y_range*0.1]);
    end
    
    title(axBase, sprintf('Baseline %s', datLabels{ri}), 'FontSize', 9, 'FontWeight', 'bold', 'Color', COL.accent);
    ylabel(axBase, 'SBR', 'FontSize', 8);
    set(axBase, 'XTick', 1:4, 'XTickLabel', groupLabelsShort, 'FontSize', 7, 'TickDir', 'out', 'Box', 'off');
    xtickangle(axBase, 30);
    formatAxis(axBase, COL);

    % =========================================================
    % BOTTOM ROW: Longitudinal Trajectories
    % =========================================================
    axLong = subplot(2, 2, ri + 2); hold(axLong, 'on');
    ph = gobjects(nGrp, 1);
    
    for g = 1:nGrp
        mY_long = nan(numel(years), 1);
        sY_long = nan(numel(years), 1);
        nY_long = nan(numel(years), 1);
        
        for yi = 1:numel(years)
            yIdx = Tpd.group4 == g & Tpd.YEAR == years(yi);
            vals = Tpd.(datRegions{ri})(yIdx);
            vals = vals(~isnan(vals));
            
            nY_long(yi) = numel(vals);
            if nY_long(yi) >= 5
                mY_long(yi) = mean(vals);
                sY_long(yi) = std(vals)/sqrt(nY_long(yi));
            end
        end
        
        vYrs = nY_long >= 5;
        if any(vYrs)
            ph(g) = plotTrajectory(axLong, years, mY_long, sY_long, nY_long, vYrs, COL.groups(g,:), 'endpoints');
        end
    end
    
    title(axLong, sprintf('Longitudinal %s', datLabels{ri}), 'FontSize', 9, 'FontWeight', 'bold', 'Color', COL.accent);
    xlabel(axLong, 'Year', 'FontSize', 8);
    ylabel(axLong, 'SBR', 'FontSize', 8);
    set(axLong, 'XTick', years, 'FontSize', 7, 'TickDir', 'out', 'Box', 'off');
    xlim(axLong, [-0.3, max(years)+0.3]);
    formatAxis(axLong, COL);
    
    if ri == 1
        valid_ph = ph(isgraphics(ph));
        valid_lg = groupLabelsShort(isgraphics(ph));
        lg = legend(axLong, valid_ph, valid_lg, 'Location', 'southwest', 'FontSize', 6, 'Box', 'off');
        lg.ItemTokenSize = [10 10];
    end
end

sgtitle('Striatal Dopaminergic Denervation (DaTscan) by Phenotype', 'FontSize', 12, 'FontWeight', 'bold', 'Color', COL.accent);
exportFigure(figDat, fullfile(outputDir, 'Fig_DaTscan_Analysis'));
fprintf('    -> Fig_DaTscan_Analysis saved.\n');

%%
%% ═══════════════════════════════════════════════════════════════════════
%% [11] COGNITION VS DATSCAN ASSOCIATIONS (CONTINUOUS INTERACTIONS)
%% ═══════════════════════════════════════════════════════════════════════
fprintf('\n[11] Analyzing Associations between DaTscan SBR and CFA Cognition...\n');

% --- 1. PREPARE VARIABLES & COVARIATES ---
cfa_targets = {'bifactor_g', 'bifactor_att_exec'};
for cv = 1:numel(cfa_targets)
    if ~ismember(cfa_targets{cv}, bl.Properties.VariableNames)
        if exist('Tpd_lme', 'var') && ismember(cfa_targets{cv}, Tpd_lme.Properties.VariableNames)
            bl_patno = str2double(string(bl.PATNO));
            lme_patno = str2double(string(Tpd_lme.PATNO));
            lme_year = double(safeCol(Tpd_lme, 'YEAR'));
            
            bl.(cfa_targets{cv}) = nan(height(bl), 1);
            for i = 1:height(bl)
                idx = (lme_patno == bl_patno(i)) & (lme_year == 0);
                if any(idx)
                    bl.(cfa_targets{cv})(i) = Tpd_lme.(cfa_targets{cv})(find(idx, 1));
                end
            end
        end
    end
end

cogVars  = {'bifactor_g', 'bifactor_att_exec'};
cogNames = {'Global Cognition (g)', 'Attention/Executive'};
datVars  = {'dat_caudate', 'dat_putamen'};
datNames = {'Caudate SBR', 'Putamen SBR'};

% Continuous Behavioral Predictors
blApathy = safeCol(bl, 'apathy_composite');
blIcb    = safeCol(bl, 'quip');

% Core Demographics & Clinical Covariates
blAge  = safeCol(bl, 'age_at_visit');
blSex  = safeCol(bl, 'SEX');
blGrp  = bl.group4; % Kept ONLY for visualization colors
blEduc = safeCol(bl, 'EDUCYRS');
blDep  = safeCol(bl, 'depression_factor');
blHY   = safeCol(bl, 'hy');
blLedd = safeCol(bl, 'LEDD');
blDur  = safeCol(bl, 'duration_yrs');
if all(isnan(blDur)), blDur = safeCol(bl, 'duration_yrs'); end

% --- 2. STATISTICAL INTERACTION MODELS (CONTINUOUS) ---
fprintf('    Linear Models (Cognition ~ DaT * Apathy * ICB + Covariates):\n');
for cv = 1:numel(cogVars)
    for dv = 1:numel(datVars)
        yVar = safeCol(bl, cogVars{cv});
        xVar = safeCol(bl, datVars{dv});
        
        okM  = ~isnan(yVar) & ~isnan(xVar) & ~isnan(blApathy) & ~isnan(blIcb) & ...
               ~isnan(blAge) & ~isnan(blSex) & ~isnan(blEduc) & ~isnan(blDep) & ...
               ~isnan(blDur) & ~isnan(blHY) & ~isnan(blLedd);
        
        if sum(okM) > 50
            tFit = table(yVar(okM), xVar(okM), blApathy(okM), blIcb(okM), blAge(okM), ...
                         blSex(okM), blEduc(okM), blDep(okM), blDur(okM), blHY(okM), blLedd(okM), ...
                'VariableNames', {'cog', 'dat', 'apathy', 'icb', 'age', 'sex', 'educ', 'dep', 'duration', 'hy', 'ledd'});
            
            % Continuous 3-way interaction model
            mdl = fitlm(tFit, 'dat ~ apathy*icb + age + sex + educ + dep + duration + hy + ledd');
            aov = anova(mdl);
            
            % Safely extract interaction terms (MATLAB sometimes flips the order like 'apathy:dat')
            rNames = aov.Properties.RowNames;
            idx_dat        = find(strcmp(rNames, 'dat'));
            idx_dat_ap     = find(contains(rNames, 'dat') & contains(rNames, 'apathy') & ~contains(rNames, 'icb'));
            idx_dat_icb    = find(contains(rNames, 'dat') & contains(rNames, 'icb') & ~contains(rNames, 'apathy'));
            idx_dat_ap_icb = find(contains(rNames, 'dat') & contains(rNames, 'apathy') & contains(rNames, 'icb'));
            
            fprintf('\n      %s vs %s (Adjusted n = %d):\n', cogNames{cv}, datNames{dv}, sum(okM));
            if ~isempty(idx_dat)
                fprintf('        Main Effect (DaT):          F = %5.2f, p = %.4f %s\n', ...
                    aov.F(idx_dat), aov.pValue(idx_dat), pToStars(aov.pValue(idx_dat)));
            end
            if ~isempty(idx_dat_ap)
                fprintf('        Interaction (DaT*Apathy):   F = %5.2f, p = %.4f %s\n', ...
                    aov.F(idx_dat_ap), aov.pValue(idx_dat_ap), pToStars(aov.pValue(idx_dat_ap)));
            end
            if ~isempty(idx_dat_icb)
                fprintf('        Interaction (DaT*ICB):      F = %5.2f, p = %.4f %s\n', ...
                    aov.F(idx_dat_icb), aov.pValue(idx_dat_icb), pToStars(aov.pValue(idx_dat_icb)));
            end
            if ~isempty(idx_dat_ap_icb)
                fprintf('        3-Way (DaT*Apathy*ICB):     F = %5.2f, p = %.4f %s\n', ...
                    aov.F(idx_dat_ap_icb), aov.pValue(idx_dat_ap_icb), pToStars(aov.pValue(idx_dat_ap_icb)));
            end
        end
    end
end

% --- 3. VISUALIZATION (2x2 SCATTER GRID) ---
fprintf('\n    Generating DaTscan vs Cognition Scatter Plots...\n');
figDatCog = figure('Units','centimeters','Position',[2 2 20 16],'Color',COL.bg);
plotIdx = 1;

for cv = 1:numel(cogVars)
    for dv = 1:numel(datVars)
        ax = subplot(2, 2, plotIdx); hold(ax, 'on');
        
        yVar = safeCol(bl, cogVars{cv});
        xVar = safeCol(bl, datVars{dv});
        % We use blGrp here just to color the scatter points intuitively
        okAll = ~isnan(yVar) & ~isnan(xVar) & ~isnan(blGrp); 
        
        ph = gobjects(4, 1); 
        
        for g = 1:4
            idx = okAll & (blGrp == g);
            if sum(idx) < 5, continue; end
            
            xG = xVar(idx);
            yG = yVar(idx);
            
            % Scatter plot
            ph(g) = scatter(ax, xG, yG, 18, COL.groups(g,:), 'filled', 'MarkerFaceAlpha', 0.4, 'MarkerEdgeColor', 'none');
            
            % Group-specific regression line for visual aid
            b = polyfit(xG, yG, 1);
            xLine = [min(xG), max(xG)];
            yLine = polyval(b, xLine);
            plot(ax, xLine, yLine, '-', 'Color', COL.groups(g,:), 'LineWidth', 2.5);
            
            % Calculate within-group Pearson correlation for console output
            [r, p] = corrcoef(xG, yG);
            if p(1,2) < 0.05
                fprintf('        Visual Group %d (%s): r = %+.3f, p = %.4f %s\n', ...
                    g, groupLabelsShort{g}, r(1,2), p(1,2), pToStars(p(1,2)));
            end
        end
        
        % Formatting
        xlabel(ax, datNames{dv}, 'FontSize', 8, 'FontWeight', 'bold');
        ylabel(ax, cogNames{cv}, 'FontSize', 8, 'FontWeight', 'bold', 'Interpreter', 'none');
        title(ax, sprintf('%s vs %s', cogNames{cv}, datNames{dv}), 'FontSize', 9, 'Color', COL.accent);
        
        % Plot reference line for Cognitive Mean
        yline(ax, 0, '--', 'Color', [0.7 0.7 0.7], 'LineWidth', 1);
        
        formatAxis(ax, COL);
        
        % Add legend to the first plot
        if plotIdx == 1
            valid_ph = ph(isgraphics(ph));
            valid_lg = groupLabelsShort(isgraphics(ph));
            lg = legend(ax, valid_ph, valid_lg, 'Location', 'southeast', 'FontSize', 6, 'Box', 'off');
            lg.ItemTokenSize = [12 12];
        end
        
        plotIdx = plotIdx + 1;
    end
end

sgtitle('Associations: DaTscan SBR vs CFA Cognition (Continuous Models)', 'FontSize', 12, 'FontWeight', 'bold', 'Color', COL.accent);
exportFigure(figDatCog, fullfile(outputDir, 'Fig_DaTscan_vs_Cognition'));
fprintf('    -> Fig_DaTscan_vs_Cognition saved.\n');

%%
%% ═══════════════════════════════════════════════════════════════════════
%% [24] GBA MUTATION ANALYSIS (BASELINE & LONGITUDINAL LME)
%% ═══════════════════════════════════════════════════════════════════════
fprintf('\n[24] Analyzing GBA mutation carrier status...\n');

% --- 1. DEFINE GBA STATUS (BASELINE) ---
% Bypass safeCol to prevent it from accidentally turning text into NaNs
if ismember('subgroup', bl.Properties.VariableNames)
    raw_sg = bl.subgroup;
elseif ismember('SUBGROUP', bl.Properties.VariableNames)
    raw_sg = bl.SUBGROUP;
else
    error('The variable "subgroup" was not found in the baseline table.');
end

% Convert robustly to string (handles cell arrays, categoricals, and char arrays)
sg_str_bl = string(raw_sg);
gba_bl = nan(height(bl), 1);

% Identify GBA carriers (1)
gba_bl(contains(sg_str_bl, 'GBA', 'IgnoreCase', true)) = 1;

% Identify valid non-carriers (0)
valid_text = (sg_str_bl ~= "") & ~ismissing(sg_str_bl) & (sg_str_bl ~= "NaN");
gba_bl(~contains(sg_str_bl, 'GBA', 'IgnoreCase', true) & valid_text) = 0;

hasGBA = sum(~isnan(gba_bl)) > 20;

if hasGBA
    n_gba1 = sum(gba_bl==1); 
    n_gba0 = sum(gba_bl==0);
    fprintf('    Extraction Check: Found %d valid text entries in subgroup.\n', sum(valid_text));
    fprintf('    GBA Carriers: %d  |  Non-carriers: %d\n', n_gba1, n_gba0);

    % =====================================================================
    % PART A: BASELINE CROSS-SECTIONAL (Prevalence & Group Differences)
    % =====================================================================
    figGBA = figure('Units','centimeters','Position',[2 2 17.4 7],'Color',COL.bg);
    
    % (a) Group prevalence by GBA status
    ax1 = subplot(1,3,1); set(ax1,'Color',COL.bg); hold on;
    for gs = 0:1
        gsIdx = gba_bl == gs & ~isnan(bl.group4);
        nGS = sum(gsIdx);
        pct = arrayfun(@(g) 100*sum(bl.group4==g & gsIdx)/nGS, 1:4);
        xpos = (1:4) + (gs-0.5)*0.35;
        
        if gs == 0
            bColor = [0.6 0.6 0.6]; % Gray for non-carriers
        else
            bColor = COL.apathy;    % Highlight color for carriers
        end
        
        bar(xpos, pct, 0.30, 'FaceColor', bColor, 'EdgeColor', 'none', 'FaceAlpha', 0.8);
    end
    set(ax1,'XTick',1:4,'XTickLabel',groupLabelsShort,'FontSize',7,'Box','off', ...
        'YGrid','on','GridColor',COL.grid,'TickDir','out');
    xtickangle(30); ylabel('Prevalence (%)','FontSize',8);
    lg = legend({sprintf('GBA- (n=%d)',n_gba0), sprintf('GBA+ (n=%d)',n_gba1)}, ...
        'Location','northeast','Box','off','FontSize',6);
    lg.ItemTokenSize=[8 8];
    title('(a) Group by GBA Status','FontSize',8,'FontWeight','bold','Color',COL.accent);
    
    % (b) Apathy composite by GBA x group
    ax2 = subplot(1,3,2); set(ax2,'Color',COL.bg); hold on;
    for gs = 0:1
        for g = 1:4
            idx = gba_bl==gs & bl.group4==g;
            vals = safeCol(bl,'apathy_composite'); vals = vals(idx); vals = vals(~isnan(vals));
            if isempty(vals), continue; end
            xp = g + (gs-0.5)*0.3;
            jit = 0.06*(rand(numel(vals),1)-0.5);
            
            if gs == 0
                mcol = COL.groups(g,:) * 0.6 + 0.4; % Lighter/washed out for non-carriers
            else
                mcol = COL.groups(g,:); % Full color for carriers
            end
            
            scatter(xp+jit, vals, 10, mcol, 'filled', 'MarkerFaceAlpha',0.3);
            plot([xp-0.08 xp+0.08], [median(vals) median(vals)], '-', 'Color',mcol,'LineWidth',2);
        end
    end
    set(ax2,'XTick',1:4,'XTickLabel',groupLabelsShort,'FontSize',7,'Box','off', ...
        'YGrid','on','GridColor',COL.grid,'TickDir','out');
    xtickangle(30); ylabel('Apathy Composite','FontSize',8);
    title('(b) Apathy \times GBA','FontSize',8,'FontWeight','bold','Color',COL.accent);
    
    % (c) CFA Global Cognition by GBA x group
    ax3 = subplot(1,3,3); set(ax3,'Color',COL.bg); hold on;
    for gs = 0:1
        for g = 1:4
            idx = gba_bl==gs & bl.group4==g;
            vals = safeCol(bl,'bifactor_g'); vals = vals(idx); vals = vals(~isnan(vals));
            if isempty(vals), continue; end
            xp = g + (gs-0.5)*0.3;
            jit = 0.06*(rand(numel(vals),1)-0.5);
            
            if gs == 0
                mcol = COL.groups(g,:) * 0.6 + 0.4;
            else
                mcol = COL.groups(g,:);
            end
            
            scatter(xp+jit, vals, 10, mcol, 'filled', 'MarkerFaceAlpha',0.3);
            plot([xp-0.08 xp+0.08], [median(vals) median(vals)], '-', 'Color',mcol,'LineWidth',2);
        end
    end
    yline(ax3, 0, '--', 'Color', [0.7 0.7 0.7], 'LineWidth', 0.5);
    set(ax3,'XTick',1:4,'XTickLabel',groupLabelsShort,'FontSize',7,'Box','off', ...
        'YGrid','on','GridColor',COL.grid,'TickDir','out');
    xtickangle(30); ylabel('Global Cognition (bifactor_g)','FontSize',8,'Interpreter','none');
    title('(c) CFA Cognition \times GBA','FontSize',8,'FontWeight','bold','Color',COL.accent);
    
    sgtitle('GBA Mutation Carrier Status Analysis','FontSize',10,'FontWeight','bold','Color',COL.accent);
    exportFigure(figGBA, fullfile(outputDir, 'Fig_GBA_Baseline_Analysis'));
    fprintf('    -> Fig_GBA_Baseline_Analysis saved.\n');
    
    % Statistical tests (Mann-Whitney U)
    fprintf('\n    GBA effects (Mann-Whitney U):\n');
    gbaTestVars = {'apathy_composite','quip','bifactor_g','depression_factor','updrs3_score','stai_trait'};
    gbaTestNames = {'Apathy','QUIP','CFA Global Cog','Depression','Motor','Anxiety'};
    for vi = 1:numel(gbaTestVars)
        v1 = safeCol(bl, gbaTestVars{vi});
        g_neg = v1(gba_bl==0 & ~isnan(v1)); 
        g_pos = v1(gba_bl==1 & ~isnan(v1));
        if numel(g_neg)>5 && numel(g_pos)>5
            p = ranksum(g_neg, g_pos);
            fprintf('      %-15s  GBA−: %5.2f±%.2f  GBA+: %5.2f±%.2f  p=%.4f %s\n', ...
                gbaTestNames{vi}, mean(g_neg),std(g_neg), mean(g_pos),std(g_pos), p, pToStars(p));
        end
    end
    
    % Odds ratio: GBA -> apathy
    a11 = sum(gba_bl==1 & safeCol(bl,'apathy_composite')>=THRESH.apathy_comp);
    a10 = sum(gba_bl==1 & safeCol(bl,'apathy_composite')<THRESH.apathy_comp);
    a01 = sum(gba_bl==0 & safeCol(bl,'apathy_composite')>=THRESH.apathy_comp);
    a00 = sum(gba_bl==0 & safeCol(bl,'apathy_composite')<THRESH.apathy_comp);
    if all([a11 a10 a01 a00] > 0)
        OR = (a11*a00)/(a10*a01);
        logOR = log(OR); SE_logOR = sqrt(1/a11+1/a10+1/a01+1/a00);
        CI = exp([logOR-1.96*SE_logOR, logOR+1.96*SE_logOR]);
        fprintf('\n    GBA -> Apathy: OR = %.2f (95%% CI: %.2f-%.2f)\n', OR, CI(1), CI(2));
    end

    % =====================================================================
    % PART B: LONGITUDINAL LME MODELS (GBA * Apathy * ICB * Time)
    % =====================================================================
    fprintf('\n    ============================================================\n');
    fprintf('    LONGITUDINAL LME MODELS CONTROLLING FOR GBA STATUS\n');
    fprintf('    ============================================================\n');
    
    % Define GBA safely in the longitudinal table (bypassing safeCol)
    if ismember('subgroup', Tpd_lme.Properties.VariableNames)
        raw_sg_lme = Tpd_lme.subgroup;
    elseif ismember('SUBGROUP', Tpd_lme.Properties.VariableNames)
        raw_sg_lme = Tpd_lme.SUBGROUP;
    else
        % Fallback if somehow not in the longitudinal table
        raw_sg_lme = string(NaN(height(Tpd_lme), 1)); 
    end
    
    sg_str_lme = string(raw_sg_lme);
    lmeT.gba = nan(height(Tpd_lme), 1);
    lmeT.gba(contains(sg_str_lme, 'GBA', 'IgnoreCase', true)) = 1;
    valid_text_lme = (sg_str_lme ~= "") & ~ismissing(sg_str_lme) & (sg_str_lme ~= "NaN");
    lmeT.gba(~contains(sg_str_lme, 'GBA', 'IgnoreCase', true) & valid_text_lme) = 0;
    
    cogOutcomes = {'bifactor_g', 'bifactor_att_exec', 'bifactor_memory', 'bifactor_visuo'};
    cogLabels   = {'General Cognition (g)', 'Attention/Executive', 'Memory', 'Visuospatial'};
    allFE_gba = {};
    
    for ci = 1:numel(cogOutcomes)
        outVar = cogOutcomes{ci};
        outLabel = cogLabels{ci};
        
        okM = okBase & ~isnan(lmeT.(outVar)) & ~isnan(lmeT.gba);
        okAdj = okM & ~isnan(lmeT.educ) & ~isnan(lmeT.nmiss) & ...
                ~isnan(lmeT.duration) & ~isnan(lmeT.ledd) & ~isnan(lmeT.hy) & ...
                ~isundefined(lmeT.sex) & ~isnan(lmeT.depfactor);
        lmeD_adj = lmeT(okAdj, :);
        
        try
            % Adjusted Model (All interactions up to 3-way with GBA)
            formula_adj = sprintf('%s ~ apathy*icb*time + gba*apathy*time + gba*icb*time + gba*apathy*icb + sex + educ + nmiss + duration + ledd + hy + depfactor + (1 + time | PATNO)', outVar);
            mdl_adj = fitlme(lmeD_adj, formula_adj, 'FitMethod','REML', 'DummyVarCoding', 'effects');
            fe_adj = coeffToTable(mdl_adj);
            
            fprintf('\n    Model: %s\n', outLabel);
            for ri = 1:height(fe_adj)
                nm = fe_adj.Name{ri};
                if contains(nm,'time') || contains(nm,'apathy') || contains(nm,'icb') || contains(nm,'gba')
                    fprintf('      %-40s  b=%+7.3f  p=%.4f %s\n', nm, fe_adj.Estimate(ri), fe_adj.pValue(ri), pToStars(fe_adj.pValue(ri)));
                end
            end
            
            fe_adj.Domain = repmat({outLabel}, height(fe_adj), 1);
            allFE_gba{ci} = fe_adj;
            
        catch ME
            fprintf('    ERROR in LME for %s: %s\n', outLabel, ME.message);
        end
    end

    % --- Forest Plot for GBA LME Results ---
    fprintf('\n    Generating GBA LME Forest Plot...\n');
    try
        figGBA_LME = figure('Units','centimeters','Position',[2 2 17.4 16],'Color',COL.bg);
        domainCols = [COL.groups(1,:); COL.apathy; COL.icb; COL.overlap];
        
        for ci = 1:4
            ax = subplot(2, 2, ci); set(ax, 'Color', COL.bg, 'FontSize', 6.5); hold on;
            if ci <= numel(allFE_gba) && ~isempty(allFE_gba{ci})
                fe = allFE_gba{ci};
                keepIdx = find(contains(fe.Name,'apathy') | contains(fe.Name,'icb') | strcmp(fe.Name,'time') | contains(fe.Name,'gba'));
                
                nK = numel(keepIdx);
                colSig = domainCols(ci,:);
                colNS = [0.65 0.65 0.65];
                
                for ki = 1:nK
                    ri = keepIdx(ki);
                    est = fe.Estimate(ri);
                    ci_lo = est - 1.96*fe.SE(ri);
                    ci_hi = est + 1.96*fe.SE(ri);
                    p = fe.pValue(ri);
                    yp = nK - ki + 1;
                    
                    mc = colNS;
                    if p < 0.05, mc = colSig; end
                    if p < 0.001, mc = colSig * 0.8; end
                    
                    plot(ax, [ci_lo ci_hi], [yp yp], '-', 'Color', mc, 'LineWidth', 2);
                    plot(ax, est, yp, 'o', 'Color', mc, 'MarkerFaceColor', mc, 'MarkerSize', 5.5, 'LineWidth', 0.5);
                    
                    if p < 0.05
                        text(ax, ci_hi + 0.01 * abs(ci_hi - ci_lo + 0.1), yp, pToStars(p), 'FontSize', 7, 'Color', mc, 'VerticalAlignment', 'middle', 'FontWeight', 'bold');
                    end
                end
                
                xline(ax, 0, '-', 'Color', [0.75 0.75 0.75], 'LineWidth', 0.5);
                
                tickLabels = fe.Name(keepIdx);
                tickLabels = strrep(tickLabels, 'gba_1', 'GBA');
                tickLabels = strrep(tickLabels, 'gba', 'GBA');
                tickLabels = strrep(tickLabels, 'apathy_1', 'Apathy');
                tickLabels = strrep(tickLabels, 'apathy', 'Apathy');
                tickLabels = strrep(tickLabels, 'icb_1', 'ICB');
                tickLabels = strrep(tickLabels, 'icb', 'ICB');
                tickLabels = strrep(tickLabels, 'time', 'Time');
                tickLabels = strrep(tickLabels, ':', ' \times ');
                
                panelLetter = char('A' + ci - 1);
                set(ax, 'YTick', 1:nK, 'YTickLabel', flip(tickLabels), 'TickDir', 'out', 'Box', 'off', 'YLim', [0.3 nK+0.7], 'XGrid', 'on', 'GridColor', COL.grid, 'GridAlpha', 0.25);
                
                xl = xlim(ax); xlim(ax, [xl(1), xl(2) + diff(xl)*0.15]);
                xlabel('Estimate (95% CI)', 'FontSize', 7);
                title(sprintf('%s   %s', panelLetter, cogLabels{ci}), 'FontSize', 8.5, 'FontWeight', 'bold', 'Color', domainCols(ci,:), 'HorizontalAlignment', 'left', 'Units', 'normalized', 'Position', [0 1.06 0]);
                formatAxis(ax, COL);
            end
        end
        
        sgtitle('LME Fixed Effects by Domain (Adjusted for GBA Mutation Status)', 'FontSize', 10, 'FontWeight', 'bold', 'Color', COL.accent);
        exportFigure(figGBA_LME, fullfile(outputDir, 'Fig_LME_GBA_Domain_Forest'));
        fprintf('    -> Fig_LME_GBA_Domain_Forest saved.\n');
    catch ME
        fprintf('    GBA Forest plot skipped: %s\n', ME.message);
    end
end
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

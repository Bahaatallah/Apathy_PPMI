function Figure4_render(featDir, runDir, slopesCSV, outStem)
% Draws manuscript Figure 4 (three panels: a, differential FC x global cognitive slope; b, structural
% correlates; c, left hippocampal atrophy) from the feature folder featDir, the output folder runDir
% of Figure4.m and Figure4_Structural.m, and the participant-slope CSV slopesCSV; writes <outStem>.pdf,
% <outStem>.png and <outStem>_printed_values.csv (every printed annotation string).
% Every annotated value is read from the run tables; the plotted points and the bootstrap density are
% recomputed from the exported features and slopes. Before export the script checks that no text box
% overlaps another or leaves the figure, that every font prints at 6 pt or more at the 140 mm inclusion
% width, and that no participant identifier appears in any text.
if nargin < 4, error('featDir, runDir, slopesCSV, outStem'); end
COL.bg=[1 1 1]; COL.txt=[0.15 0.15 0.18]; COL.grid=[0.85 0.85 0.85];
cSig=[0.20 0.40 0.65]; cUnc=[0.45 0.60 0.80]; cNS=[0.78 0.78 0.78]; cBG=[0.97 0.97 0.97]; cVerm=[0.80 0.40 0.00];
FS=8;
FSann=6.5;      % annotation font
FStit=7;      % panel title font; fits inside the plot box
FSlet=10;       % panel letter

% ---------- data (unchanged)
F  = readtable(fullfile(featDir,'nbs_suprathreshold_fc.csv'));
S  = readtable(slopesCSV);
sid = string(regexprep(string(F.subject_id),'^sub-','')); ssid = string(regexprep(string(S.subject_id),'^sub-',''));
[tf,loc] = ismember(sid, ssid); slope = nan(height(F),1); slope(tf) = S.slope_g(loc(tf));
fcDiff = F.diff_apathy_gt_quip_neg_comp1_suprathresh_fc;
ap = F.apathy_composite; FD = F.mean_FD; TIV = F.TIV; age = F.age; educ = F.education_years;
sexS = upper(strtrim(string(F.sex))); sexN = double(startsWith(sexS,'M')); sexN(~(startsWith(sexS,'M')|startsWith(sexS,'F'))) = NaN;
sens = readtable(fullfile(runDir,'Table_Fig4_Sensitivity.csv')); sens = sens(strcmp(sens.Test,'FCxGlobal'),:);
med  = readtable(fullfile(runDir,'Table_Fig4_Mediation.csv'));
st   = readtable(fullfile(runDir,'Table_Structural_Correlates.csv'));
Tl   = readtable(fullfile(runDir,'Table_Hippocampus_Longitudinal_Data.csv'));
Th   = readtable(fullfile(runDir,'Table_Hippocampus_Slope_x_FC.csv'));   % NEW read (p_FDR for panel d)
% pFDR across the FWE-significant components of the run (as the submitted figure: BH over all components)
q = bh(sens.Spearman_p); kd = find(startsWith(sens.Component,'Diff(-)1: Ap-dom'),1);
rTab = sens.Spearman_r(kd); qTab = q(kd); pTab = sens.Spearman_p(kd);   % raw p for the panel a annotation

fig = figure('Units','centimeters','Position',[1 1 17.5 6.8],'Color',COL.bg,'Visible','off');
set(fig,'DefaultAxesFontName','Helvetica','DefaultTextFontName','Helvetica','DefaultAxesFontWeight','normal','DefaultTextFontWeight','normal','DefaultAxesTitleFontWeight','normal');
tl = tiledlayout(1,20,'TileSpacing','compact','Padding','compact'); %#ok<NASGU>
cmap = interp1(linspace(0,1,5),[0.92 0.92 0.94; 0.72 0.76 0.82; 0.48 0.56 0.68; 0.28 0.38 0.56; 0.14 0.22 0.42],linspace(0,1,256),'pchip');
ANN = {};   % every printed annotation: {panel, item, string}
LAY = {};   % annotations placed after layout: {ax, lines, objectsXY, radius_pt, tag}

% ---------- (a) pooled scatter, differential FC x cognitive slope
axA = nexttile(1,[1,6]); ax = axA; set(ax,'Color',COL.bg); hold(ax,'on');
Z = [FD TIV]; ok = ~isnan(fcDiff)&~isnan(slope)&all(~isnan(Z),2);
Zi = [ones(sum(ok),1) Z(ok,:)];
xR = fcDiff(ok) - Zi*(Zi\fcDiff(ok)); yR = slope(ok) - Zi*(Zi\slope(ok));
try, dens = mvksdensity([xR yR],[xR yR],'Bandwidth',[std(xR)*0.35, std(yR)*0.35]); catch, dens = ones(size(xR)); end
scatter(ax,xR,yR,30,dens,'filled','MarkerFaceAlpha',0.75,'MarkerEdgeColor','none'); colormap(ax,cmap);
b = robustfit(xR,yR); xF = linspace(min(xR),max(xR),100); plot(ax,xF,b(1)+b(2)*xF,'-','Color',cVerm,'LineWidth',2.2);
yline(ax,0,'--','Color',COL.grid,'LineWidth',0.5);
xlabel(ax,'Apathy > ICB FC (resid)','FontSize',FS); ylabel(ax,'Global Cog Slope (resid)','FontSize',FS);
formatAxis(ax,COL);
setDataLimits(ax, xR, yR);
ylim(ax,[-0.4 0.2]);                            % fixed y-range as in the submitted figure
nOutA = sum(yR < -0.4 | yR > 0.2);
fprintf('[chain] panel a: %d of %d residualised slopes lie outside the fixed y-range [-0.4, 0.2] (included in the statistics, not displayed)\n', nOutA, numel(yR));
sA1 = sprintf('\\rho = %+.3f; p %s', rTab, eqP(pTab));   % two shorter lines for the narrow panel
sA2 = sprintf('p_{\\fontsize{6}FDR} %s; n = %d', eqP(qTab), sum(ok));
ANN(end+1,:) = {'a','rho; p', sA1}; ANN(end+1,:) = {'a','p_FDR; n', sA2};
LAYA = {ax, {sA1, sA2}, [xR yR], sqrt(30)/2, 'a'};   % panel a: annotation at the bottom-left (fixed y-range; no top band)
panelTitle(ax, 'a', 'Cognitive decline', FStit, FSlet, COL);
[rChk,pChk] = partialcorr(fcDiff(ok),slope(ok),Z(ok,:),'Type','Spearman');

% (the mediation analysis is reported in the Supplement)
% ---------- (c) structural forest (M1, differential), as submitted
axC = nexttile(7,[1,8]); ax = axC; set(ax,'Color',COL.bg); hold(ax,'on');
m = startsWith(st.FC_Component,'Diff(-)1: Ap-dom') & strcmp(st.Model,'M1_Primary'); T = st(m,:);   % the NBS-defined differential only (at 343 a bridge feature also starts with Diff(-)1)
rR=T.rho; pR=T.p_unc; qR=T.p_FDR; nm=T.ROI; nR=numel(rR);
isF = ~isnan(qR)&qR<0.05; isU = ~isnan(pR)&pR<0.01&~isF; keep=find(isF|isU);
if numel(keep)<15, rem_=setdiff((1:nR)',keep); [~,o]=sort(abs(rR(rem_)),'descend'); keep=[keep; rem_(o(1:min(15-numel(keep),numel(rem_))))]; end
[~,sK]=sort(abs(rR(keep)),'ascend'); sh=keep(sK); nS=numel(sh);
for ri=1:nS
    oi=sh(ri); y=ri; r=rR(oi);
    if mod(ri,2)==0, fill(ax,[-0.35 0.35 0.35 -0.35],[y-0.45 y-0.45 y+0.45 y+0.45],cBG,'EdgeColor','none','FaceAlpha',0.5); end
    if isF(oi), plot(ax,[0 r],[y y],'-','Color',cSig,'LineWidth',1.8); plot(ax,r,y,'o','MarkerSize',7,'MarkerFaceColor',cSig,'MarkerEdgeColor','w'); t_=text(ax,r+sign(r)*0.028,y,sprintf('%.2f',r),'FontSize',6.5,'FontWeight','bold','Color',cSig,'VerticalAlignment','middle','HorizontalAlignment',ternary(r>=0,'left','right')); ANN(end+1,:) = {'c',['rho ' char(nm(oi))], t_.String}; %#ok<AGROW>
    elseif isU(oi), plot(ax,[0 r],[y y],'-','Color',cUnc,'LineWidth',1.2); plot(ax,r,y,'o','MarkerSize',6,'MarkerFaceColor','w','MarkerEdgeColor',cUnc,'LineWidth',1.4); t_=text(ax,r+sign(r)*0.028,y,sprintf('%.2f',r),'FontSize',6,'Color',[0.35 0.35 0.35],'VerticalAlignment','middle','HorizontalAlignment',ternary(r>=0,'left','right')); ANN(end+1,:) = {'c',['rho ' char(nm(oi))], t_.String}; %#ok<AGROW>
    else, plot(ax,[0 r],[y y],'-','Color',[0.85 0.85 0.85],'LineWidth',0.6); plot(ax,r,y,'o','MarkerSize',4,'MarkerFaceColor',cNS,'MarkerEdgeColor','none'); end
end
ytl = cellfun(@roiLabel, nm(sh), 'UniformOutput', false);
for ri=1:nS, if isU(sh(ri))||isF(sh(ri)), ytl{ri}=['\bf\it\color[rgb]{0.15 0.15 0.18}' ytl{ri}]; else, ytl{ri}=['\rm\color[rgb]{0.45 0.45 0.45}' ytl{ri}]; end, end
xline(ax,0,'-','Color',[0.6 0.6 0.6]); set(ax,'YTick',1:nS,'YTickLabel',ytl,'YLim',[0.4 nS+0.6],'XLim',[-0.35 0.35]);
xlabel(ax,'\rho_{\fontsize{6}Sp} (partial)','FontSize',FS);
formatAxis(ax,COL); set(ax,'FontSize',7); ax.YAxis.FontWeight='normal';
sC = sprintf('%d of %d ROIs shown', nS, nR);
ANN(end+1,:) = {'c','ROIs shown', sC};
panelTitle(ax, 'b', 'Structural correlates', FStit, FSlet, COL);

% ---------- (d) longitudinal left hippocampus, as submitted
axD = nexttile(15,[1,6]); ax = axD; set(ax,'Color',COL.bg); hold(ax,'on');
pat = double(str2double(sid)); [~,li] = ismember(pat, double(Tl.PATNO));
slL=nan(height(F),1); vB=nan(height(F),1); mk=li>0; slL(mk)=Tl.slope_left(li(mk)); vB(mk)=Tl.vBL_left(li(mk));
Zg=[age sexN FD educ vB]; okg=~isnan(slL)&~isnan(fcDiff)&all(~isnan(Zg),2);
Zr=zeros(sum(okg),size(Zg,2)); for cc=1:size(Zg,2), Zr(:,cc)=tiedrank(Zg(okg,cc)); end; Zr=[ones(sum(okg),1) Zr];
xr=tiedrank(fcDiff(okg)); yr=tiedrank(slL(okg)); xr=xr-Zr*(Zr\xr); yr=yr-Zr*(Zr\yr);
try, dens=mvksdensity([xr yr],[xr yr],'Bandwidth',[std(xr)*0.35, std(yr)*0.35]); catch, dens=ones(size(xr)); end
scatter(ax,xr,yr,28,dens,'filled','MarkerFaceAlpha',0.75,'MarkerEdgeColor','none'); colormap(ax,cmap);
b=robustfit(xr,yr); xF=linspace(min(xr),max(xr),100); plot(ax,xF,b(1)+b(2)*xF,'-','Color',cVerm,'LineWidth',1.8);
[rg,pg]=partialcorr(fcDiff(okg),slL(okg),Zg(okg,:),'Type','Spearman');
% p_FDR of the same test, from the run's own table (FDR within side, across components)
hr = Th(strcmp(Th.Model,'M1_Primary') & strcmp(Th.Hippocampus,'Left') & startsWith(Th.FC_Component,'Diff(-)1: Ap-dom'),:);
assert(height(hr)==1, 'panel d: expected one M1/Left/Diff(-)1 row in Table_Hippocampus_Slope_x_FC.csv');
assert(abs(hr.rho - rg) < 1e-9 && abs(hr.p_unc - pg) < 1e-9 && hr.n == sum(okg), ...
       'panel d: table rho/p/n differ from the plotted test');
qg = hr.p_FDR_within_side;
xlabel(ax,{'Apathy > ICB FC','(resid)'},'FontSize',FS); ylabel(ax,{'Left hippocampus \Delta/yr','(resid)'},'FontSize',FS);
formatAxis(ax,COL);
setDataLimits(ax, xr, yr);
sD1 = sprintf('\\rho = %+.3f; p %s', rg, eqP(pg));
sD2 = sprintf('p_{\\fontsize{6}FDR} %s; n = %d (left)', eqP(qg), sum(okg));
ANN(end+1,:) = {'c','rho; p', sD1}; ANN(end+1,:) = {'c','p_FDR; n', sD2};
LAY(end+1,:) = {ax, {sD1, sD2}, [xr yr], sqrt(28)/2, 'c'};
panelTitle(ax, 'c', 'Hippocampal atrophy', FStit, FSlet, COL);

fprintf('[chain] path b diff: annotated rho %+.4f pFDR %.4f (recomputed %+.4f p %.4g)\n', rTab, qTab, rChk, pChk);
fprintf('[chain] longitudinal left rho %+.4f p %.4g pFDR %.4g n %d\n', rg, pg, qg, sum(okg));
set(fig,'Renderer','painters');
% ---------- annotations after the tiled layout has settled: text stacked from
% the top-left corner, then the y-limit raised just enough that every plotted
% object (with its marker radius) sits below the annotation block
drawnow;
H = cell(size(LAY,1),1);
for k = 1:size(LAY,1)
    H{k} = annTop(LAY{k,1}, LAY{k,2}, FSann, COL);
end
drawnow;
for k = 1:size(LAY,1), fitTopBand(LAY{k,1}, H{k}, LAY{k,3}, LAY{k,4}); end
drawnow;
for k = 1:size(LAY,1), fitTopBand(LAY{k,1}, H{k}, LAY{k,3}, LAY{k,4}); end   % second pass after tick relabelling
drawnow;
HA = annBottom(LAYA{1}, LAYA{2}, FSann, COL);   % panel a annotation, bottom-left
drawnow;
placeLetters(fig);
% (c) note under the ROI-name column, outside the plotting area (was inside,
% over the bottom row): right-aligned just left of the plot box, below it
text(axC, -0.04, -0.035, sC, 'Units','normalized', 'FontSize', FSann, 'Color', [0.35 0.35 0.35], ...
     'HorizontalAlignment','right', 'VerticalAlignment','top');
drawnow;
% ---------- checks: every point inside its axes, no point/object under an annotation, no text collisions
for k = 1:size(LAY,1)
    runOverlapCheck(LAY{k,1}, H{k}, LAY{k,3}, LAY{k,4}, LAY{k,5});
end
runOverlapCheck(LAYA{1}, HA, LAYA{3}, LAYA{4}, LAYA{5});
checkTextOverlaps(fig);
titleFit(fig);
fontAudit(fig);
noIdentifiers(fig, F.subject_id);
exportgraphics(fig,[outStem '.pdf'],'ContentType','vector'); exportgraphics(fig,[outStem '.png'],'Resolution',600);
Tann = cell2table(ANN, 'VariableNames', {'panel','item','printed'});
Tann.plain = cellfun(@plainTeX, Tann.printed, 'UniformOutput', false);
writetable(Tann, [outStem '_printed_values.csv']);
fprintf('[chain] wrote %s.pdf/.png and %s_printed_values.csv\n', outStem, outStem);
end

% ======================================================================
function q = bh(p)
    n=numel(p); [ps,o]=sort(p(:)); q=zeros(n,1); prev=1;
    for k=n:-1:1, prev=min(prev, ps(k)*n/k); q(o(k))=prev; end
end
function out=ternary(c,a,b), if c,out=a;else,out=b;end, end
function s=fmtP(p), if isnan(p), s='NaN'; elseif p<0.001, s='< 0.001'; else, s=sprintf('%.3f',p); end, end
function s=eqP(p)
    % "= 0.035" or "< 0.001" (same rounding as fmtP of the 2026-09-13 version)
    s = fmtP(p); if ~startsWith(s,'<'), s = ['= ' s]; end
end
function formatAxis(ax,COL)
    set(ax,'FontSize',8,'Box','off','TickDir','out','XGrid','off','YGrid','off','GridColor',COL.grid,'GridAlpha',0.4,'FontWeight','bold');
    try, set(ax.XLabel,'FontWeight','bold'); set(ax.YLabel,'FontWeight','bold'); catch, end
end
function setDataLimits(ax, x, y)
    % limits from the plotted values, 4 % padding on each side (the top is
    % raised later by fitTopBand to make room for the annotation)
    rx = max(x) - min(x); ry = max(y) - min(y);
    xlim(ax, [min(x) - 0.04*rx, max(x) + 0.04*rx]);
    ylim(ax, [min(y) - 0.04*ry, max(y) + 0.04*ry]);
end
function fitTopBand(ax, hTxt, P, rPt)
    % raise YLim(2) so that the highest plotted object (plus its marker
    % radius and a 1.5 % gap) lies below the lowest annotation line
    e = get(hTxt(end), 'Extent');                     % normalized, lowest line
    fB = e(2) - 0.015;
    u = ax.Units; ax.Units = 'points'; pp = ax.Position; ax.Units = u;
    rN = rPt / pp(4);
    yl = ax.YLim; yTop = max(P(:,2));
    yl2 = yl(1) + (yTop - yl(1)) / (fB - rN);
    ax.YLim = [yl(1) max(yl2, yl(2))];
end
function h = annBottom(ax, lines, fs, COL)
    % annotation lines stacked upward from the bottom-left corner of the axes
    n = numel(lines); h = gobjects(n,1); y = 0.03;
    for i = n:-1:1
        h(i) = text(ax, 0.025, y, lines{i}, 'Units','normalized', 'FontSize', fs, ...
            'VerticalAlignment','bottom', 'HorizontalAlignment','left', 'Color', COL.txt, 'Margin', 0.5);
        e = get(h(i), 'Extent'); y = e(2) + e(4) + 0.004;
    end
end
function h = annTop(ax, lines, fs, COL)
    % annotation lines stacked from the top-left corner of the axes
    h = gobjects(numel(lines),1);
    y = 0.985;
    for i = 1:numel(lines)
        h(i) = text(ax, 0.025, y, lines{i}, 'Units','normalized', 'FontSize', fs, ...
            'VerticalAlignment','top', 'HorizontalAlignment','left', 'Color', COL.txt, 'Margin', 0.5);
        e = get(h(i), 'Extent'); y = e(2) - 0.004;
    end
end
function panelTitle(ax, letter, ttl, fsT, fsL, COL, sub)
    % title left-aligned on the plot box; bold panel letter to its left,
    % outside the plot box, on the same baseline
    if nargin < 7 || isempty(sub)
        title(ax, ttl, 'FontSize', fsT, 'FontWeight', 'normal', 'Color', COL.txt);
    else
        [t, s] = title(ax, ttl, sub, 'FontSize', fsT, 'FontWeight', 'normal', 'Color', COL.txt); %#ok<ASGLU>
        s.FontSize = 6.5; s.Color = [0.35 0.35 0.35]; s.FontWeight = 'normal';
    end
    ax.TitleHorizontalAlignment = 'left';
    ax.UserData = struct('letter', letter, 'fsL', fsL);
end
function runOverlapCheck(ax, hTxt, P, rPt, tag)
    % points (or object vertices) P [x y] in data units; rPt = marker radius
    % in points. Counts P outside the axes and P within any annotation box.
    xl = ax.XLim; yl = ax.YLim;
    inAx = P(:,1) >= xl(1) & P(:,1) <= xl(2) & P(:,2) >= yl(1) & P(:,2) <= yl(2);
    u = ax.Units; ax.Units = 'points'; pp = ax.Position; ax.Units = u;
    rx = rPt / pp(3) * diff(xl); ry = rPt / pp(4) * diff(yl);
    nHit = 0;
    for i = 1:numel(hTxt)
        e = get(hTxt(i), 'Extent');                         % normalized
        bx = xl(1) + [e(1), e(1)+e(3)] * diff(xl);
        by = yl(1) + [e(2), e(2)+e(4)] * diff(yl);
        hit = P(:,1) >= bx(1)-rx & P(:,1) <= bx(2)+rx & P(:,2) >= by(1)-ry & P(:,2) <= by(2)+ry;
        nHit = nHit + sum(hit);
    end
    fprintf('[check] panel %s: %d of %d plotted positions inside the axes; %d under an annotation\n', ...
        tag, sum(inAx), size(P,1), nHit);
end
function placeLetters(fig)
    % bold panel letter at the far left of each tile's y-axis decorations,
    % bottom-aligned with the (left-aligned) title
    ax = findall(fig, 'Type', 'axes');
    for i = 1:numel(ax)
        U = ax(i).UserData;
        if isstruct(U) && isfield(U, 'letter')
            ti = ax(i).TightInset; u = ax(i).Units; ax(i).Units = 'normalized'; p = ax(i).Position; ax(i).Units = u;
            t = ax(i).Title; tu = t.Units; t.Units = 'normalized'; tp = t.Position; te = t.Extent; t.Units = tu;
            xL = -(ti(1)/p(3)) * 0.98;                        % far left of the y-axis decorations
            text(ax(i), xL, te(2), U.letter, 'Units','normalized', 'FontSize', U.fsL, ...
                'FontWeight','bold', 'VerticalAlignment','bottom', 'HorizontalAlignment','left', ...
                'Color',[0.10 0.10 0.12], 'Tag','panelLetter'); %#ok<NASGU>
            if tp(1) < -0.005, warning('title of panel %s starts left of the plot box (%.4f)', U.letter, tp(1)); end
        end
    end
end
function checkTextOverlaps(fig)
    % every pair of titles, subtitles, letters, axis labels and in-axes
    % annotations must be disjoint and inside the figure (figure units)
    T = findall(fig, 'Type', 'text');
    B = zeros(numel(T), 4); lab = strings(numel(T),1); keepT = false(numel(T),1);
    for i = 1:numel(T)
        s = string(T(i).String); if isempty(s) || all(strlength(s)==0), continue; end
        a = ancestor(T(i), 'axes'); if isempty(a), continue; end
        u = T(i).Units; T(i).Units = 'normalized'; e = T(i).Extent; T(i).Units = u;
        ua = a.Units; a.Units = 'normalized'; p = a.Position; a.Units = ua;
        B(i,:) = [p(1) + e(1)*p(3), p(2) + e(2)*p(4), e(3)*p(3), e(4)*p(4)];
        lab(i) = strjoin(s, ' '); keepT(i) = true;
    end
    idx = find(keepT); nOv = 0;
    for i = 1:numel(idx)
        bi = B(idx(i),:);
        if bi(1) < -0.001 || bi(1)+bi(3) > 1.001 || bi(2) < -0.001 || bi(2)+bi(4) > 1.001
            fprintf('[check] OFF-FIGURE: "%s"\n', lab(idx(i)));
            nOv = nOv + 1;
        end
        for j = i+1:numel(idx)
            bj = B(idx(j),:);
            ox = min(bi(1)+bi(3), bj(1)+bj(3)) - max(bi(1), bj(1));
            oy = min(bi(2)+bi(4), bj(2)+bj(4)) - max(bi(2), bj(2));
            if ox > 0 && oy > 0
                nOv = nOv + 1;
                fprintf('[check] TEXT OVERLAP: "%s" x "%s"\n', lab(idx(i)), lab(idx(j)));
            end
        end
    end
    fprintf('[check] text boxes (titles, letters, labels, annotations): %d overlaps / off-figure\n', nOv);
end
function s = plainTeX(s)
    s = regexprep(char(s), '\\fontsize\{[^}]*\}', '');
    s = strrep(s, '\rho', 'rho'); s = strrep(s, '\Delta', 'Delta');
    s = regexprep(s, '_\{([^}]*)\}', '_$1');
end
function s = roiLabel(r)
    % readable names for the 42 bilaterally averaged FreeSurfer ROIs
    M = containers.Map( ...
        {'Accumbens','Amygdala','Brain Stem','Caudate','Cereb Cortex','Hippocampus','Pallidum','Putamen','Thalamus', ...
         'caudalanteriorcingulate','caudalmiddlefrontal','cuneus','entorhinal','frontalpole','fusiform', ...
         'inferiorparietal','inferiortemporal','insula','isthmuscingulate','lateraloccipital','lateralorbitofrontal', ...
         'lingual','medialorbitofrontal','middletemporal','paracentral','parahippocampal','parsopercularis', ...
         'parsorbitalis','parstriangularis','pericalcarine','postcentral','posteriorcingulate','precentral', ...
         'precuneus','rostralanteriorcingulate','rostralmiddlefrontal','superiorfrontal','superiorparietal', ...
         'superiortemporal','supramarginal','temporalpole','transversetemporal','bankssts'}, ...
        {'Accumbens','Amygdala','Brain stem','Caudate','Cerebellar cortex','Hippocampus','Pallidum','Putamen','Thalamus', ...
         'Caudal anterior cingulate','Caudal middle frontal','Cuneus','Entorhinal','Frontal pole','Fusiform', ...
         'Inferior parietal','Inferior temporal','Insula','Isthmus cingulate','Lateral occipital','Lateral orbitofrontal', ...
         'Lingual','Medial orbitofrontal','Middle temporal','Paracentral','Parahippocampal','Pars opercularis', ...
         'Pars orbitalis','Pars triangularis','Pericalcarine','Postcentral','Posterior cingulate','Precentral', ...
         'Precuneus','Rostral anterior cingulate','Rostral middle frontal','Superior frontal','Superior parietal', ...
         'Superior temporal','Supramarginal','Temporal pole','Transverse temporal','Banks of the STS'});
    r = char(r);
    if isKey(M, r), s = M(r); else, error('roiLabel: no readable name for "%s"', r); end
end
function fontAudit(fig)
    % smallest nominal font size of every text object (incl. \fontsize{}
    % overrides), tick labels and legends; flags TeX sub/superscripts left at
    % the interpreter default size (about 0.76 x the base size)
    fs = []; lab = strings(0,1); nDef = 0;
    T = findall(fig, 'Type', 'text');
    for i = 1:numel(T)
        s = strjoin(txtOf(T(i).String)', ' '); if strlength(strtrim(s)) == 0, continue; end
        fs(end+1) = T(i).FontSize; lab(end+1) = s; %#ok<AGROW>
        o = regexp(char(s), '\\fontsize\{([\d.]+)\}', 'tokens');
        for k = 1:numel(o), fs(end+1) = str2double(o{k}{1}); lab(end+1) = s; end %#ok<AGROW>
        if ~isempty(regexp(char(s), '[_^]\{(?!\\fontsize)', 'once')), nDef = nDef + 1;
            fprintf('[check] default-size sub/superscript: "%s"\n', s); end
    end
    A = findall(fig, 'Type', 'axes');
    for i = 1:numel(A), fs(end+1) = A(i).FontSize; lab(end+1) = "[tick labels]"; end %#ok<AGROW>
    L = findall(fig, 'Type', 'legend');
    for i = 1:numel(L), fs(end+1) = L(i).FontSize; lab(end+1) = "[legend]"; end %#ok<AGROW>
    [m, im] = min(fs);
    fprintf('[check] fonts: smallest nominal %.2f pt ("%s"); %d default-size sub/superscripts\n', m, lab(im), nDef);
end
function noIdentifiers(fig, idsIn)
    % abort before export if a participant identifier (with or without the
    % "sub-" prefix) appears in any text object, tick label or legend entry
    ids = string(idsIn); ids = unique([ids; regexprep(ids, '^sub-', '')]); ids = ids(strlength(ids) >= 4);
    txt = strings(0,1);
    T = findall(fig, 'Type', 'text'); for i = 1:numel(T), txt = [txt; txtOf(T(i).String)]; end %#ok<AGROW>
    A = findall(fig, 'Type', 'axes');
    for i = 1:numel(A), txt = [txt; txtOf(A(i).XTickLabel); txtOf(A(i).YTickLabel)]; end %#ok<AGROW>
    L = findall(fig, 'Type', 'legend'); for i = 1:numel(L), txt = [txt; txtOf(L(i).String)]; end %#ok<AGROW>
    all_ = strjoin(txt, ' | ');
    hit = ids(arrayfun(@(s) contains(all_, s), ids));
    assert(isempty(hit), 'participant identifier found in the figure text (%d hits); not exported', numel(hit));
    assert(~contains(all_, 'sub-'), 'a "sub-" string is present in the figure text; not exported');
    fprintf('[check] identifiers: none of %d participant identifiers in %d text strings\n', numel(unique(regexprep(ids,'^sub-',''))), numel(txt));
end
function s = txtOf(v)
    % text content as a column of strings, one per line (a char row vector is
    % one string, not one string per character)
    if isempty(v), s = strings(0,1);
    elseif ischar(v), s = string(cellstr(v));
    elseif iscell(v), s = string(v(:));
    else, s = string(v(:));
    end
end
function titleFit(fig)
    % every (left-aligned) title must end inside its plot box
    A = findall(fig, 'Type', 'axes'); worst = -Inf;
    for i = 1:numel(A)
        t = A(i).Title; if isempty(t.String), continue; end
        u = t.Units; t.Units = 'normalized'; e = t.Extent; t.Units = u;
        ua = A(i).Units; A(i).Units = 'centimeters'; p = A(i).Position; A(i).Units = ua;
        over = (e(1) + e(3) - 1) * p(3) * 10;                 % mm beyond the right end of the plot box
        worst = max(worst, over);
        if over > 0, fprintf('[check] TITLE OVERHANG %.2f mm: "%s"\n', over, strjoin(txtOf(t.String)', ' ')); end
    end
    fprintf('[check] titles: largest overhang beyond the plot box %.2f mm (<= 0 means inside)\n', worst);
    assert(worst <= 0, 'a panel title runs beyond its plot box');
end

%% myComparePlotIntensity_less50.m
% Figure-S6-style capillary-density comparison for the 10, 20, 30, 40 and 50%
% residual-intensity compartments.
%
% Required inputs (made by findROI_intensity_less50.m):
%   IntensityLeft10_SIvalid5_ROI_WTAD_full.mat ...
%   IntensityLeft50_SIvalid5_ROI_WTAD_full.mat
%
% Outputs are written beside this script:
%   FigureS6Style_CapillaryDensity_Intensity_Left10to50.png/.pdf/.eps
%   FigureS6Style_CapillaryDensity_Intensity_Slopes.png/.pdf/.eps
%   Scatter_WT_AD_Intensity_Left<percent>.png/.pdf
%   FigureS6Style_CapillaryDensity_Intensity_Left10to50_stats.csv/.mat

clc; clear; close all;

%% -------------------------------- settings --------------------------------
thisDir = fileparts(mfilename('fullpath'));
vessDir = fullfile(thisDir,'data');   % inputs in data/
addpath(thisDir);  % layoutScatterKey.m

xlsxFile = fullfile(vessDir,'1-s2.0-S0896627321000805-mmc2.xlsx');
leftPercents = 10:10:50;
capillaryColumn = 'Var17';
analysisSuffix = '_SIvalid5';
minimumAnimalsPerROI = 6;  % at least half of each 11-mouse group

commonXLimits = [0 2.0];
% Same axis as myComparePlotIntensity_marked.m, and the same scale in every
% panel so threshold changes are not exaggerated by rescaling.
fixedYLimits = [-1 2];

cWT = [0.0000 0.4470 0.7410];
cAD = [0.8500 0.3250 0.0980];

% Match myComparePlotIntensity_marked.m, except that all ROIs remain plain
% background dots and no highlighted-marker key is added.
BG_SIZE  = 30;
BG_ALPHA = 0.35;
FIG_W    = 400;
FIG_H    = 400;
FS_AXIS  = 16;
FS_FIT   = 12;

CP = readtable(xlsxFile);
if ~ismember(capillaryColumn,CP.Properties.VariableNames)
    error('IntensityLess50:MissingCapillaryColumn', ...
        'Column %s is not present in %s.',capillaryColumn,xlsxFile);
end

%% -------------------------- load data and fit -------------------------------
nThreshold = numel(leftPercents);
Panel = repmat(struct(),nThreshold,1);

wtSlope = nan(nThreshold,1); wtIntercept = nan(nThreshold,1);
wtR = nan(nThreshold,1); wtP = nan(nThreshold,1); wtNROI = nan(nThreshold,1);
adSlope = nan(nThreshold,1); adIntercept = nan(nThreshold,1);
adR = nan(nThreshold,1); adP = nan(nThreshold,1); adNROI = nan(nThreshold,1);
wtMouseMeanSlope = nan(nThreshold,1); wtMouseSEM = nan(nThreshold,1); wtN = nan(nThreshold,1);
adMouseMeanSlope = nan(nThreshold,1); adMouseSEM = nan(nThreshold,1); adN = nan(nThreshold,1);
WTMouseSlopes = cell(nThreshold,1);
ADMouseSlopes = cell(nThreshold,1);
allDisplayedY = [];

for t = 1:nThreshold
    leftPercent = leftPercents(t);
    newROIFile = fullfile(vessDir,sprintf( ...
        'IntensityLeft%d%s_ROI_WTAD_full.mat',leftPercent,analysisSuffix));
    if ~exist(newROIFile,'file')
        error('IntensityLess50:MissingROITable', ...
            ['Missing %s. Run findROI_intensity_less50.m after the WT ' ...
             'segmentation and random-group analysis finish.'],newROIFile);
    end

    % All five thresholds come from the same NaN-aware tables
    % (IntensityLeft<percent>_SIvalid5).
    MouseData = load(newROIFile,'WTAnimalROI','ADAnimalROI','WTNames','ADNames');
    pooledROIFile = newROIFile;
    PooledData = load(pooledROIFile,'ROI');
    ROI = PooledData.ROI;
    ROI.SourceRow = (1:height(ROI))';
    merged = innerjoin(CP,ROI,'LeftKeys','Acronym','RightKeys','Name');

    x = str2double(string(merged.(capillaryColumn)));
    withinX = isfinite(x) & x >= commonXLimits(1) & x <= commonXLimits(2);
    if ismember('WTn',merged.Properties.VariableNames)
        sufficientWT = merged.WTn >= minimumAnimalsPerROI;
        sufficientAD = merged.ADn >= minimumAnimalsPerROI;
    else
        sufficientWT = true(height(merged),1);
        sufficientAD = true(height(merged),1);
    end
    validWT = withinX & isfinite(merged.WTmean) & sufficientWT;
    validAD = withinX & isfinite(merged.ADmean) & sufficientAD;

    fitWT = fitRegionalLine(x,merged.WTmean,validWT,commonXLimits);
    fitAD = fitRegionalLine(x,merged.ADmean,validAD,commonXLimits);
    wtSlope(t) = fitWT.slope; wtIntercept(t) = fitWT.intercept;
    wtR(t) = fitWT.r; wtP(t) = fitWT.p; wtNROI(t) = fitWT.n;
    adSlope(t) = fitAD.slope; adIntercept(t) = fitAD.intercept;
    adR(t) = fitAD.r; adP(t) = fitAD.p; adNROI(t) = fitAD.n;

    sourceRows = double(merged.SourceRow);
    WTMouseSlopes{t} = fitMouseLines(x,validWT,MouseData.WTAnimalROI(:,sourceRows));
    ADMouseSlopes{t} = fitMouseLines(x,validAD,MouseData.ADAnimalROI(:,sourceRows));
    wtMouseMeanSlope(t) = mean(WTMouseSlopes{t},'omitnan');
    adMouseMeanSlope(t) = mean(ADMouseSlopes{t},'omitnan');
    wtN(t) = sum(isfinite(WTMouseSlopes{t}));
    adN(t) = sum(isfinite(ADMouseSlopes{t}));
    wtMouseSEM(t) = std(WTMouseSlopes{t},0,'omitnan') / sqrt(wtN(t));
    adMouseSEM(t) = std(ADMouseSlopes{t},0,'omitnan') / sqrt(adN(t));

    Panel(t).xWT = x(validWT);
    Panel(t).yWT = merged.WTmean(validWT);
    Panel(t).xAD = x(validAD);
    Panel(t).yAD = merged.ADmean(validAD);
    Panel(t).fitWT = fitWT;
    Panel(t).fitAD = fitAD;
    allDisplayedY = [allDisplayedY; Panel(t).yWT; Panel(t).yAD; ...
        fitWT.y(:); fitAD.y(:)]; %#ok<AGROW>

    fprintf(['%2d%% left: WT slope %.5f, r %.3f, p %.3g (ROI n=%d); ' ...
             'AD slope %.5f, r %.3f, p %.3g (ROI n=%d)\n'], ...
        leftPercent,fitWT.slope,fitWT.r,fitWT.p,fitWT.n, ...
        fitAD.slope,fitAD.r,fitAD.p,fitAD.n);
end

if isempty(fixedYLimits)
    commonYLimits = paddedLimits(allDisplayedY,0.06);
else
    commonYLimits = fixedYLimits;
end

%% ------------------------------ five-panel figure ---------------------------
fig = figure('Color','w','Position',[50 100 1800 400]);
tl = tiledlayout(fig,1,5,'TileSpacing','compact','Padding','compact');

for t = 1:nThreshold
    ax = nexttile(tl);
    drawPanel(ax,Panel(t),cWT,cAD, ...
        commonXLimits,commonYLimits,BG_SIZE,BG_ALPHA,FS_AXIS);
    title(ax,sprintf('%d%% residual voxels',leftPercents(t)), ...
        'FontWeight','normal','FontSize',16);
    xlim('tight')
    % addFitText(ax,Panel(t).fitWT,Panel(t).fitAD,cWT,cAD);
end

xlabel(tl,'Capillary Length Density (m/mm^3)', ...
    'Interpreter','none','FontSize',17);
ylabel(tl,'Mean Normalized Signal Intensity (a.u.)','FontSize',17);
title(tl,'Residual-intensity sensitivity (5% minimum valid block coverage)', ...
    'FontSize',20,'FontWeight','normal');

outBase = fullfile(thisDir, ...
    ['FigureS6Style_CapillaryDensity_Intensity_Left10to50' analysisSuffix]);
exportgraphics(fig,[outBase '.png'],'Resolution',300);
exportgraphics(fig,[outBase '.pdf'],'ContentType','vector');
exportgraphics(fig,[outBase '.eps'],'ContentType','vector');

%% ---------------------------- individual panels -----------------------------
for t = 1:nThreshold
    leftPercent = leftPercents(t);
    oneFig = figure('Color','w','Position',[100 100 FIG_W FIG_H]);
    ax = axes(oneFig);
    [hWT,hWTFit,hAD,hADFit] = drawPanel(ax,Panel(t),cWT,cAD, ...
        [],commonYLimits,BG_SIZE,BG_ALPHA,FS_AXIS);
    xlabel(ax,'Capillary Length Density (m/mm^3)','Interpreter','none');
    ylabel(ax,{'Normalized Signal Intensity' '(a.u.)'});

    wtLegend = sprintf('WT Fit: y = %.3fx + %.3f (r=%.3f, p=%.3g)', ...
        Panel(t).fitWT.slope,Panel(t).fitWT.intercept,Panel(t).fitWT.r,Panel(t).fitWT.p);
    adLegend = sprintf('AD Fit: y = %.3fx + %.3f (r=%.3f, p=%.3g)', ...
        Panel(t).fitAD.slope,Panel(t).fitAD.intercept,Panel(t).fitAD.r,Panel(t).fitAD.p);
    lgFit = legend(ax,[hWT hWTFit hAD hADFit], ...
        {'WTmean Data',wtLegend,'ADmean Data',adLegend}, ...
        'Location','southoutside','Box','off');
    fontsize(oneFig,FS_AXIS,'points');
    lgFit.FontSize = FS_FIT;
    layoutScatterKey(ax,lgFit,[],[]);

    oneBase = fullfile(thisDir,sprintf( ...
        'Scatter_WT_AD_Intensity_Left%d%s',leftPercent,analysisSuffix));
    exportgraphics(oneFig,[oneBase '.png'],'Resolution',300);
    exportgraphics(oneFig,[oneBase '.pdf'],'ContentType','vector');
end

%% ---------------------- slope-versus-threshold figure -----------------------
slopeFig = figure('Color','w','Position',[100 100 800 400]);
ax = axes(slopeFig);
hold(ax,'on');
% yline(ax,0,':','Color',[0.45 0.45 0.45],'HandleVisibility','off');
hWT = errorbar(ax,leftPercents,wtSlope,wtMouseSEM,'-o', ...
    'Color',cWT,'MarkerFaceColor',cWT,'MarkerSize',8, ...
    'LineWidth',2,'CapSize',11);
hAD = errorbar(ax,leftPercents,adSlope,adMouseSEM,'-s', ...
    'Color',cAD,'MarkerFaceColor',cAD,'MarkerSize',8, ...
    'LineWidth',2,'CapSize',11);
hold(ax,'off');

set(ax,'XTick',leftPercents,'XLim',[5 55],'FontSize',14, ...
    'LineWidth',1,'Box','off','Layer','top');
grid(ax,'off');
xlabel(ax,'Residual voxel fraction (%)','FontSize',17);
ylabel(ax,'Capillary-density regression slope','FontSize',17);
% title(ax,'Pooled regional slope by residual compartment', ...
%     'FontSize',20,'FontWeight','normal');
legend(ax,[hWT hAD],{'WT','AD'},'Location','best','Box','off','FontSize',15);

% note = sprintf(['Points: pooled ROI fits. Error bars: SEM of mouse-level slopes ' ...
%                 '(WT n=%s; AD n=%s).'],countLabel(wtN),countLabel(adN));
% annotation(slopeFig,'textbox',[0.08 0.015 0.84 0.07], ...
%     'String',note,'EdgeColor','none','HorizontalAlignment','center', ...
%     'VerticalAlignment','middle','FontSize',11,'Color',[0.30 0.33 0.38]);

slopeBase = fullfile(thisDir, ...
    ['FigureS6Style_CapillaryDensity_Intensity_Slopes' analysisSuffix]);
exportgraphics(slopeFig,[slopeBase '.png'],'Resolution',300);
exportgraphics(slopeFig,[slopeBase '.pdf'],'ContentType','vector');
exportgraphics(slopeFig,[slopeBase '.eps'],'ContentType','vector');

%% --------------------------- result table/cache -----------------------------
ResidualPercent = leftPercents';
WT_PooledSlope = wtSlope;
WT_Intercept = wtIntercept;
WT_R = wtR;
WT_P = wtP;
WT_NROI = wtNROI;
WT_MouseMeanSlope = wtMouseMeanSlope;
WT_MouseSlopeSEM = wtMouseSEM;
WT_N = wtN;
AD_PooledSlope = adSlope;
AD_Intercept = adIntercept;
AD_R = adR;
AD_P = adP;
AD_NROI = adNROI;
AD_MouseMeanSlope = adMouseMeanSlope;
AD_MouseSlopeSEM = adMouseSEM;
AD_N = adN;

Results = table(ResidualPercent, ...
    WT_PooledSlope,WT_Intercept,WT_R,WT_P,WT_NROI, ...
    WT_MouseMeanSlope,WT_MouseSlopeSEM,WT_N, ...
    AD_PooledSlope,AD_Intercept,AD_R,AD_P,AD_NROI, ...
    AD_MouseMeanSlope,AD_MouseSlopeSEM,AD_N);

statsBase = fullfile(thisDir, ...
    ['FigureS6Style_CapillaryDensity_Intensity_Left10to50_stats' analysisSuffix]);
writetable(Results,[statsBase '.csv']);
save([statsBase '.mat'],'Results','WTMouseSlopes','ADMouseSlopes', ...
    'commonXLimits','commonYLimits','capillaryColumn','analysisSuffix', ...
    'minimumAnimalsPerROI');

fprintf('\nSaved the combined, individual and slope figures in %s.\n',thisDir);
disp(Results);

%% -------------------------------- functions --------------------------------
function fitResult = fitRegionalLine(x,y,use,~)
    xUse = double(x(use));
    yUse = double(y(use));
    if numel(xUse) < 3 || numel(unique(xUse)) < 2
        error('IntensityLess50:TooFewROIs','At least three usable ROIs are required.');
    end
    coefficient = polyfit(xUse,yUse,1);
    [R,P] = corrcoef(xUse,yUse);
    fitResult.slope = coefficient(1);
    fitResult.intercept = coefficient(2);
    fitResult.r = R(1,2);
    fitResult.p = P(1,2);
    fitResult.n = numel(xUse);
    fitResult.x = linspace(min(xUse),max(xUse),200);
    fitResult.y = polyval(coefficient,fitResult.x);
end

function slopes = fitMouseLines(x,pooledUse,animalY)
    nMouse = size(animalY,1);
    slopes = nan(nMouse,1);
    for m = 1:nMouse
        use = pooledUse & isfinite(animalY(m,:)');
        if sum(use) >= 3 && numel(unique(x(use))) >= 2
            coefficient = polyfit(x(use),animalY(m,use)',1);
            slopes(m) = coefficient(1);
        end
    end
end

function [hWT,hWTFit,hAD,hADFit] = drawPanel(ax,panel,cWT,cAD, ...
        xLimits,yLimits,bgSize,bgAlpha,fontSize)
    hold(ax,'on');
    hWT = scatter(ax,panel.xWT,panel.yWT,bgSize,'filled', ...
        'MarkerFaceColor',cWT,'MarkerEdgeColor','none','MarkerFaceAlpha',bgAlpha);
    hWTFit = plot(ax,panel.fitWT.x,panel.fitWT.y,'-', ...
        'Color',cWT,'LineWidth',2);
    hAD = scatter(ax,panel.xAD,panel.yAD,bgSize,'filled', ...
        'MarkerFaceColor',cAD,'MarkerEdgeColor','none','MarkerFaceAlpha',bgAlpha);
    hADFit = plot(ax,panel.fitAD.x,panel.fitAD.y,'-', ...
        'Color',cAD,'LineWidth',2);
    hold(ax,'off');
    if isempty(xLimits)
        xlim(ax,'tight');
    else
        xlim(ax,xLimits);
    end
    if ~isempty(yLimits)
        ylim(ax,yLimits);
    end
    grid(ax,'off');
    set(ax,'FontSize',fontSize,'LineWidth',1,'Box','off','Layer','top');
end

function addFitText(ax,fitWT,fitAD,cWT,cAD)
    wtText = sprintf(['WT: y = %.3fx %+.3f\n' ...
        'r = %.2f, p = %.2g'], ...
        fitWT.slope,fitWT.intercept,fitWT.r,fitWT.p);
    adText = sprintf(['AD: y = %.3fx %+.3f\n' ...
        'r = %.2f, p = %.2g'], ...
        fitAD.slope,fitAD.intercept,fitAD.r,fitAD.p);

    text(ax,0.03,0.16,wtText, ...
        'Units','normalized','Color',cWT,'FontSize',9, ...
        'VerticalAlignment','bottom');
    text(ax,0.03,0.025,adText, ...
        'Units','normalized','Color',cAD,'FontSize',9, ...
        'VerticalAlignment','bottom');
end

function lim = paddedLimits(values,fraction)
    values = values(isfinite(values));
    if isempty(values)
        lim = [-1 1];
        return;
    end
    low = min(values);
    high = max(values);
    span = high-low;
    if span == 0
        span = max(abs(low),1);
    end
    lim = [low-fraction*span high+fraction*span];
end

function out = countLabel(n)
    uniqueN = unique(n(isfinite(n)));
    if isscalar(uniqueN)
        out = sprintf('%d',uniqueN);
    else
        out = mat2str(n(:)');
    end
end

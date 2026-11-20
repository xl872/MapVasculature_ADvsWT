%% myComparePlotIntensity_marked.m
%  Fig S7: regional residual signal intensity (WT, AD) against the capillary
%  length density of the Ji et al. 2021 atlas, with the ROIs of interest drawn
%  as their own markers. Intensity counterpart of myComparePlotDensity_marked.m.
%
%  x  : the paper's capillary length density (mergedTable.Var17).
%  y  : for the marked ROIs it comes from the hand-drawn intensity ROI masks
%       (Intensity1_ROIsel_WTAD.mat, made by myROIforCompare_selected_intensity.m),
%       so a marker is the signal intensity of the patch findROI_intensity.m
%       drew, not the whole-structure atlas average.
%  The fits, r and p are computed over every valid atlas ROI of the
%  hierarchy-filled Intensity1_ROI_WTAD_full.
%
%  Run first:  myAggregateHigherROIs.m            (writes Intensity1_ROI_WTAD_full)
%              myROIforCompare_selected_intensity.m (writes Intensity1_ROIsel_WTAD)

clc;clear all;close all;
DATA_DIR = fullfile(fileparts(mfilename('fullpath')),'data');
CP=readtable(fullfile(DATA_DIR,'1-s2.0-S0896627321000805-mmc2.xlsx'));

%% ------------------------------- SETTINGS --------------------------------
tag = 1;                 % the intensity data only exists for tag 1

% x value of the marked ROIs.
% 'weighted' -> the spreadsheet capillary density averaged over the atlas labels
%               the mask actually covers, weighted by voxel count
%               (ROIsel_xvalue.mat, from myROIxLengthDensity_selected.m).
% 'acronym'  -> the single spreadsheet row for hiAcr.
X_SOURCE = 'weighted';
load(fullfile(DATA_DIR,'ROIsel_xvalue'))      % -> XSEL

SHOW_LABELS = false;      % only two markers here, so the names fit on the plot
BG_SIZE  = 30;           % plain dots: size and opacity
BG_ALPHA = 0.35;
HI_SIZE  = 9;            % highlighted markers

FIG_W = 500; FIG_H = 400;   % figure size in pixels
FS_AXIS = 16;               % axis labels and ticks
FS_FIT  = 12;               % legend with the fit equations
FS_KEY  = 12;               % ROI marker key (shrinks itself if it will not fit)

cWT = '#0072BD';
cAD = '#D95319';

% Shapes: 7 easily separated forms, filled first then hollow.
SHAPES = {'o','s','d','^','v','>','<'};

% load(fullfile(DATA_DIR,['Intensity',num2str(tag),'_ROI_WTAD_full']))
load(fullfile(DATA_DIR,'IntensityLeft50_SIvalid5_ROI_WTAD_full.mat'))

load(fullfile(DATA_DIR,['Intensity',num2str(tag),'_ROIsel_WTAD']))
mergedTable = innerjoin(CP, ROI, 'LeftKeys', 'Acronym', 'RightKeys', 'Name');

% the ROIs that get their own marker, and where their y value comes from
hiName = ROIsel.Name;
hiAcr  = ROIsel.Acronym;
hiWT   = ROIsel.WTmean;
hiAD   = ROIsel.ADmean;
nHi    = numel(hiAcr);

hiMk   = cell(nHi,1);  hiFilled = false(nHi,1);
for i = 1:nHi
    hiMk{i}     = SHAPES{mod(i-1,numel(SHAPES))+1};
    hiFilled(i) = i <= numel(SHAPES);
end

%%
% Var17 = "Capillary length density (m/mm^3), Mean" in the spreadsheet
target_col_name = 'Var17';
x_raw = mergedTable.(target_col_name);
x_num = str2double(x_raw);
y_wt_raw = mergedTable.WTmean;
y_ad_raw = mergedTable.ADmean;

% WTmean
valid_idx_wt = ~isnan(x_num) & (y_wt_raw ~= 0);
x_wt = x_num(valid_idx_wt);
y_wt = y_wt_raw(valid_idx_wt);
a_wt = mergedTable.Acronym(valid_idx_wt);      % acronym of every plotted WT point

% ADmean
valid_idx_ad = ~isnan(x_num) & (y_ad_raw ~= 0);
x_ad = x_num(valid_idx_ad);
y_ad = y_ad_raw(valid_idx_ad);
a_ad = mergedTable.Acronym(valid_idx_ad);      % acronym of every plotted AD point

% WTmean
[R_wt, P_wt] = corrcoef(x_wt, y_wt);
r_wt = R_wt(1, 2);
p_wt = P_wt(1, 2);
coeff_wt = polyfit(x_wt, y_wt, 1);
x_fit_wt = linspace(min(x_wt), max(x_wt), 100);
y_fit_wt = polyval(coeff_wt, x_fit_wt);

% ADmean
[R_ad, P_ad] = corrcoef(x_ad, y_ad);
r_ad = R_ad(1, 2);
p_ad = P_ad(1, 2);
coeff_ad = polyfit(x_ad, y_ad, 1);
x_fit_ad = linspace(min(x_ad), max(x_ad), 100);
y_fit_ad = polyval(coeff_ad, x_fit_ad);

%% ------------------------------- figure ---------------------------------
fig = figure('Position',[100 100 FIG_W FIG_H]);
ax  = axes(fig);
hold(ax,'on');

isHiWt = ismember(a_wt, hiAcr);
isHiAd = ismember(a_ad, hiAcr);

% everything that is not highlighted stays a plain dot
h1 = scatter(ax, x_wt(~isHiWt), y_wt(~isHiWt), BG_SIZE, 'filled', ...
             'MarkerFaceColor', cWT, 'MarkerEdgeColor', 'none', 'MarkerFaceAlpha', BG_ALPHA);
h2 = plot(ax, x_fit_wt, y_fit_wt, '-', 'Color', cWT, 'LineWidth', 2);

h3 = scatter(ax, x_ad(~isHiAd), y_ad(~isHiAd), BG_SIZE, 'filled', ...
             'MarkerFaceColor', cAD, 'MarkerEdgeColor', 'none', 'MarkerFaceAlpha', BG_ALPHA);
h4 = plot(ax, x_fit_ad, y_fit_ad, '-', 'Color', cAD, 'LineWidth', 2);

% the ROIs of interest: x from the paper, y from the hand-drawn intensity ROI
missing = {};
for i = 1:nHi
    xi = NaN;
    switch X_SOURCE
        case 'weighted'
            q = find(strcmp(XSEL.Name, hiName{i}), 1);
            if ~isempty(q), xi = XSEL.CLD_est(q); end
        case 'acronym'
            j = find(strcmp(mergedTable.Acronym, hiAcr{i}), 1);
            if ~isempty(j), xi = x_num(j); end
        otherwise
            error('X_SOURCE must be ''weighted'' or ''acronym''');
    end
    if isnan(xi)
        missing{end+1} = hiName{i}; continue; %#ok<SAGROW>
    end

    if hiFilled(i)
        argWT = {'MarkerFaceColor', cWT, 'MarkerEdgeColor', 'k'};
        argAD = {'MarkerFaceColor', cAD, 'MarkerEdgeColor', 'k'};
    else
        argWT = {'MarkerFaceColor', 'none', 'MarkerEdgeColor', cWT, 'LineWidth', 1.8};
        argAD = {'MarkerFaceColor', 'none', 'MarkerEdgeColor', cAD, 'LineWidth', 1.8};
    end

    plot(ax, xi, hiWT(i), 'LineStyle','none', 'Marker', hiMk{i}, ...
         'MarkerSize', HI_SIZE, argWT{:}, 'HandleVisibility','off');
    plot(ax, xi, hiAD(i), 'LineStyle','none', 'Marker', hiMk{i}, ...
         'MarkerSize', HI_SIZE, argAD{:}, 'HandleVisibility','off');

    if SHOW_LABELS
        text(ax, xi, hiWT(i), ['  ' hiName{i}], 'FontSize', FS_KEY-2, ...
             'Color', [0.25 0.25 0.25], 'VerticalAlignment','middle', ...
             'HorizontalAlignment','left', 'Clipping','on');
    end
end
if ~isempty(missing)
    fprintf('no x value for %s\n', strjoin(missing,', '));
end

hold(ax,'off');
xlim('tight');
ylim([-1 1.5]);
xlabel(ax, 'Capillary Length Density (m/mm^3)', 'Interpreter', 'none');
ylabel(ax, {'Normalized Signal Intensity' '(a.u.)'});
% title(ax, 'Linear Regression: WTmean and ADmean');
grid(ax,'off');

str_wt_data = 'WTmean Data';
str_wt_fit  = sprintf('WT Fit: y = %.3fx + %.3f (r=%.3f, p=%.3g)', coeff_wt(1), coeff_wt(2), r_wt, p_wt);
str_ad_data = 'ADmean Data';
str_ad_fit  = sprintf('AD Fit: y = %.3fx + %.3f (r=%.3f, p=%.3g)', coeff_ad(1), coeff_ad(2), r_ad, p_ad);

lgFit = legend(ax, [h1 h2 h3 h4], {str_wt_data, str_wt_fit, str_ad_data, str_ad_fit}, ...
       'Location', 'southoutside', 'Box', 'off');

% font sizes first - the layout below measures the legends, so they have to be
% at their final size before it runs
fontsize(fig, FS_AXIS, "points")
lgFit.FontSize = FS_FIT;

axKey = axes(fig, 'Position', ax.Position, 'Color','none', ...
             'XColor','none', 'YColor','none', 'XLim',[0 1], 'YLim',[0 1]);
hold(axKey,'on');
hKey = gobjects(nHi,1);
for i = 1:nHi
    if hiFilled(i)
        hKey(i) = plot(axKey, NaN, NaN, 'LineStyle','none', 'Marker', hiMk{i}, ...
                       'MarkerSize', 9, 'MarkerFaceColor','k', 'MarkerEdgeColor','k');
    else
        hKey(i) = plot(axKey, NaN, NaN, 'LineStyle','none', 'Marker', hiMk{i}, ...
                       'MarkerSize', 9, 'MarkerFaceColor','none', 'MarkerEdgeColor','k', ...
                       'LineWidth', 1.5);
    end
end
hold(axKey,'off');
lgKey = legend(axKey, hKey, hiName, 'Box','off', 'FontSize', FS_KEY);
title(lgKey, 'ROI');

% works out the axes box and both legend positions from their measured sizes,
% so the figure can be resized and the font changed without hand-tuning
layoutScatterKey(ax, lgFit, axKey, lgKey);

%%
exportgraphics(fig, ['Scatter_WT_AD_Intensity_marked',num2str(tag),'.png'], 'Resolution', 300);
exportgraphics(fig, ['Scatter_WT_AD_Intensity_marked',num2str(tag),'.eps'], 'ContentType', 'vector');

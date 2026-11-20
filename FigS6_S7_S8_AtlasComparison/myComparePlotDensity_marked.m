%% myComparePlotDensity_marked.m
%  Fig S6 A-E: regional vessel density (WT, AD) against the vessel length
%  density of the Ji et al. 2021 atlas. With SHOW_MARKED = true (Fig S6 A) the
%  ROIs of interest are drawn as their own markers with a key; with false
%  (Fig S6 B-E) every ROI is a plain dot.
%
%  x  : the paper's value (mergedTable.ThisStudy) for the acronym.
%  y  : for the marked ROIs it comes from the hand-drawn bar-plot ROI masks
%       (Density<tag>_ROIsel_WTAD.mat, made by myROIforCompare_selected.m)
%       measured on this tag's own data - so a marker is the vessel density of
%       the region the box plots use, not the whole-structure atlas average, and
%       it moves with tag just like every other point on the figure.
%  The fits, r and p are computed over every valid atlas ROI from
%  Density<tag>_ROI_WTAD_full.
%
%  Run myROIforCompare_selected.m once before this.

clc;clear all;close all;
DATA_DIR = fullfile(fileparts(mfilename('fullpath')),'data');
CP=readtable(fullfile(DATA_DIR,'1-s2.0-S0896627321000805-mmc2.xlsx'));

%% ------------------------------- SETTINGS --------------------------------
% 'tag'  -> the bar-plot ROI masks measured on this tag's own Density<tag>WTAD
%           data, so the markers follow the tag the same way the table ROIs do.
%           Mean/std over the 100 repeats, same definition as the table columns.
%           For tag=1 this tracks the ROIbar.m values at r=0.9999 (WT sits a
%           uniform ~0.037 lower, AD matches to +/-0.01).
% 'bar'    -> the exact ROIbar.m numbers (22 animals, forROIbarplot.mat). That
%           data is tag=1 only, so the markers would not change with tag.
MARKER_SOURCE = 'tag';

% x value of the marked ROIs.
% 'weighted' -> the spreadsheet density averaged over the atlas labels the mask
%               actually covers, weighted by voxel count (ROIsel_xvalue.mat,
%               written by myROIxLengthDensity_selected.m). Matches how the y
%               value is built, and fixes masks that are not the structure their
%               name suggests (mEC is mostly ENTl, IC misses ICc, ACC is ~half MOs).
% 'acronym'  -> the single spreadsheet row for hiAcr.
X_SOURCE = 'weighted';
load(fullfile(DATA_DIR,'ROIsel_xvalue'))      % -> XSEL

% true  -> the selected ROIs get their own marker plus the key on the right, and
%          their plain dot is removed so each one appears once.
% false -> no markers and no key: every ROI is drawn as a plain dot, including
%          the selected ones, and the axes take the full width.
% Either way the fits, r and p are identical - they always use every valid atlas
% ROI. Output goes to ..._marked<tag> or ..._plain<tag> so the two do not
% overwrite each other.
SHOW_MARKED = true;

SHOW_LABELS = false;     % also print the ROI name next to its point (they collide
                         % in the crowded middle, the key on the right is enough)
BG_SIZE  = 30;           % plain dots: size and opacity
BG_ALPHA = 0.35;
HI_SIZE  = 8;            % highlighted markers

FIG_W = 500; FIG_H = 400;   % figure size in pixels 500 for marker 400 for plain
FS_AXIS = 16;               % axis labels and ticks
FS_FIT  = 12;               % legend with the fit equations
FS_KEY  = 12;               % ROI marker key (shrinks itself if it will not fit)

cWT = '#0072BD';
cAD = '#D95319';

% Shapes: 6 easily separated forms (no pentagram/hexagram: they look alike at
% this size), filled first, then hollow.
SHAPES = {'o','s','d','^','v','>'};

for tag=[1 6 11 16 21 -1 -6 -11 -16 -21]
load(fullfile(DATA_DIR,['Density',num2str(tag),'_ROI_WTAD_full']))
load(fullfile(DATA_DIR,['Density',num2str(tag),'_ROIsel_WTAD']))
mergedTable = innerjoin(CP, ROI, 'LeftKeys', 'Acronym', 'RightKeys', 'Name');

% the ROIs that get their own marker, and where their y value comes from
hiName = ROIsel.Name;
hiAcr  = ROIsel.Acronym;
nHi    = numel(hiAcr);
switch MARKER_SOURCE
    case 'bar', hiWT = ROIsel.WTmean_bar; hiAD = ROIsel.ADmean_bar;
    case 'tag', hiWT = ROIsel.WTmean;     hiAD = ROIsel.ADmean;
    otherwise,  error('MARKER_SOURCE must be ''bar'' or ''tag''');
end
if strcmp(MARKER_SOURCE,'bar') && tag ~= 1
    fprintf(['tag %d: markers show the tag=1 bar-plot values (forROIbarplot.mat), ' ...
             'while the cloud and the fit are tag %d.\n'], tag, tag);
end

hiMk   = cell(nHi,1);  hiFilled = false(nHi,1);
for i = 1:nHi
    hiMk{i}     = SHAPES{mod(i-1,numel(SHAPES))+1};
    hiFilled(i) = i <= numel(SHAPES);
end

%%
target_col_name = 'ThisStudy';
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

% with the markers off, nothing is held back from the plain scatter
if SHOW_MARKED
    isHiWt = ismember(a_wt, hiAcr);
    isHiAd = ismember(a_ad, hiAcr);
else
    isHiWt = false(size(a_wt));
    isHiAd = false(size(a_ad));
end

% everything that is not highlighted stays a plain dot
h1 = scatter(ax, x_wt(~isHiWt), y_wt(~isHiWt), BG_SIZE, 'filled', ...
             'MarkerFaceColor', cWT, 'MarkerEdgeColor', 'none', 'MarkerFaceAlpha', BG_ALPHA);
h2 = plot(ax, x_fit_wt, y_fit_wt, '-', 'Color', cWT, 'LineWidth', 2);

h3 = scatter(ax, x_ad(~isHiAd), y_ad(~isHiAd), BG_SIZE, 'filled', ...
             'MarkerFaceColor', cAD, 'MarkerEdgeColor', 'none', 'MarkerFaceAlpha', BG_ALPHA);
h4 = plot(ax, x_fit_ad, y_fit_ad, '-', 'Color', cAD, 'LineWidth', 2);

% the ROIs of interest: x from the paper, y from the bar-plot ROI
missing = {};
for i = 1:nHi*SHOW_MARKED
    xi = NaN;
    switch X_SOURCE
        case 'weighted'
            q = find(strcmp(XSEL.Name, hiName{i}), 1);
            if ~isempty(q), xi = XSEL.VLD_est(q); end
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
        text(ax, xi, hiWT(i), ['  ' hiName{i}], 'FontSize', 9, ...
             'Color', [0.25 0.25 0.25], 'VerticalAlignment','middle', ...
             'HorizontalAlignment','left', 'Clipping','on');
    end
end
if ~isempty(missing)
    fprintf('tag %d: no x value for %s\n', tag, strjoin(missing,', '));
end

hold(ax,'off');
xlim('tight')
xlabel(ax, 'Vessel Length Density (m/mm^3)', 'Interpreter', 'none');
ylabel(ax, 'Vessel Density');
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

if SHOW_MARKED
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
else
    axKey = [];
    lgKey = [];
end

% works out the axes box and both legend positions from their measured sizes,
% so the figure can be resized and the font changed without hand-tuning
layoutScatterKey(ax, lgFit, axKey, lgKey);
ylim([0 0.6])
%%
if SHOW_MARKED, suffix = 'marked'; else, suffix = 'plain'; end
exportgraphics(fig, ['Scatter_WT_AD_Density_',suffix,num2str(tag),'.png'], 'Resolution', 300);
exportgraphics(fig, ['Scatter_WT_AD_Density_',suffix,num2str(tag),'.eps'], 'ContentType', 'vector');
end

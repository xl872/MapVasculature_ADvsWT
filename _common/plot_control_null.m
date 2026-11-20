function [f, P] = plot_control_null(nullFile, titleStr, varargin)
%PLOT_CONTROL_NULL  One panel of Fig S4 / S5 from data/S4_null_<ROI>.mat or
%data/S5_null_<ROI>.mat.
%
%   plot_control_null('data/S4_null_Cg.mat','ACC')
%   plot_control_null('data/S5_null_HP.mat','HP')
%   [f, P] = plot_control_null(file, title, 'Position', pos)   % figure size (pixels)
%   plot_control_null(file, title, 'Direction', 'WT-AD')       % as stored
%
% Same drawing as findROI.m ("distribution_stat_nonzero", Fig S4) and
% findROI_intensity.m ("distribution", Fig S5): normalised histogram of the
% control-permutation values in the ROI (black), gauss1 fit (red curve),
% the one-sided 95% bound of the fit (red line, b1 +/- 1.96*c1/sqrt(2))
% and the observed group difference (star). The stored N/edges/gauss/xp/xlim
% are those of findROI.m / findROI_intensity.m; nothing is refitted here.
%
% Sign: the .mat files store WT - AD (as computed by findROI*.m). The panel
% is drawn as AD - WT by default, like Fig 4-6: the x-axis is mirrored at
% plot time (edges, N, observed, xp, the Gaussian centre and xlim), the
% .mat files are not changed. 'Direction','WT-AD' draws the stored values.
%
% P returns the plotted values (observed, xp, gauss, xlim, edges, N).
o = struct('Direction','AD-WT', 'Position',[100 100 300 250], 'XLabel','');
for k = 1:2:numel(varargin), o.(varargin{k}) = varargin{k+1}; end
S = load(nullFile);
if nargin < 2 || isempty(titleStr), titleStr = S.name; end
if strcmpi(o.Direction,'AD-WT')            % stored WT - AD -> AD - WT
    S.edges = -fliplr(S.edges); S.N = fliplr(S.N);
    S.observed = -S.observed; S.xp = -S.xp; S.gauss(2) = -S.gauss(2);
    S.xlim = -fliplr(S.xlim);
    d = '(AD - WT)';
else
    d = '(WT - AD)';
end
if isempty(o.XLabel)                       % S4: vessel density; S5: residual signal
    [~, fn] = fileparts(nullFile);
    if startsWith(fn,'S5_'), q = 'Residual signal intensity difference';
    else,                    q = 'Vessel density differential'; end
    o.XLabel = {q, d};
end
ctr = (S.edges(1:end-1)+S.edges(2:end))/2;
f = figure(Position=o.Position);
plot(ctr,S.N,LineWidth=2,Color='k'); hold on
yl = ylim;
scatter(S.observed, yl(2)/25, 250, 'filled', 'pentagram')
xx = linspace(S.xlim(1),S.xlim(2),1000);
plot(xx, S.gauss(1)*exp(-((xx-S.gauss(2))/S.gauss(3)).^2), 'r');
yl = ylim;
line([S.xp S.xp], yl, 'Color','r','LineWidth',1); ylim(yl)
xlim(S.xlim)
xlabel(o.XLabel); ylabel('Distribution of control')
title(titleStr,'FontWeight','normal'); fontsize(13,'points')   % normal weight
P = struct('observed',S.observed, 'xp',S.xp, 'gauss',S.gauss, 'xlim',S.xlim, ...
           'edges',S.edges, 'N',S.N, 'direction',o.Direction);
end

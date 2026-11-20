function plot_ROI_boxplots(valuesFile, roiNames, outPrefix)
%PLOT_ROI_BOXPLOTS  Box plots of Fig 4B / Fig 6B-C from data/*_ROI_values.mat
%(the MdenWT / MdenAD matrices of ROIbar.m / ROIbarIntensity.m).
%
%   plot_ROI_boxplots('data/Fig4B_ROI_values.mat', {'Cg','Septal','PIR','SS','CM','DG', ...
%        'Amygdalar','LEn','MGN','MEn','IC','PAG5'})
%   plot_ROI_boxplots('data/Fig6BC_ROI_values.mat', {'VISp','HP'})
%
% Same drawing/statistics code as ROIbar.m: one box per group, two-sample
% t-test over the rows (22 rows per group = 11 animals x 2 hemisphere copies),
% stars at p<0.05/0.01/0.001.
S = load(valuesFile);
if nargin < 2 || isempty(roiNames), roiNames = cellstr(S.names); end
if nargin < 3, outPrefix = ''; end
for r = 1:numel(roiNames)
    i = find(S.names == string(roiNames{r}));
    if isempty(i), warning('ROI %s not in %s', roiNames{r}, valuesFile); continue; end
    figure(Position=[100 100 150 250]);
    h = boxplot([S.MdenWT(:,i),S.MdenAD(:,i)],'Labels',{'WT','AD'}, ...
        'colors',[0, 0.4470, 0.7410;0.8500, 0.3250, 0.0980]);
    set(h,'linewidth',1)
    xt = get(gca,'XTick'); yl = ylim; yl = yl(1)+0.95*(yl(2)-yl(1));
    [~,p] = ttest2(S.MdenWT(:,i),S.MdenAD(:,i));
    hold on
    if p<0.001,     text(mean(xt([1 2])),yl,'***',HorizontalAlignment='center')
    elseif p<0.01,  text(mean(xt([1 2])),yl,'**', HorizontalAlignment='center')
    elseif p<0.05,  text(mean(xt([1 2])),yl,'*',  HorizontalAlignment='center')
    end
    ylabel(S.ylab); title(roiNames{r}); fontsize(10,"points")
    fprintf('%-12s n=%d/%d  p=%.3g\n', roiNames{r}, size(S.MdenWT,1), size(S.MdenAD,1), p);
    if ~isempty(outPrefix), print(gcf,[outPrefix roiNames{r} '.eps'],'-depsc','-r300'); end
end
end

function plot_FigS5(outdir)
%PLOT_FIGS5  Fig S5 panels: control-permutation null distribution of the
%residual signal intensity in HP and VISp, drawn with
%../_common/plot_control_null.m (findROI_intensity.m drawing).
%Shown as AD - WT, like Fig 6: the S5_null_*.mat files store WT - AD and the
%x-axis is mirrored at plot time. The quantity is the z-scored residual
%signal intensity. Only the figure size is set (see figpos): 2 panels of
%2.424 in, the width at which Fig S5 is placed in the manuscript (4.85 in).
%Writes FigS5_values_AD-WT.csv with the plotted observed value and 95% line.
here = fileparts(mfilename('fullpath'));
addpath(fullfile(here,'..','_common'));
dd = fullfile(here,'..','Fig6_FigS5_CapillaryIntensityMaps','data');
R = {'HP','HP'; 'VISp','VISp'};                 % panel order
n = size(R,1); obs = zeros(n,1); xp = obs; side = strings(n,1);
for k = 1:n
    [f, P] = plot_control_null(fullfile(dd,['S5_null_' R{k,1} '.mat']), R{k,2}, ...
                               'Position', figpos(2.424,2.14,13));
    exportpanel(f, outdir, sprintf('FigS5_%d_%s', k, R{k,2})); close(f);
    obs(k) = P.observed; xp(k) = P.xp;
    if sign(P.observed) == sign(P.xp), side(k) = "same"; else, side(k) = "OPPOSITE"; end
    fprintf('  FigS5 %-7s observed %+.5f  95%% line %+.5f  (%s side)\n', R{k,2}, obs(k), xp(k), side(k));
end
T = table(string(R(:,1)), string(R(:,2)), obs, xp, side, ...
    'VariableNames', {'roi_file','label','observed_ADminusWT','line95_ADminusWT','star_line_side'});
writetable(T, fullfile(outdir,'FigS5_values_AD-WT.csv'));
end

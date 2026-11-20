function plot_FigS4(outdir)
%PLOT_FIGS4  Fig S4 panels: control-permutation null distribution of each of
%the 12 Fig 4B ROIs, drawn with ../_common/plot_control_null.m (findROI.m:
%black histogram curve, red gauss1 fit, red 95% line, star = observed).
%Drawn as AD - WT, like Fig 4: the S4_null_*.mat files store WT - AD and the
%x-axis is mirrored at plot time. Only the figure size is set (see figpos):
%4 columns of 1.625 in in the 6.5-in-wide figure, in the order of R below.
%Writes FigS4_values_AD-WT.csv with the plotted observed value and 95% line.
here = fileparts(mfilename('fullpath'));
addpath(fullfile(here,'..','_common'));
dd = fullfile(here,'..','Fig4_Fig5_FigS4_VesselDensityMaps','data');
R = {'Cg','ACC'; 'Septal','MS'; 'PIR','Pir'; 'SS','S2';          % row 1
     'CM','CM'; 'DG','DG'; 'Amygdalar','Amy-HP'; 'MGN','MGN';     % row 2
     'LEn','lEC'; 'MEn','mEC'; 'IC','IC'; 'PAG5','PAG'};          % row 3
n = size(R,1); obs = zeros(n,1); xp = obs; side = strings(n,1);
for k = 1:n
    [f, P] = plot_control_null(fullfile(dd,['S4_null_' R{k,1} '.mat']), R{k,2}, ...
                               'Position', figpos(1.625,1.47,13));
    exportpanel(f, outdir, sprintf('FigS4_%02d_%s', k, R{k,2})); close(f);
    obs(k) = P.observed; xp(k) = P.xp;
    if sign(P.observed) == sign(P.xp), side(k) = "same"; else, side(k) = "OPPOSITE"; end
    fprintf('  FigS4 %-7s observed %+.5f  95%% line %+.5f  (%s side)\n', R{k,2}, obs(k), xp(k), side(k));
end
T = table(string(R(:,1)), string(R(:,2)), obs, xp, side, ...
    'VariableNames', {'roi_file','label','observed_ADminusWT','line95_ADminusWT','star_line_side'});
writetable(T, fullfile(outdir,'FigS4_values_AD-WT.csv'));
end

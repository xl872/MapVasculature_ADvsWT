function plot_Fig6(outdir)
%PLOT_FIG6  Fig 6A: residual-signal (capillary) differential map in z-scored
%intensity, drawn as AD - WT like Fig 4/5 (data/FD_intensity_tag1.mat stores
%WT - AD, so the sign is flipped): red/yellow = higher residual signal in AD,
%i.e. lower capillary density in AD; blue = lower residual signal in AD, i.e.
%higher capillary density in AD. 18 coronal slices as plot_FD_slices.m, with a
%horizontal colour bar (-0.5 left, 0.5 right, as Fig 5). Fig 6B/C: VISp and HP
%box plots with the drawing code of ROIbarIntensity.m plus the ROI name as
%title; only the figure size is set (see figpos).
S = figstyle();
here = fileparts(mfilename('fullpath'));
dd = fullfile(here,'..','Fig6_FigS5_CapillaryIntensityMaps','data');
D = load(fullfile(dd,'FD_intensity_tag1.mat'),'slices','info');
sl = D.slices; crop = mapcrop(sl);
gap = 0.015; tileW = (4.50-5*gap)/6;
[f, ~] = mappanel(sl, D.info.climB, -1, D.info.climF, sl.yIdx, 6, tileW, gap, S, crop);
mapH = f.Position(4);                      % mappanel figures are in inches
exportpanel(f,outdir,'Fig6A_signal_map_AD-WT','raster'); close(f);
f = cbpanel(2.2,0.45,D.info.climF,[-0.5 0.5],'Vessel Signal Intensity (AD-WT), p<0.01','horizontal','normal');
exportpanel(f,outdir,'Fig6A_colorbar'); close(f);

V = load(fullfile(dd,'Fig6BC_ROI_values.mat'));
names = {'VISp','HP'};
for r = 1:2
    i = find(string(V.names) == names{r});
    MdenWT = V.MdenWT; MdenAD = V.MdenAD;
    f = figure(Position=figpos(0.98,0.45+mapH,10));
    h = boxplot([MdenWT(:,i),MdenAD(:,i)],'Labels',{'WT','AD'}, ...
        'colors', [0, 0.4470, 0.7410;0.8500, 0.3250, 0.0980]);
    set(h,'linewidth',1)
    xt = get(gca, 'XTick');
    yl=ylim;
    yl=yl(1)+0.95*(yl(2)-yl(1));
    [~,p]=ttest2(MdenWT(:,i),MdenAD(:,i));
    hold on
    if p<0.001
        text(mean(xt([1 2])), yl, '***',HorizontalAlignment='center')
    elseif p<0.01
        text(mean(xt([1 2])), yl, '**',HorizontalAlignment='center')
    elseif p<0.05
        text(mean(xt([1 2])), yl, '*',HorizontalAlignment='center')
    end
    if r == 1, ylabel('Normalized image intensity (a.u.)'), end
    title(names{r},'FontWeight','normal')
    fontsize(10,"points")
    exportpanel(f,outdir,sprintf('Fig6%s_%s_boxplot',char('B'+r-1),names{r})); close(f);
    fprintf('  Fig6 %s: n=%d/%d p=%.3g\n',names{r},sum(~isnan(MdenWT(:,i))),sum(~isnan(MdenAD(:,i))),p);
end
end

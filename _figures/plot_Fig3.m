function plot_Fig3(outdir)
%PLOT_FIG3  Fig 3 panels with the drawing code of vesselsize.m (A/B profiles)
%and newfigures.m (C/D distributionPlot, E boxplotGroup): same colours, plot
%types, line widths and font sizes; only the figure size (Position) is set
%for the assembled figure (see figpos).
%FWHM is converted from voxels to mm with the voxel size of the registered
%grid (S.voxel_mm = 0.030), and the profile x-axis uses the same 30 um spacing.
%The x-axis of C-E is labelled 'Segmentation threshold (top % of brain
%voxels)', on two lines in E.
S = figstyle();
here = fileparts(mfilename('fullpath'));
dd = fullfile(here,'..','Fig1C_Fig3_VesselProfiles','data');
addpath(fullfile(here,'..','_third_party','distributionPlot'), ...
        fullfile(here,'..','_third_party','boxplotGroup'));
C=[0 0.4470 0.7410; 0.8500 0.3250 0.0980; 0.9290 0.6940 0.1250; 0.4940 0.1840 0.5560; 0.4660 0.6740 0.1880];
index=[1 11 16 21 25];
G = {'WT','AD'};

%% A, B: profiles (vesselsize.m, normalized) ----------------------------
for g = 1:2
    gene = G{g};
    load(fullfile(dd,['vsall',gene,'.mat']),'VC');
    time_all=1:15;time_all=(time_all-8)*S.voxel_mm*1000;
    ymean=squeeze(mean(VC,2));
    for i=1:size(VC,1)           % = [ymean(i,:),setting]=mapminmax(ymean(i,:),-1,1); VC(i,:,:)=mapminmax(VC(i,:,:),setting)
        mn=min(ymean(i,:)); mx=max(ymean(i,:));
        ymean(i,:)=2*(ymean(i,:)-mn)/(mx-mn)-1;
        VC(i,:,:)=2*(VC(i,:,:)-mn)/(mx-mn)-1;
    end
    f = figure(Position=figpos(2.10,1.62,11));
    plot(time_all,ymean(index,:)','linewidth',2);
    ystd=squeeze(std(VC,0,2))/sqrt(size(VC,2));
    hold on
    for i=1:5
        h=fill([time_all,time_all(end:-1:1)],[ymean(index(i),:)-ystd(index(i),:),ymean(index(i),end:-1:1)+ystd(index(i),end:-1:1)],C(i,:));
        set(h,'FaceColor',C(i,:),'FaceAlpha',0.2,'EdgeColor','none');
    end
    xlim("tight")
    xlabel(['Radial distance (' char(181) 'm)'])
    ylabel('Normalized signal intensity (a.u.)')
    title(gene)
    fontsize(11,'points')
    leg=legend('50%','30%','20%','10%','2%');set(leg,'Box','off');leg.Location='southeast';
    exportpanel(f,outdir,['Fig3' char('A'+g-1) '_profile_' gene]); close(f);
end

%% C, D: FWHM distributions (newfigures.m, distributionPlot) -------------
for g = 1:2
    gene = G{g};
    load(fullfile(dd,['vsall',gene,'.mat']),'VW');
    VW(VW==0)=nan;
    VW=VW*S.voxel_mm;                       % voxels -> mm
    f = figure(Position=figpos(4.25,2.46,10));
    hd = distributionPlot(reshape(VW,size(VW,1),size(VW,2)*size(VW,3))','histOpt',2); % histOpt=2 works better for uniform distributions than the default
    set(gca,'XTick',1:25)
    set(gca,'XTickLabel',split(num2str(50:-2:2))')
    ylabel('FWHM (mm)','FontSize',12)
    xlabel('Segmentation threshold (top % of brain voxels)')
    title(gene)
    fontsize(10,"points")
    if g == 1                                % legend above the top-right corner of the axes
        ax = gca;
        lg = legend([hd{2}(1) hd{2}(2)],{'mean','median'},'Box','off','Orientation','horizontal');
        lg.Position(1:2) = [ax.Position(1)+ax.Position(3)-lg.Position(3), ax.Position(2)+ax.Position(4)];
    end
    exportpanel(f,outdir,['Fig3' char('C'+g-1) '_FWHM_distribution_' gene]); close(f);
end

%% E: WT vs AD (newfigures.m, boxplotGroup) ------------------------------
load(fullfile(dd,'vsallAD.mat'),'VW');
VW(VW==0)=nan; VW=VW*S.voxel_mm;
ADVW=reshape(VW,[size(VW,1),size(VW,2)*size(VW,3)]);
load(fullfile(dd,'vsallWT.mat'),'VW');
VW(VW==0)=nan; VW=VW*S.voxel_mm;
WTVW=reshape(VW,[size(VW,1),size(VW,2)*size(VW,3)]);
f = figure(Position=figpos(2.10,1.62,10));
h = boxplotGroup({WTVW(index,:)', ADVW(index,:)'}, 'PrimaryLabels', {'',''}, ...
    'SecondaryLabels', {'50%','30%','20%','10%','2%'}, ...
    'groupLabelType', 'horizontal', 'groupLines', true, ...
    'Colors', [0, 0.4470, 0.7410;0.8500, 0.3250, 0.0980], 'GroupType', 'betweenGroups');
% boxplotGroup sets the y-limits to the data range +-5 %; linking its second
% (label) axes resets them to [0 1] in R2023b, so set them again
v = [WTVW(index,:) ADVW(index,:)]; lo = min(v(:)); hi = max(v(:));
ylim(h.axis, [lo hi] + [-1 1]*0.05*(hi-lo));
ylabel(h.axis,'FWHM (mm)')
xl = {'Segmentation threshold','(top % of brain voxels)'};   % two lines: too wide for one at 2.1 in
if isfield(h,'axis2'), xlabel(h.axis2,xl), else, xlabel(xl), end
hold(h.axis,'on')
l1 = plot(h.axis,nan,nan,'-','Color',[0, 0.4470, 0.7410]);
l2 = plot(h.axis,nan,nan,'-','Color',[0.8500, 0.3250, 0.0980]);
fontsize(f,10,"points")
% boxplotGroup fixes the axes position (linked to its label axes): leave room
% for the y label on the left, the two-line x label below and the legend
% above the top-right corner
ax = h.axis;
ax.InnerPosition = [0.20 0.28 0.77 0.58];
lg = legend(h.axis,[l1 l2],{'WT','AD'},'Box','off','Orientation','horizontal');
lg.Position(1:2) = [sum(ax.InnerPosition([1 3]))-lg.Position(3), sum(ax.InnerPosition([2 4]))];
exportpanel(f,outdir,'Fig3E_FWHM_box_WTvsAD'); close(f);
for k = 1:5
    [~,p] = ttest2(WTVW(index(k),:),ADVW(index(k),:));
    fprintf('  Fig3E %2d%%: median WT %.3f AD %.3f mm, t-test over vessels p = %.3g\n', ...
        50-2*(index(k)-1), median(WTVW(index(k),:),'omitnan'), median(ADVW(index(k),:),'omitnan'), p);
end
end

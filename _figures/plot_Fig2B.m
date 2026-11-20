function plot_Fig2B(outdir)
%PLOT_FIG2B  Fig 2B: voxel-intensity histogram of the grayscale-inverted image
%that is thresholded for the segmentation (myseg_tissue_top2to50.m:
%V1 = exp(-(V-B)/1000+1) inside the brain mask, NLM-filtered), WT06, with the
%top 50/30/20/10/2 % thresholds.
%
%The counts/edges are histcounts(raw(mask1==1)) of that script
%(data/Fig2B_histograms_all_animals.mat), i.e. exactly what
%histogram(raw(mask1==1)) draws. Drawing code as the Fig 2B section of
%newfigures.m (black histogram, log y, red dashed lines, linewidth 2); only the
%figure size is set.
here = fileparts(mfilename('fullpath'));
S = load(fullfile(here,'..','Fig2_VesselSegmentation','data','Fig2B_histograms_all_animals.mat'));
H = S.H(find(contains({S.H.file},'WT06'),1));
x = H.thr([1 11 16 21 25]);                   % top 50, 30, 20, 10, 2 %
pct = {'50%','30%','20%','10%','2%'};

f = figure(Position=figpos(3.06,1.17,10));
histogram('BinEdges',H.edges,'BinCounts',H.N,'FaceColor','k');set(gca,'YScale','log');hold on;
% red dashed threshold lines (xline(x,'--r',linewidth=2) in newfigures.m),
% stopped below the labels
ytop = 10^5.6;
for k = 1:5
    plot(x(k)*[1 1],[1 ytop],'--r',linewidth=2);
end
xlabel('Inverted image intensity (a.u.)')
ylabel('Voxels (log)')
fontsize(10,"points")
ylim([1 10^7])
% labels spread out above the lines, each joined to its line by a thin leader
xl = xlim; xs = mean(x) + (-2:2)*0.11*diff(xl);
for k = 1:5
    plot([x(k) xs(k)],[ytop 10^6.05],'-','Color',[0.4 0.4 0.4],'LineWidth',0.5);
    text(xs(k), 10^6.12, pct{k}, 'FontSize',10,'HorizontalAlignment','center','VerticalAlignment','bottom');
end
hold off
fprintf('  Fig2B (%s): thresholds %s\n', H.file, mat2str(x,3));
exportpanel(f,outdir,'Fig2B_histogram_inverted'); close(f);
end

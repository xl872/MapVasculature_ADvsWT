function plot_Fig2B_histogram(subject)
%PLOT_FIG2B_HISTOGRAM  Histogram of the normalised + NLM-filtered image inside
%the brain mask with the segmentation thresholds (top 50% and top 2%) marked,
%from data/Fig2B_histograms_all_animals.mat (one entry per animal; the counts
%and bin edges are exactly histcounts(raw(mask1==1)) of myseg_tissue_top2to50.m,
%thr(k) is the threshold for percent(k), k=1..25 <-> top 50% .. top 2%).
%
%   plot_Fig2B_histogram('AD02')
S = load(fullfile(fileparts(mfilename('fullpath')),'data','Fig2B_histograms_all_animals.mat'));
H = S.H;
i = find(contains({H.file}, subject), 1);
if isempty(i), error('no animal matching %s', subject); end
ctr = (H(i).edges(1:end-1)+H(i).edges(2:end))/2;
figure; plot(ctr, H(i).N, 'k', 'LineWidth', 1.5); hold on
yl = ylim;
for k = [1 numel(H(i).thr)]
    line(H(i).thr(k)*[1 1], yl, 'Color','r','LineStyle','--','LineWidth',1.5);
end
xlabel('Normalized signal intensity'); ylabel('Voxel count');
title(H(i).file, 'Interpreter','none');
end

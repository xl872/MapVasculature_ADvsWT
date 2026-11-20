%% FigS1_plot_from_data.m
% Redraws Fig S1 from data/ (no imaging data needed).
%   data/FigS1_block_means.csv : per animal, 8 block-wise mean intensities before (raw 75 um)
%                               and after SR (37.5 um), as computed by FigS1_raw_vs_SR.m
%   data/FigS1_linear_fit.csv  : per-animal fit SR = a*raw + b
clc; clear; close all;
T = readtable(fullfile('data','FigS1_block_means.csv'));
F = readtable(fullfile('data','FigS1_linear_fit.csv'));
for g = {'AD','WT'}
    Tg = T(strcmp(T.group,g{1}),:);
    animals = unique(Tg.animal,'stable');
    C = jet(256); C = C(round(linspace(1,256,numel(animals))),:);
    figure('Position',[100 100 350 300]); hold on
    for i = 1:numel(animals)
        k = strcmp(Tg.animal,animals{i});
        scatter(Tg.raw_75um_mean(k), Tg.SR_37p5um_mean(k), 40, C(i,:), 'filled');
    end
    xlim([0 10000]); ylim([0 1600]);
    xlabel('Image intensity before super-resolution');
    ylabel('Image intensity after super-resolution');
    title(g{1}); fontsize(12,'points')
end
for v = {'a','b'}
    x = F.(v{1}); grp = F.group;
    figure('Position',[100 100 300 250]);
    boxchart(categorical(grp), x); hold on
    swarmchart(categorical(grp), x, 20, 'filled');
    [~,p] = ttest2(x(strcmp(grp,'WT')), x(strcmp(grp,'AD')));
    ylabel(v{1}); title(sprintf('p = %.3f', p)); fontsize(12,'points')
end

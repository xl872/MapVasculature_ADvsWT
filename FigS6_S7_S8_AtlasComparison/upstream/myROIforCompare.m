%% density
clc;clear all;close all;
% Mean and std of every Allen atlas label (labels.txt) over the 100 permutation
% group maps: Density<tag>_ROI_WTAD.mat from Density<tag>WTAD.mat (vessel
% density, tags 1 6 11 16 21 -1 -6 -11 -16 -21) and Intensity1_ROI_WTAD.mat
% from Density1WTAD_seg_rawz.mat (residual signal intensity). The atlas
% Allen_label_AMBMC30+orig is resampled onto the 450-um block grid of the maps.
% Inputs and outputs in WORK_DIR.
% WORK_DIR = folder holding AD/, WT/ (registered volumes + segmentation .mat),
% the AMBMC/Allen BRIKs (HarvardDataverse/AMBMC) and the intermediate results.
WORK_DIR = '/path/to/vess';   % <-- set this
cd(WORK_DIR);                 % BrikLoad needs the HEAD/BRIK pair in the cwd
for tag=[1 6 11 16 21 -1 -6 -11 -16 -21]
sr=15;
mask= BrikLoad ('ambmcmask+orig.BRIK');
mask=imresize3(mask,1/sr);
mask=imresize3(mask,sr);
mask=permute(mask,[3,2,1]);
[U] = BrikLoad ('Allen_label_AMBMC30+orig.BRIK');
U=permute(U,[3,2,1]);U=imresize3(U,size(mask),"Method","nearest");
U=U(:,1:end-165,:);
mask=mask(:,1:end-165,:);
ROI = readtable('labels.txt');
load(['Density',num2str(tag),'WTAD.mat'])
%%
WT0=cellfun(@(x) permute(imresize3(mean(x, 4),sr),[3,2,1]), WT00, 'UniformOutput', false);
AD0=cellfun(@(x) permute(imresize3(mean(x, 4),sr),[3,2,1]), AD00, 'UniformOutput', false);
%%
alllabel=unique(U);
alllabel=alllabel(2:end);
idxL = arrayfun(@(l) find(U==l), alllabel, 'UniformOutput', false);   % voxel indices of each label, once per tag: x(idx) = x(U==label)
for li=1:numel(alllabel)
    label=alllabel(li); idx=idxL{li};
    Dt=cellfun(@(x) mean(x(idx)), WT0, 'UniformOutput', false);
    ROI{ROI{:, 1} == label, 3} = mean(cell2mat(Dt));
    ROI{ROI{:, 1} == label, 4} = std(cell2mat(Dt),0,'all');
    Dt=cellfun(@(x) mean(x(idx)), AD0, 'UniformOutput', false);
    ROI{ROI{:, 1} == label, 5} = mean(cell2mat(Dt));
    ROI{ROI{:, 1} == label, 6} = std(cell2mat(Dt),0,'all');
end
%%
ROI.Properties.VariableNames = {'label', 'Name', 'WTmean','WTstd','ADmean','ADstd'};
save(['Density',num2str(tag),'_ROI_WTAD'], 'ROI')
end
%% intensity
clc;clear all;close all;
for tag=[1 ]
sr=15;
mask= BrikLoad ('ambmcmask+orig.BRIK');
mask=imresize3(mask,1/sr);
mask=imresize3(mask,sr);
mask=permute(mask,[3,2,1]);
[U] = BrikLoad ('Allen_label_AMBMC30+orig.BRIK');
U=permute(U,[3,2,1]);U=imresize3(U,size(mask),"Method","nearest");
U=U(:,1:end-165,:);
mask=mask(:,1:end-165,:);
ROI = readtable('labels.txt');
load(['Density',num2str(tag),'WTAD_seg_rawz.mat'])
%
WT0=cellfun(@(x) permute(imresize3(mean(x, 4),sr),[3,2,1]), WT00, 'UniformOutput', false);
AD0=cellfun(@(x) permute(imresize3(mean(x, 4),sr),[3,2,1]), AD00, 'UniformOutput', false);
%
alllabel=unique(U);
alllabel=alllabel(2:end);
idxL = arrayfun(@(l) find(U==l), alllabel, 'UniformOutput', false);   % voxel indices of each label, once per tag: x(idx) = x(U==label)
for li=1:numel(alllabel)
    label=alllabel(li); idx=idxL{li};
    Dt=cellfun(@(x) mean(x(idx)), WT0, 'UniformOutput', false);
    ROI{ROI{:, 1} == label, 3} = mean(cell2mat(Dt));
    ROI{ROI{:, 1} == label, 4} = std(cell2mat(Dt),0,'all');
    Dt=cellfun(@(x) mean(x(idx)), AD0, 'UniformOutput', false);
    ROI{ROI{:, 1} == label, 5} = mean(cell2mat(Dt));
    ROI{ROI{:, 1} == label, 6} = std(cell2mat(Dt),0,'all');
end
%
ROI.Properties.VariableNames = {'label', 'Name', 'WTmean','WTstd','ADmean','ADstd'};
save(['Intensity',num2str(tag),'_ROI_WTAD'], 'ROI')   % Fig S7 table; also the template of findROI_intensity_less50.m
end
%% Random_group_result_double_intensity_z_less50.m
% Residual-intensity sensitivity analysis for the 10, 20, 30, 40 and 50%
% compartments produced by myseg_tissue_less50.m.
%
% Vout{k} marks voxels above the kth quantile (k = 1:5 corresponds to
% 10:10:50%).  Those marked voxels are removed, so tag 1 leaves the lowest
% 10% of voxels, tag 2 leaves 20%, and so on.
%
% Outputs: IntensityLeft<percent>_SIvalid5. Blocks with less than 5% retained
% voxels stay NaN rather than becoming a measured zero signal. Each epoch draws
% from all 11 animals of each group.

clc; close all; clear;

%% -------------------------------- settings --------------------------------
thisFile = mfilename('fullpath');
if isempty(thisFile)
    vessDir = pwd;
else
    vessDir = fileparts(thisFile);
end
startDir = pwd;
cleanupDir = onCleanup(@() cd(startDir));
cd(vessDir);

tags             = 1:5;
leftPercents     = 10:10:50;
epochs           = 100;
animalsPerEpoch  = 10;       % 10 of 11 per group, as in Random_group_result_double_intensity_z.m
expectedNWT      = 11;
expectedNAD      = 11;
sr               = 15;
voxelPThreshold  = 0.05;
randomSeed       = 20260817;
minimumBlockValidFraction = 0.05;
analysisSuffix   = '_SIvalid5';

rng(randomSeed,'twister');

%% ------------------------- brain mask and inputs ---------------------------
assertSegmentationCount('AD',expectedNAD);
assertSegmentationCount('WT',expectedNWT);

mask = BrikLoad('ambmcmask+orig.BRIK');
shiftedMask = zeros(size(mask),'like',mask);
shiftedMask(1:339,:,1:243) = mask(1:339,:,7:249);
mask = shiftedMask .* mask;
mask = mask(:,1:end-165,:);
clear shiftedMask

fprintf(['Precomputing z-scored residual-intensity maps with a %.1f%% ' ...
    'minimum block coverage. Missing blocks remain NaN.\n'], ...
    100*minimumBlockValidFraction);
[ADMaps,ADNames] = loadGroupMaps('AD',mask,tags,leftPercents,sr, ...
    expectedNAD,minimumBlockValidFraction);
[WTMaps,WTNames] = loadGroupMaps('WT',mask,tags,leftPercents,sr, ...
    expectedNWT,minimumBlockValidFraction);

nAD = numel(ADNames);
nWT = numel(WTNames);
if animalsPerEpoch > min(nAD,nWT)
    error('IntensityLess50:SampleSize', ...
        'animalsPerEpoch=%d exceeds the available cohort size.',animalsPerEpoch);
end

% Use the same resampled cohorts at every residual threshold. This makes the
% threshold-to-threshold sensitivity comparison paired rather than adding a
% second source of random variation.
ADIndex = makeSelections(nAD,animalsPerEpoch,epochs);
WTIndex = makeSelections(nWT,animalsPerEpoch,epochs);

%% -------------------------- threshold analyses -----------------------------
for t = 1:numel(tags)
    tag = tags(t);
    leftPercent = leftPercents(t);
    fprintf('\nResidual compartment: %d%% left (Vout tag %d)\n',leftPercent,tag);

    ADMouse = bilateralMouseMaps(ADMaps{t});
    WTMouse = bilateralMouseMaps(WTMaps{t});

    WT00 = cell(1,epochs);
    AD00 = cell(1,epochs);
    lowSize = size(WTMaps{t},1:3);
    map  = zeros([epochs lowSize],'single');
    mapt = zeros([epochs lowSize],'single');
    p    = ones([epochs lowSize],'single');

    for epoch = 1:epochs
        AD0 = selectedHemispheres(ADMaps{t},ADIndex{epoch});
        WT0 = selectedHemispheres(WTMaps{t},WTIndex{epoch});
        AD00{epoch} = AD0;
        WT00{epoch} = WT0;

        [~,pVDV] = ttest2(AD0,WT0,'Dim',4);
        pVDV(~isfinite(pVDV)) = 1;
        groupDifference = mean(WT0,4,'omitnan') - mean(AD0,4,'omitnan');
        significant = pVDV < voxelPThreshold & isfinite(groupDifference);
        significantDifference = zeros(size(groupDifference),'like',groupDifference);
        significantDifference(significant) = groupDifference(significant);

        map(epoch,:,:,:)  = single(groupDifference);
        mapt(epoch,:,:,:) = single(significantDifference);
        p(epoch,:,:,:)    = single(pVDV);

        if epoch == 1 || mod(epoch,10) == 0
            fprintf('  epoch %3d/%d\n',epoch,epochs);
        end
    end

    [Mh,ph] = ttest(mapt,0,'Dim',1);
    Mh = squeeze(Mh);
    ph = squeeze(ph);
    ph(~isfinite(ph)) = 1;

    outStem = sprintf('IntensityLeft%d%s',leftPercent,analysisSuffix);
    save([outStem 'WTAD.mat'], ...
        'WT00','AD00','WTMouse','ADMouse','WTNames','ADNames', ...
        'WTIndex','ADIndex','tag','leftPercent','epochs','animalsPerEpoch', ...
        'sr','randomSeed','minimumBlockValidFraction','analysisSuffix','-v7.3');
    save([outStem '.mat'], ...
        'map','mapt','p','tag','leftPercent','sr','voxelPThreshold', ...
        'minimumBlockValidFraction','analysisSuffix','-v7.3');
    save([outStem 'final.mat'], ...
        'Mh','ph','tag','leftPercent','sr','voxelPThreshold', ...
        'minimumBlockValidFraction','analysisSuffix','-v7.3');

    fprintf('  saved %sWTAD.mat, %s.mat and %sfinal.mat\n', ...
        outStem,outStem,outStem);
    clear WT00 AD00 WTMouse ADMouse map mapt p Mh ph
end

fprintf('\nAll residual-intensity thresholds completed.\n');

%% -------------------------------- functions --------------------------------
function assertSegmentationCount(groupName,expectedN)
% Stop before loading the 53-MB mask or any 330-MB image if a cohort is still
% being segmented.
    files = dir(fullfile(groupName,'*resultinrangeless50.mat'));
    files = files(~[files.isdir]);
    if numel(files) ~= expectedN
        error('IntensityLess50:IncompleteSegmentation', ...
            ['Expected %d %s *resultinrangeless50.mat files but found %d. ' ...
             'Wait for myseg_tissue_less50.m to finish before running this script.'], ...
             expectedN,groupName,numel(files));
    end
end

function [maps,names] = loadGroupMaps(groupName,mask,tags,leftPercents,sr, ...
        expectedN,minimumBlockValidFraction)
% Load one Vout cell at a time from each v7.3 segmentation file. Loading the
% complete five-cell Vout variable would require about 1.65 GB per mouse.
    files = dir(fullfile(groupName,'*resultinrangeless50.mat'));
    files = files(~[files.isdir]);
    [~,order] = sort({files.name});
    files = files(order);

    if numel(files) ~= expectedN
        error('IntensityLess50:IncompleteSegmentation', ...
            ['Expected %d %s *resultinrangeless50.mat files but found %d. ' ...
             'Wait for myseg_tissue_less50.m to finish before running this script.'], ...
             expectedN,groupName,numel(files));
    end

    names = {files.name}';
    maps = cell(1,numel(tags));
    for animal = 1:numel(files)
        segFile = fullfile(files(animal).folder,files(animal).name);
        seg = matfile(segFile);
        percent = double(seg.percent);
        if numel(percent) < max(tags) || ...
                any(abs(percent(tags) - leftPercents/100) > 1e-10)
            error('IntensityLess50:ThresholdMismatch', ...
                '%s does not contain the expected 10:10:50%% thresholds.',segFile);
        end

        rawFile = rawImageForSegmentation(files(animal));
        raw = BrikLoad(rawFile);
        raw = raw(:,1:end-165,:);
        if ~isequal(size(raw),size(mask))
            error('IntensityLess50:RawGridMismatch', ...
                'Raw image %s has size %s; mask size is %s.', ...
                rawFile,mat2str(size(raw)),mat2str(size(mask)));
        end

        brainValues = raw(mask == 1);
        rawSD = std(brainValues,0);
        if ~isfinite(rawSD) || rawSD == 0
            error('IntensityLess50:InvalidRawScale', ...
                'Cannot z-score %s because its brain-mask SD is zero/nonfinite.',rawFile);
        end
        raw = (raw - mean(brainValues)) ./ rawSD;
        clear brainValues

        fprintf('  %s mouse %2d/%d: %s\n',groupName,animal,numel(files),files(animal).name);
        for t = 1:numel(tags)
            oneCell = seg.Vout(1,tags(t));
            Tout = oneCell{1};
            if ~isequal(size(Tout),size(mask))
                error('IntensityLess50:SegmentationGridMismatch', ...
                    'Vout{%d} in %s has size %s; mask size is %s.', ...
                    tags(t),segFile,mat2str(size(Tout)),mat2str(size(mask)));
            end
            small = residualIntensity(mask,Tout,raw,sr,minimumBlockValidFraction);
            if animal == 1
                maps{t} = zeros([size(small) numel(files)],'like',small);
            end
            maps{t}(:,:,:,animal) = small;
            clear oneCell Tout small
        end
        clear raw seg
    end
end

function rawFile = rawImageForSegmentation(fileEntry)
% The first 19 characters are AZ_Vasculature_<mouse ID> for every source
% naming variant in this data set (including WT02.sr and the 38-um scans).
    if numel(fileEntry.name) < 19
        error('IntensityLess50:UnexpectedFilename','Unexpected file name: %s',fileEntry.name);
    end
    subjectStem = fileEntry.name(1:19);
    if contains(fileEntry.name,'38um')
        rawName = [subjectStem '_38um.nii_reg1+orig.BRIK'];
    else
        rawName = [subjectStem '_75um.nii_regraw+orig.BRIK'];
    end
    rawFile = fullfile(fileEntry.folder,rawName);
    if ~exist(rawFile,'file')
        error('IntensityLess50:MissingRawImage','Missing raw image %s.',rawFile);
    end
end

function out = residualIntensity(mask,Tout,raw,sr,minimumBlockValidFraction)
% Remove the vessel-weighted/high-intensity compartment and average only the
% residual voxels in each sr-by-sr-by-sr block.
    out = raw;
    out(Tout == 1 | mask == 0) = NaN;
    out = mysizedown_SI(out,sr,minimumBlockValidFraction);
end

function selections = makeSelections(nAnimals,nSelected,nEpoch)
    selections = cell(1,nEpoch);
    for e = 1:nEpoch
        selections{e} = randperm(nAnimals,nSelected);
    end
end

function hemispheres = selectedHemispheres(animalMaps,selection)
% Put the map and its left-right flipped copy next to one another for every
% selected mouse, as in Random_group_result_double_intensity_z.m.
    nSelected = numel(selection);
    mapSize = size(animalMaps,1:3);
    hemispheres = zeros([mapSize 2*nSelected],'like',animalMaps);
    for q = 1:nSelected
        one = animalMaps(:,:,:,selection(q));
        hemispheres(:,:,:,2*q-1) = one;
        hemispheres(:,:,:,2*q) = flip(one,1);
    end
end

function mouseMaps = bilateralMouseMaps(animalMaps)
% One independent, atlas-aligned observation per mouse for ROI statistics.
% If only one mirrored side has adequate support, retain that side rather than
% turning the bilateral value into NaN.
    mouseMaps = mean(cat(5,animalMaps,flip(animalMaps,1)),5,'omitnan');
end

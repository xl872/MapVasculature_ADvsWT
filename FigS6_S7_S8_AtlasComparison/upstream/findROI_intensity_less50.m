%% findROI_intensity_less50.m
% Batch ROI processing for the residual-intensity sensitivity analysis.
%
% Run Random_group_result_double_intensity_z_less50.m first. For each 10,
% 20, 30, 40 and 50% residual compartment, this script:
%   1. measures every mouse separately in every Allen atlas label;
%   2. calculates WT/AD group mean and SD;
%   3. fills atlas parent structures using the same hierarchy/voxel weights as
%      the 50% analysis (Intensity1_ROI_WTAD_full); and
%   4. writes direct and hierarchy-filled MAT/CSV tables (Fig S8).
%
% It uses NaN-aware resizing and valid-voxel weights, requires at least 5% valid
% support in an atlas ROI for each animal and does not turn missing
% measurements into zero.
%
% Unlike findROI_intensity.m, this is a batch step: it does not require
% manually choosing one slice or drawing a single ROI.

clc; clear; close all;

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

leftPercents = 10:10:50;
sr = 15;
cropTail = 165;
analysisSuffix = '_SIvalid5';
minimumROIValidFraction = 0.05;

labelFile = fullfile(vessDir,'labels.txt');
atlasFile = fullfile(vessDir,'Allen_label_AMBMC30+orig.BRIK');
maskFile = fullfile(vessDir,'ambmcmask+orig.BRIK');
templateFile = fullfile(vessDir,'Intensity1_ROI_WTAD_full.mat');

if ~exist(templateFile,'file')
    error('IntensityLess50:MissingHierarchyTemplate', ...
        ['Missing %s. Run myAggregateHigherROIs.m once for the original ' ...
         'Intensity1 analysis before this script.'],templateFile);
end

%% ----------------------- atlas and hierarchy template ----------------------
fprintf('Building the atlas grid used for the original intensity ROI table ...\n');
mask = BrikLoad(maskFile);
mask = imresize3(mask,1/sr);
mask = imresize3(mask,sr);
mask = permute(mask,[3 2 1]);

U = BrikLoad(atlasFile);
U = permute(U,[3 2 1]);
U = imresize3(U,size(mask),'Method','nearest');
U = U(:,1:end-cropTail,:);
clear mask

atlasSub = double(U(:)) + 1;
nAtlasLabel = max(atlasSub);
atlasCount = accumarray(atlasSub,1,[nAtlasLabel 1],@sum,0);
clear U

baseROI = readtable(labelFile);
baseROI = baseROI(:,1:2);
baseROI.Properties.VariableNames = {'label','Name'};
nBaseROI = height(baseROI);

T = load(templateFile,'ROI');
templateROI = T.ROI;
clear T
requiredTemplateColumns = {'Nvox','SubtreeVox','IsAggregated'};
if ~all(ismember(requiredTemplateColumns,templateROI.Properties.VariableNames))
    error('IntensityLess50:InvalidHierarchyTemplate', ...
        '%s lacks Nvox, SubtreeVox or IsAggregated.',templateFile);
end
if height(templateROI) < nBaseROI || ...
        ~isequal(double(templateROI.label(1:nBaseROI)),double(baseROI.label))
    error('IntensityLess50:HierarchyOrderMismatch', ...
        'The hierarchy template and labels.txt are not in the same label order.');
end

%% ---------------------------- process thresholds ---------------------------
for leftPercent = leftPercents
    inputFile = fullfile(vessDir,sprintf('IntensityLeft%d%sWTAD.mat', ...
        leftPercent,analysisSuffix));
    if ~exist(inputFile,'file')
        error('IntensityLess50:MissingRandomGroupResult', ...
            ['Missing %s. Finish the segmentation and run ' ...
             'Random_group_result_double_intensity_z_less50.m first.'],inputFile);
    end

    fprintf('\nMeasuring atlas ROIs for the %d%% residual compartment ...\n',leftPercent);
    S = load(inputFile,'WTMouse','ADMouse','WTNames','ADNames','sr', ...
        'leftPercent','minimumBlockValidFraction');
    if ~isfield(S,'WTMouse') || ~isfield(S,'ADMouse')
        error('IntensityLess50:MissingMouseMaps', ...
            '%s does not contain WTMouse and ADMouse.',inputFile);
    end
    if S.sr ~= sr
        error('IntensityLess50:ScaleMismatch', ...
            '%s uses sr=%d, but this script is configured for sr=%d.',inputFile,S.sr,sr);
    end

    [WTAnimalDirect,WTAnimalDirectWeight] = atlasAnimalValues( ...
        S.WTMouse,baseROI.label,atlasSub,atlasCount,sr,'WT', ...
        minimumROIValidFraction);
    [ADAnimalDirect,ADAnimalDirectWeight] = atlasAnimalValues( ...
        S.ADMouse,baseROI.label,atlasSub,atlasCount,sr,'AD', ...
        minimumROIValidFraction);

    % Direct-label table (same columns as Intensity1_ROI_WTAD).
    ROI = baseROI;
    [ROI.WTmean,ROI.WTstd,ROI.WTn] = summarizeAnimals(WTAnimalDirect);
    [ROI.ADmean,ROI.ADstd,ROI.ADn] = summarizeAnimals(ADAnimalDirect);

    directBase = fullfile(vessDir,sprintf( ...
        'IntensityLeft%d%s_ROI_WTAD',leftPercent,analysisSuffix));
    WTAnimalROI = WTAnimalDirect;
    ADAnimalROI = ADAnimalDirect;
    WTAnimalValidWeight = WTAnimalDirectWeight;
    ADAnimalValidWeight = ADAnimalDirectWeight;
    WTNames = S.WTNames;
    ADNames = S.ADNames;
    save([directBase '.mat'],'ROI','WTAnimalROI','ADAnimalROI', ...
        'WTAnimalValidWeight','ADAnimalValidWeight','WTNames','ADNames', ...
        'leftPercent','sr','analysisSuffix','minimumROIValidFraction');
    writetable(ROI,[directBase '.csv']);

    % Hierarchy-filled table. The template records exactly which parent rows
    % lacked a direct atlas value in the finalized 50% analysis.
    [WTAnimalROI,WTAnimalValidWeight] = fillHierarchy( ...
        WTAnimalDirect,WTAnimalDirectWeight,templateROI,nBaseROI);
    [ADAnimalROI,ADAnimalValidWeight] = fillHierarchy( ...
        ADAnimalDirect,ADAnimalDirectWeight,templateROI,nBaseROI);
    ROI = templateROI;
    [ROI.WTmean,ROI.WTstd,ROI.WTn] = summarizeAnimals(WTAnimalROI);
    [ROI.ADmean,ROI.ADstd,ROI.ADn] = summarizeAnimals(ADAnimalROI);

    fullBase = fullfile(vessDir,sprintf( ...
        'IntensityLeft%d%s_ROI_WTAD_full',leftPercent,analysisSuffix));
    save([fullBase '.mat'],'ROI','WTAnimalROI','ADAnimalROI', ...
        'WTAnimalValidWeight','ADAnimalValidWeight','WTNames','ADNames', ...
        'leftPercent','sr','analysisSuffix','minimumROIValidFraction');
    writetable(ROI,[fullBase '.csv']);

    fprintf('  saved %s(.mat/.csv) and %s(.mat/.csv)\n',directBase,fullBase);
    clear S ROI WTAnimalDirect ADAnimalDirect WTAnimalROI ADAnimalROI
end

fprintf('\nAll 10-50%% residual-intensity ROI tables completed.\n');

%% -------------------------------- functions --------------------------------
function [values,validWeights] = atlasAnimalValues(mouseMaps,roiLabels, ...
        atlasSub,atlasCount,sr,groupName,minimumROIValidFraction)
% One atlas-label mean per independent mouse. Mouse maps have already been
% averaged across their original/flipped hemispheres by the random-group step.
    nMouse = size(mouseMaps,4);
    nROI = numel(roiLabels);
    roiSub = double(roiLabels) + 1;
    labelExists = roiSub >= 1 & roiSub <= numel(atlasCount);
    candidate = find(labelExists);
    labelExists(candidate) = atlasCount(roiSub(candidate)) > 0;
    values = nan(nMouse,nROI);
    validWeights = zeros(nMouse,nROI);

    fprintf('  %s mice',groupName);
    for m = 1:nMouse
        [high,support] = resizeFiniteMean(mouseMaps(:,:,:,m),sr);
        high = permute(high,[3 2 1]);
        support = permute(support,[3 2 1]);
        if numel(high) ~= numel(atlasSub)
            error('IntensityLess50:AtlasGridMismatch', ...
                '%s mouse %d map size %s does not match the atlas grid.', ...
                groupName,m,mat2str(size(high)));
        end

        finite = isfinite(high(:)) & support(:) > 0;
        labelWeight = accumarray(atlasSub(finite),double(support(finite)), ...
            size(atlasCount),@sum,0);
        labelSum = accumarray(atlasSub(finite), ...
            double(high(finite)).*double(support(finite)), ...
            size(atlasCount),@sum,0);
        labelMean = labelSum ./ labelWeight;
        labelCoverage = labelWeight ./ atlasCount;
        labelMean(labelWeight == 0 | ...
            labelCoverage < minimumROIValidFraction) = NaN;
        values(m,labelExists) = labelMean(roiSub(labelExists))';
        validWeights(m,labelExists) = labelWeight(roiSub(labelExists))';
        validWeights(m,~isfinite(values(m,:))) = 0;
        fprintf('.');
    end
    fprintf(' done (n=%d)\n',nMouse);
end

function [high,support] = resizeFiniteMean(low,sr)
% Resize values and their validity mask separately. Dividing the resized
% numerator by resized support prevents NaNs from spreading to adjacent blocks.
    valid = isfinite(low);
    numerator = double(low);
    numerator(~valid) = 0;
    numerator = imresize3(numerator,sr,'linear');
    support = imresize3(double(valid),sr,'linear');
    high = numerator ./ support;
    high(support <= eps) = NaN;
end

function [fullValues,fullWeights] = fillHierarchy( ...
        directValues,directWeights,templateROI,nBaseROI)
% Reproduce the hierarchy fill used by myAggregateHigherROIs.m, including
% spelling aliases appended after the base labels.txt rows.
    nMouse = size(directValues,1);
    nROI = height(templateROI);
    fullValues = nan(nMouse,nROI);
    fullWeights = zeros(nMouse,nROI);
    fullValues(:,1:nBaseROI) = directValues;
    fullWeights(:,1:nBaseROI) = directWeights;

    for r = 1:nBaseROI
        if ~templateROI.IsAggregated(r)
            continue;
        end

        targetVoxels = double(templateROI.SubtreeVox(r));
        cumulativeVoxels = cumsum(double(templateROI.Nvox(r:nBaseROI)));
        relativeEnd = find(cumulativeVoxels >= targetVoxels,1);
        if isempty(relativeEnd) || cumulativeVoxels(relativeEnd) ~= targetVoxels
            error('IntensityLess50:SubtreeMismatch', ...
                'Could not recover the descendant block for ROI %s.',string(templateROI.Name(r)));
        end

        rows = r:(r+relativeEnd-1);
        atlasWeights = double(templateROI.Nvox(rows))';
        useWeight = atlasWeights > 0;
        rows = rows(useWeight);

        v = directValues(:,rows);
        weights = directWeights(:,rows);
        weights(~isfinite(v)) = 0;
        v(~isfinite(v)) = 0;
        numerator = sum(v .* weights,2);
        denominator = sum(weights,2);
        fullValues(:,r) = numerator ./ denominator;
        fullValues(denominator == 0,r) = NaN;
        fullWeights(:,r) = denominator;
    end

    % myAggregateHigherROIs.m appends spelling aliases after labels.txt.
    for r = nBaseROI+1:nROI
        source = find(double(templateROI.label(1:nBaseROI)) == ...
            double(templateROI.label(r)),1);
        if isempty(source)
            error('IntensityLess50:AliasMismatch', ...
                'No base ROI found for alias %s.',string(templateROI.Name(r)));
        end
        fullValues(:,r) = fullValues(:,source);
        fullWeights(:,r) = fullWeights(:,source);
    end
end

function [groupMean,groupSD,groupN] = summarizeAnimals(animalValues)
    groupMean = mean(animalValues,1,'omitnan')';
    groupSD = std(animalValues,0,1,'omitnan')';
    groupN = sum(isfinite(animalValues),1)';
    groupMean(groupN == 0) = NaN;
    groupSD(groupN < 2) = NaN;
end

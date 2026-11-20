%% myPlotFigureS6SlopesSTE.m
% Fig S6 F: slope of the pooled regional fit (vessel density against the
% vessel length density of the Ji et al. 2021 atlas) at the top 50, 40, 30, 20
% and 10 % thresholds, with an across-mouse standard error (STE/SEM).
%
% The marker centers are intentionally the same pooled regional fits used in
% myComparePlotDensity_marked.m: one fit of the group-average density across
% 345 ROIs at each threshold.  They are NOT replaced by averages of animal
% slopes.  The error bars are calculated separately by fitting the same
% regional relationship for each mouse and then taking
%
%                  STE = std(mouse slopes) / sqrt(number of mice).
%
% Density<tag>WTAD.mat contains left and right hemispheres in adjacent pairs.
% The two hemispheres are averaged first so that the independent unit for the
% error bar is a mouse, not a hemisphere or one of the 100 resampling repeats.
% The repeats store the same mouse at a fixed pair of indices; a mouse left out
% of a repeat is represented by an all-zero pair.  recoverMouseMaps() takes the
% first nonzero occurrence of every mouse and therefore counts each mouse once.
%
% Outputs (in this script's directory):
%   FigureS6_pooled_slopes_with_STE.png/.pdf/.eps
%   FigureS6_pooled_slopes_with_STE.csv
%   FigureS6_pooled_slopes_with_STE.mat

clc; clear; close all;

%% ------------------------------- settings -------------------------------
thisDir = fileparts(mfilename('fullpath'));
vessDir = fullfile(thisDir,'data');      % small inputs in data/
ATLAS_DIR = '/path/to/HarvardDataverse/AMBMC';   % <-- set this: Allen_label_AMBMC30+orig.* and ambmcmask+orig.*

xlsxFile  = fullfile(vessDir,'1-s2.0-S0896627321000805-mmc2.xlsx');
labelFile = fullfile(vessDir,'labels.txt');
atlasFile = fullfile(ATLAS_DIR,'Allen_label_AMBMC30+orig.BRIK');
maskFile  = fullfile(ATLAS_DIR,'ambmcmask+orig.BRIK');

tags       = [1 6 11 16 21];             % top 50, 40, 30, 20 and 10 percent
thresholds = [50 40 30 20 10];
sr         = 15;                          % same as myROIforCompare.m
cropTail   = 165;                         % same as myROIforCompare.m

expectedNWT = 12;                         % a warning is printed if the data hold another count
expectedNAD = 11;

cWT = [0.0000 0.4470 0.7410];
cAD = [0.8500 0.3250 0.0980];

%% -------------------- reference values and atlas grid -------------------
CP = readtable(xlsxFile);

% Build the atlas exactly as myROIforCompare.m does.  atlasSub is a linear
% label index used by accumarray, allowing all ROI means to be obtained in one
% pass through each upsampled mouse map.
fprintf('Building the atlas grid used by myROIforCompare.m ...\n');
% The installed AFNI MATLAB reader only recognizes a .HEAD/.BRIK pair when
% the file is in MATLAB's current folder, even if an absolute path is passed.
% Read from ATLAS_DIR, then return to the caller's folder.
startDir = pwd;
cd(ATLAS_DIR);
mask = BrikLoad(maskFile);
mask = imresize3(mask,1/sr);
mask = imresize3(mask,sr);
mask = permute(mask,[3 2 1]);

U = BrikLoad(atlasFile);
cd(startDir);
U = permute(U,[3 2 1]);
U = imresize3(U,size(mask),'Method','nearest');
U = U(:,1:end-cropTail,:);
clear mask

atlasSub = double(U(:)) + 1;              % accumarray indices must start at 1
nAtlasLabel = max(atlasSub);
atlasCount = accumarray(atlasSub,1,[nAtlasLabel 1],@sum,0);
clear U

LAB = readtable(labelFile);
nBaseROI = height(LAB);                    % excludes spelling aliases appended later
clear LAB

%% ---------------- pooled and per-mouse slope at every threshold ---------
nTag = numel(tags);
pooledWT = nan(nTag,1); pooledAD = nan(nTag,1);
meanMouseWT = nan(nTag,1); meanMouseAD = nan(nTag,1);
steWT = nan(nTag,1); steAD = nan(nTag,1);
nWT = nan(nTag,1); nAD = nan(nTag,1);
nRoiWT = nan(nTag,1); nRoiAD = nan(nTag,1);
mouseSlopesWT = cell(nTag,1); mouseSlopesAD = cell(nTag,1);

for t = 1:nTag
    tag = tags(t);
    fprintf('\nTag %d (top %d%%)\n',tag,thresholds(t));

    % ----- the pooled point: identical definition to the marked script -----
    roiFile = fullfile(vessDir,sprintf('Density%d_ROI_WTAD_full.mat',tag));
    S = load(roiFile,'ROI');
    ROI = S.ROI;
    ROI.SourceRow = (1:height(ROI))';
    mergedTable = innerjoin(CP,ROI,'LeftKeys','Acronym','RightKeys','Name');

    x = str2double(string(mergedTable.ThisStudy));
    validWT = isfinite(x) & isfinite(mergedTable.WTmean) & mergedTable.WTmean ~= 0;
    validAD = isfinite(x) & isfinite(mergedTable.ADmean) & mergedTable.ADmean ~= 0;
    nRoiWT(t) = sum(validWT);
    nRoiAD(t) = sum(validAD);

    fitWT = polyfit(x(validWT),mergedTable.WTmean(validWT),1);
    fitAD = polyfit(x(validAD),mergedTable.ADmean(validAD),1);
    pooledWT(t) = fitWT(1);
    pooledAD(t) = fitAD(1);

    % ----- recover one bilateral density map per independent mouse ----------
    % Density%dWTAD_mice.mat: the group maps of each mouse (left-right pairs), as in
    % Density%dWTAD.mat; recoverMouseMaps() takes each mouse once.
    densityFile = fullfile(vessDir,sprintf('Density%dWTAD_mice.mat',tag));
    D = load(densityFile,'WT00','AD00');
    wtMaps = recoverMouseMaps(D.WT00,'WT');
    adMaps = recoverMouseMaps(D.AD00,'AD');
    clear D

    nWT(t) = size(wtMaps,4);
    nAD(t) = size(adMaps,4);

    % One value per atlas ROI per mouse, including the same hierarchy-filled
    % parent regions used in Density<tag>_ROI_WTAD_full.mat.
    animalRoiWT = mouseRoiValues(wtMaps,ROI,nBaseROI,atlasSub,atlasCount,sr,'WT');
    animalRoiAD = mouseRoiValues(adMaps,ROI,nBaseROI,atlasSub,atlasCount,sr,'AD');
    clear wtMaps adMaps

    row = double(mergedTable.SourceRow);
    mouseSlopesWT{t} = fitMouseSlopes(x,validWT,animalRoiWT(:,row));
    mouseSlopesAD{t} = fitMouseSlopes(x,validAD,animalRoiAD(:,row));
    clear animalRoiWT animalRoiAD

    meanMouseWT(t) = mean(mouseSlopesWT{t});
    meanMouseAD(t) = mean(mouseSlopesAD{t});
    steWT(t) = std(mouseSlopesWT{t},0) / sqrt(numel(mouseSlopesWT{t}));
    steAD(t) = std(mouseSlopesAD{t},0) / sqrt(numel(mouseSlopesAD{t}));

    fprintf(['  pooled slopes: WT %.6f, AD %.6f (%d/%d ROIs)\n' ...
             '  mouse slopes : WT mean %.6f, STE %.6f (n=%d); ' ...
             'AD mean %.6f, STE %.6f (n=%d)\n'], ...
             pooledWT(t),pooledAD(t),nRoiWT(t),nRoiAD(t), ...
             meanMouseWT(t),steWT(t),nWT(t), ...
             meanMouseAD(t),steAD(t),nAD(t));
end

% Warn if the number of mice differs from expectedNWT / expectedNAD; the STE
% uses the mice actually present.
if any(nWT ~= expectedNWT)
    warning('FigureS6:WTCount', ...
        ['The response document states %d WT mice, but the saved Density<tag>WTAD ' ...
         'files contain %s. The plotted WT STE uses the mice actually present.'], ...
         expectedNWT,mat2str(unique(nWT)'));
end
if any(nAD ~= expectedNAD)
    warning('FigureS6:ADCount', ...
        ['The response document states %d AD mice, but the saved Density<tag>WTAD ' ...
         'files contain %s. The plotted AD STE uses the mice actually present.'], ...
         expectedNAD,mat2str(unique(nAD)'));
end
if any(nRoiWT ~= 345) || any(nRoiAD ~= 345)
    warning('FigureS6:ROICount', ...
        'At least one pooled fit does not contain the expected 345 ROIs.');
end

%% -------------------------------- figure --------------------------------
fig = figure('Color','w','Position',[100 100 700 400]);
ax = axes(fig,'Position',[0.10 0.18 0.87 0.72]);
hold(ax,'on');

hWT = errorbar(ax,thresholds,pooledWT,steWT,'-o', ...
    'Color',cWT,'MarkerFaceColor',cWT,'MarkerEdgeColor',cWT, ...
    'MarkerSize',9,'LineWidth',2.2,'CapSize',12);
hAD = errorbar(ax,thresholds,pooledAD,steAD,'-s', ...
    'Color',cAD,'MarkerFaceColor',cAD,'MarkerEdgeColor',cAD, ...
    'MarkerSize',9,'LineWidth',2.2,'CapSize',12);

% Numerical labels (not drawn): placed beyond the STE bars so the values and
% uncertainty do not overlap.
labelGap = 0.0035;
% for t = 1:nTag
%     text(ax,thresholds(t),pooledWT(t)+steWT(t)+labelGap, ...
%         sprintf('%.3f',pooledWT(t)),'Color',cWT,'FontSize',11, ...
%         'HorizontalAlignment','center','VerticalAlignment','bottom');
%     text(ax,thresholds(t),pooledAD(t)-steAD(t)-labelGap, ...
%         sprintf('%.3f',pooledAD(t)),'Color',cAD,'FontSize',11, ...
%         'HorizontalAlignment','center','VerticalAlignment','top');
% end

hold(ax,'off');
set(ax,'XDir','reverse','XTick',sort(thresholds),'XLim',[5 55], ...
    'FontSize',14,'LineWidth',1,'Box','off','Layer','top');
ylimTop = max(0.14,max([pooledWT+steWT; pooledAD+steAD])+0.012);
ylim(ax,[0 ylimTop]);
grid(ax,'off');
ax.GridColor = [0.82 0.84 0.87];
ax.GridAlpha = 0.7;
ax.XMinorGrid = 'off';

% title(ax,'Pooled regional regression slope by segmentation threshold', ...
%     'FontSize',20,'FontWeight','normal');
xlabel(ax,'Segmentation threshold (top voxel fraction, %)','FontSize',17);
ylabel(ax,'Regression slope across 345 ROIs','FontSize',17);
legend(ax,[hWT hAD],{'WT','AD'},'Location','northeast','Box','off','FontSize',16);

% note = sprintf(['Points: pooled 345-ROI fits. Error bars: STE of mouse-level ' ...
%                 'slopes (WT n=%s; AD n=%s).'], ...
%                 countLabel(nWT),countLabel(nAD));
% annotation(fig,'textbox',[0.08 0.025 0.84 0.065], ...
%     'String',note,'EdgeColor','none','HorizontalAlignment','center', ...
%     'VerticalAlignment','middle','FontSize',12,'Color',[0.30 0.33 0.38]);

outBase = fullfile(thisDir,'FigureS6_pooled_slopes_with_STE');
exportgraphics(fig,[outBase '.png'],'Resolution',300);
exportgraphics(fig,[outBase '.pdf'],'ContentType','vector');
exportgraphics(fig,[outBase '.eps'],'ContentType','vector');

%% ------------------------- results table and cache -----------------------
ThresholdTopPercent = thresholds';
Tag = tags';
WT_PooledSlope = pooledWT;
WT_MouseMeanSlope = meanMouseWT;
WT_STE = steWT;
WT_N = nWT;
WT_NROI = nRoiWT;
AD_PooledSlope = pooledAD;
AD_MouseMeanSlope = meanMouseAD;
AD_STE = steAD;
AD_N = nAD;
AD_NROI = nRoiAD;

Results = table(ThresholdTopPercent,Tag, ...
    WT_PooledSlope,WT_MouseMeanSlope,WT_STE,WT_N,WT_NROI, ...
    AD_PooledSlope,AD_MouseMeanSlope,AD_STE,AD_N,AD_NROI);
writetable(Results,[outBase '.csv']);
save([outBase '.mat'],'Results','mouseSlopesWT','mouseSlopesAD');

fprintf('\nSaved:\n  %s.png/.pdf/.eps\n  %s.csv/.mat\n',outBase,outBase);
disp(Results);

%% ------------------------------ functions -------------------------------
function mouseMaps = recoverMouseMaps(repeats,groupName)
% Recover each mouse once from fixed left/right positions across repeats.
    nRepeat = numel(repeats);
    nMouse = 0;
    for e = 1:nRepeat
        nMouse = max(nMouse,floor(size(repeats{e},4)/2));
    end

    mapSize = size(repeats{1});
    mouseMaps = zeros([mapSize(1:3) nMouse],'like',repeats{1});
    found = false(1,nMouse);

    for e = 1:nRepeat
        A = repeats{e};
        for m = find(~found)
            j = 2*m-1;
            if size(A,4) < j+1
                continue
            end
            left  = A(:,:,:,j);
            right = A(:,:,:,j+1);

            % Dynamic assignment in the generating script fills an omitted
            % mouse's pair with zeros.  A real whole-brain density map always
            % contains at least one nonzero voxel.
            if any(left(:) ~= 0) || any(right(:) ~= 0)
                mouseMaps(:,:,:,m) = (left + right)/2;
                found(m) = true;
            end
        end
        if all(found), break; end
    end

    if ~all(found)
        error('FigureS6:MissingMouse', ...
            '%s mouse indices not recoverable from the repeats: %s', ...
            groupName,mat2str(find(~found)));
    end
    fprintf('  recovered %d %s mice (bilateral maps)\n',nMouse,groupName);
end

function animalROI = mouseRoiValues(mouseMaps,ROI,nBaseROI, ...
                                    atlasSub,atlasCount,sr,groupName)
% Direct atlas-label means followed by the same parent-region aggregation as
% myAggregateHigherROIs.m.  Directly labelled parents keep their own label;
% parents with no direct voxels use the voxel-weighted descendant mean.
    nMouse = size(mouseMaps,4);
    nROI = height(ROI);
    roiLabelSub = double(ROI.label) + 1;
    labelIsInAtlas = roiLabelSub <= numel(atlasCount);
    direct = nan(nMouse,nROI);

    fprintf('  measuring %s mouse ROIs',groupName);
    for m = 1:nMouse
        high = imresize3(mouseMaps(:,:,:,m),sr);
        high = permute(high,[3 2 1]);
        if numel(high) ~= numel(atlasSub)
            error('FigureS6:GridMismatch', ...
                'Mouse map %s does not match the atlas grid.',mat2str(size(high)));
        end

        labelSum = accumarray(atlasSub,double(high(:)), ...
                              size(atlasCount),@sum,0);
        labelMean = labelSum ./ atlasCount;
        labelMean(atlasCount == 0) = NaN;
        direct(m,labelIsInAtlas) = labelMean(roiLabelSub(labelIsInAtlas))';
        fprintf('.');
    end
    fprintf(' done\n');

    animalROI = direct;
    nBaseROI = min(nBaseROI,nROI);

    % labels.txt and the base ROI table are in depth-first order.  SubtreeVox
    % therefore identifies the contiguous descendant block without needing to
    % duplicate the ontology-building code from myAggregateHigherROIs.m.
    for r = 1:nBaseROI
        if ~ROI.IsAggregated(r)
            continue
        end
        target = double(ROI.SubtreeVox(r));
        cumulative = cumsum(double(ROI.Nvox(r:nBaseROI)));
        relEnd = find(cumulative >= target,1);
        if isempty(relEnd) || cumulative(relEnd) ~= target
            error('FigureS6:SubtreeMismatch', ...
                'Could not recover the descendant block for ROI %s.',string(ROI.Name(r)));
        end
        rr = r:(r+relEnd-1);
        w = double(ROI.Nvox(rr))';
        use = w > 0;
        rr = rr(use);
        w = w(use);

        vals = direct(:,rr);
        finite = isfinite(vals);
        numerator = sum(vals .* w .* finite,2,'omitnan');
        denominator = sum(finite .* w,2);
        animalROI(:,r) = numerator ./ denominator;
    end

    % myAggregateHigherROIs.m appends spelling aliases after labels.txt.
    % Copy the already-computed source row (same atlas label) for each alias.
    for r = nBaseROI+1:nROI
        src = find(double(ROI.label(1:nBaseROI)) == double(ROI.label(r)),1);
        if isempty(src)
            error('FigureS6:AliasMismatch', ...
                'No base ROI found for appended alias %s.',string(ROI.Name(r)));
        end
        animalROI(:,r) = animalROI(:,src);
    end
end

function slopes = fitMouseSlopes(x,validGroup,animalY)
% Fit x versus y independently for each mouse over the pooled fit's ROI set.
    nMouse = size(animalY,1);
    slopes = nan(nMouse,1);
    expectedNROI = sum(validGroup);
    for m = 1:nMouse
        use = validGroup & isfinite(animalY(m,:)');
        if sum(use) ~= expectedNROI
            error('FigureS6:AnimalROICount', ...
                'Mouse %d has %d/%d usable ROIs.',m,sum(use),expectedNROI);
        end
        coeff = polyfit(x(use),animalY(m,use)',1);
        slopes(m) = coeff(1);
    end
end

function out = countLabel(n)
% Compact actual-n label, also handles the unlikely case of tag-varying counts.
    u = unique(n(:));
    if isscalar(u)
        out = sprintf('%d',u);
    else
        out = mat2str(n(:)');
    end
end

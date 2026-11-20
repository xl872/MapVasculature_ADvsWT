%% myROIforCompare_selected_intensity.m
%
%  Signal intensity for the hand-drawn ROIs of findROI_intensity.m
%  (INTROI<name>.mat), computed the same way myROIforCompare.m computes the
%  atlas ROIs in its intensity section. This is the intensity counterpart of
%  myROIforCompare_selected.m, and it feeds the markers of
%  ../myComparePlotIntensity_marked.m.
%
%  As with the density masks, these are NOT whole anatomical structures:
%  findROI_intensity.m takes one coronal slice, keeps the voxels that passed the
%  significance threshold with the right sign, and saves that patch.
%
%  Space: the masks are stored in the un-permuted [345 495 255] grid, the one
%  imresize3(mean(WT00{i},4),sr) lands in, so there is no permute([3,2,1]) here.
%
%  Only tag = 1 exists for the intensity data, which is all that is needed.
%
%  Output: Intensity1_ROIsel_WTAD.mat -> table ROIsel
%          Intensity1_ROIsel_WTAD.csv

clc;clear all;close all;

%% ------------------------------- SETTINGS --------------------------------
vessDir   = '/path/to/vess';   % <-- set this: WORK_DIR with the intermediate results
sr      = 15;
tags    = 1;                       % the intensity data only exists for tag 1
dataPat = 'Density%dWTAD_seg_rawz.mat';    % holds WT00 / AD00 for the intensity
outPat  = 'Intensity%d_ROIsel_WTAD';

%  {name in the figure key , INTROI<...>.mat suffix , acronym for the x axis}
%
%  Checked against the atlas by looking up the labels each mask covers:
%    HP   -> CA1slm/so/sp/sr, CA3slu/so/sp/sr, DG-mo/po/sg   (8907 voxels) -> HIP
%    VISp -> VISp1, VISp2/3, VISp4, VISp5, VISp6a            (1488 voxels) -> VISp
SEL = { 'HP'   , 'HP'   , 'HIP'  ; ...
        'VISp' , 'VISp' , 'VISp' };

nSel     = size(SEL,1);
Name     = SEL(:,1);
MaskFile = SEL(:,2);
Acronym  = SEL(:,3);

%% ------------------------- 1. load the ROI masks --------------------------
idx  = cell(nSel,1);
Nvox = zeros(nSel,1);
for r = 1:nSel
    f = fullfile(vessDir,['INTROI',MaskFile{r},'.mat']);
    if ~exist(f,'file'), error('missing mask %s', f); end
    M       = load(f,'ROI');
    idx{r}  = find(M.ROI == 1);
    Nvox(r) = numel(idx{r});
    fprintf('%-6s <- INTROI%-8s %5d voxels\n', Name{r}, [MaskFile{r} '.mat'], Nvox(r));
end

%% ------------------------- 2. one table per tag ---------------------------
for tag = tags
    f = fullfile(vessDir, sprintf(dataPat, tag));
    if ~exist(f,'file'), warning('missing %s - skipped', f); continue; end
    load(f,'WT00','AD00');

    nRep = numel(WT00);
    Dw   = zeros(nRep,nSel);
    Da   = zeros(nRep,nSel);
    for c = 1:nRep
        tW = imresize3(mean(WT00{c},4),sr);
        tA = imresize3(mean(AD00{c},4),sr);
        for r = 1:nSel
            Dw(c,r) = mean(tW(idx{r}));
            Da(c,r) = mean(tA(idx{r}));
        end
    end

    WTmean = mean(Dw,1)';   WTstd = std(Dw,0,1)';
    ADmean = mean(Da,1)';   ADstd = std(Da,0,1)';

    ROIsel = table(Name, Acronym, MaskFile, Nvox, WTmean, WTstd, ADmean, ADstd);

    outBase = fullfile(vessDir, sprintf(outPat, tag));
    save(outBase,'ROIsel');
    writetable(ROIsel,[outBase '.csv']);

    fprintf('\ntag %d : %d ROIs -> %s.mat\n', tag, nSel, outBase);
    disp(ROIsel(:,{'Name','Acronym','Nvox','WTmean','ADmean'}));
end

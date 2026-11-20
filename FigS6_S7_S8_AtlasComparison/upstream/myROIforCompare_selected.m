%% myROIforCompare_selected.m
%
%  Vessel density for the hand-drawn ROIs that the bar plots use, i.e. the
%  ROI<name>.mat masks built by findROI.m, computed the same way
%  myROIforCompare.m computes the atlas ROIs. The point is to give
%  myComparePlotDensity_marked.m a y value for each marker that matches the
%  ROIbar.m box plots instead of the whole-structure atlas average.
%
%  Note these masks are NOT whole anatomical structures. findROI.m takes a single
%  coronal slice, keeps only the voxels that passed the significance threshold
%  (Ah) with the right sign of the effect, and saves that. So a mask is a small
%  patch inside the structure - which is exactly what the bar plots show.
%
%  Space: the masks are saved in the un-permuted [345 495 255] grid, the same one
%  imresize3(mean(WT00{i},4),sr) lands in, so unlike myROIforCompare.m there is
%  no permute([3,2,1]) here.
%
%  Output: Density<tag>_ROIsel_WTAD.mat -> table ROIsel
%          Density<tag>_ROIsel_WTAD.csv
%
%    WTmean/WTstd/ADmean/ADstd   mean and std over the 100 repeats in
%                                Density<tag>WTAD.mat  (same definition as the
%                                WTmean/ADmean columns of myROIforCompare.m)
%    WTmean_bar/WTstd_bar/...    mean and std over the 22 animals of
%                                forROIbarplot.mat - these are literally the
%                                numbers behind the ROIbar.m box plots, kept for
%                                cross-checking (identical for every tag)
%    Acronym                     acronym for looking the x value up in the paper
%                                spreadsheet

clc;clear all;close all;

%% ------------------------------- SETTINGS --------------------------------
vessDir   = '/path/to/vess';   % <-- set this: WORK_DIR with the intermediate results
sr      = 15;
tags    = [1 6 11 16 21 -1 -6 -11 -16 -21];

%  {name in the figure key , ROI<...>.mat suffix , acronym for the x axis}
%
%  The mask -> acronym column was checked by looking up which atlas labels each
%  mask actually covers:
%    Cg        -> MOs2/3, MOs5, ACAd2/3, ACAd5, ACAd6a, ACAv5   (about half MOs)
%    Septal    -> MS (494 vox) + LSr (3)
%    PIR       -> PIR1/2/3 (179) + OT1/2/3 (105)
%    SS        -> SSs1, SSs2/3, SSs4, SSs5
%    DG        -> DG-mo, DG-sg
%    CM        -> CM
%    LEn       -> ENTl1..ENTl5
%    MGN       -> MGd, MGv, MGm
%    MEn       -> mostly ENTl1/2a/2b (840) + ENTm1/ENTm2 (126)   <-- see note below
%    Amygdalar -> MEAav
%    IC        -> ICc, ICd, ICe
%    PAG5      -> PAG
%
%  MS, CM and MEAav are reported as "-" in the paper, so the x value comes from
%  the nearest structure that does have one: MSC, ILM and MEA.
SEL = { 'ACC'    , 'Cg'        , 'ACA'  ; ...
        'MS'     , 'Septal'    , 'MSC'  ; ...
        'Pir'    , 'PIR'       , 'PIR'  ; ...
        'S2'     , 'SS'        , 'SSs'  ; ...
        'DG'     , 'DG'        , 'DG'   ; ...
        'CM'     , 'CM'        , 'ILM'  ; ...
        'lEC'    , 'LEn'       , 'ENTl' ; ...
        'MGN'    , 'MGN'       , 'MG'   ; ...
        'mEC'    , 'MEn'       , 'ENTm' ; ...
        'Amy-HP' , 'Amygdalar' , 'MEA'  ; ...
        'IC'     , 'IC'        , 'IC'   ; ...
        'PAG'    , 'PAG5'      , 'PAG'  };

nSel     = size(SEL,1);
Name     = SEL(:,1);
MaskFile = SEL(:,2);
Acronym  = SEL(:,3);

%% ------------------------- 1. load the ROI masks --------------------------
idx  = cell(nSel,1);
Nvox = zeros(nSel,1);
for r = 1:nSel
    f = fullfile(vessDir,['ROI',MaskFile{r},'.mat']);
    if ~exist(f,'file'), error('missing mask %s', f); end
    M      = load(f,'ROI');
    idx{r} = find(M.ROI == 1);
    Nvox(r)= numel(idx{r});
    fprintf('%-8s <- ROI%-11s %5d voxels\n', Name{r}, [MaskFile{r} '.mat'], Nvox(r));
end

%% ---------- 2. the numbers the bar plots are actually drawn from ----------
%  Same computation as ROIbar.m: one value per animal, then averaged.
B    = load(fullfile(vessDir,'forROIbarplot.mat'));
nSub = size(B.WT.VDV,4);
barWT = zeros(nSub,nSel);
barAD = zeros(nSub,nSel);
for i = 1:nSub
    tW = imresize3(B.WT.VDV(:,:,:,i),sr);
    tA = imresize3(B.AD.VDV(:,:,:,i),sr);
    for r = 1:nSel
        barWT(i,r) = mean(tW(idx{r}));
        barAD(i,r) = mean(tA(idx{r}));
    end
end
WTmean_bar = mean(barWT,1)';   WTstd_bar = std(barWT,0,1)';
ADmean_bar = mean(barAD,1)';   ADstd_bar = std(barAD,0,1)';
fprintf('\nBar-plot reference values computed over %d animals\n\n', nSub);

%% ------------------------- 3. one table per tag ---------------------------
for tag = tags
    f = fullfile(vessDir,['Density',num2str(tag),'WTAD.mat']);
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

    ROIsel = table(Name, Acronym, MaskFile, Nvox, ...
                   WTmean, WTstd, ADmean, ADstd, ...
                   WTmean_bar, WTstd_bar, ADmean_bar, ADstd_bar);

    outBase = fullfile(vessDir,['Density',num2str(tag),'_ROIsel_WTAD']);
    save(outBase,'ROIsel');
    writetable(ROIsel,[outBase '.csv']);

    dWT = max(abs(WTmean - WTmean_bar));
    fprintf('tag %4d : %d ROIs, max|WTmean - barplot| = %.4f  -> %s.mat\n', ...
            tag, nSel, dWT, outBase);
end

fprintf(['\nNote: the mEC mask (ROIMEn.mat) is mostly lateral entorhinal by atlas\n' ...
         'label - 840/966 voxels are ENTl1/2a/2b and only 126 are ENTm. It overlaps\n' ...
         'the lEC mask (ROILEn.mat) heavily. Both are kept as you drew them, but the\n' ...
         'two markers are not independent regions.\n']);

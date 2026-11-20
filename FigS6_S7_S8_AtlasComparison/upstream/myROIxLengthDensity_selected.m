%% myROIxLengthDensity_selected.m
%
%  Estimates the reference (Ji et al. 2021) vessel length density and capillary
%  length density of the hand-drawn ROIs, instead of borrowing the value of one
%  whole structure.
%
%  Each mask covers a known set of atlas labels with known voxel counts, so the
%  spreadsheet value can be aggregated over them exactly the way the y values
%  were aggregated:
%
%        x_ROI = sum_i ( nvox_i * value_i ) / sum_i ( nvox_i )
%
%  where nvox_i is how many voxels of the mask carry label i, and value_i is the
%  spreadsheet value for that label. The paper reports "-" for many fine labels
%  (CA1slm, ENTl2a, MEAav, OT1, ...), so when a label has no value of its own the
%  nearest ancestor that does have one is used; the FracDirect column says how
%  much of the mask got a value without that fallback.
%
%  This matters most where a mask is not the structure its name suggests:
%  ROICg is about half MOs, ROIPIR is a third olfactory tubercle, and ROIMEn
%  is mostly lateral entorhinal.
%
%  Output: ROIsel_xvalue.mat -> table XSEL
%          ROIsel_xvalue.csv
%
%  Used by ../myComparePlotDensity_marked.m and
%  ../myComparePlotIntensity_marked.m when X_SOURCE = 'weighted'.

clc;clear all;close all;

%% ------------------------------- SETTINGS --------------------------------
vessDir   = '/path/to/vess';   % <-- set this: WORK_DIR with the intermediate results
xlsxFile  = fullfile(fileparts(mfilename('fullpath')),'..','data','1-s2.0-S0896627321000805-mmc2.xlsx');   % in ../data
labelFile = fullfile(vessDir,'labels.txt');
atlasBrik = fullfile(vessDir,'Allen_label_AMBMC30+orig.BRIK');
maskBrik  = fullfile(vessDir,'ambmcmask+orig.BRIK');
sr        = 15;
cropTail  = 165;

% {name in the figure key , mask file , acronym of the single spreadsheet row}
SEL = { 'ACC'    , 'ROICg.mat'        , 'ACA'  ; ...
        'MS'     , 'ROISeptal.mat'    , 'MSC'  ; ...
        'Pir'    , 'ROIPIR.mat'       , 'PIR'  ; ...
        'S2'     , 'ROISS.mat'        , 'SSs'  ; ...
        'DG'     , 'ROIDG.mat'        , 'DG'   ; ...
        'CM'     , 'ROICM.mat'        , 'ILM'  ; ...
        'lEC'    , 'ROILEn.mat'       , 'ENTl' ; ...
        'MGN'    , 'ROIMGN.mat'       , 'MG'   ; ...
        'mEC'    , 'ROIMEn.mat'       , 'ENTm' ; ...
        'Amy-HP' , 'ROIAmygdalar.mat' , 'MEA'  ; ...
        'IC'     , 'ROIIC.mat'        , 'IC'   ; ...
        'PAG'    , 'ROIPAG5.mat'      , 'PAG'  ; ...
        'HP'     , 'INTROIHP.mat'     , 'HIP'  ; ...
        'VISp'   , 'INTROIVISp.mat'   , 'VISp' };

nSel     = size(SEL,1);
Name     = SEL(:,1);
MaskFile = SEL(:,2);
Acronym  = SEL(:,3);

%% --------------------- 1. atlas on the mask grid --------------------------
fprintf('Building the atlas on the mask grid...\n');
mask  = BrikLoad(maskBrik);
mask  = imresize3(mask,1/sr);
mask  = imresize3(mask,sr);
mask  = permute(mask,[3,2,1]);
mskSz = size(mask);
clear mask

U = BrikLoad(atlasBrik);
U = permute(U,[3,2,1]);
U = imresize3(U,mskSz,'Method','nearest');
U = U(:,1:end-cropTail,:);
Uu = permute(U,[3,2,1]);        % -> [345 495 255], the grid the masks live in
clear U

%% ------------------------- 2. Allen hierarchy -----------------------------
T = allenTree(labelFile, xlsxFile);
if ~T.ok, warning('the hierarchy failed a consistency check - results may be off'); end

%% ------------- 3. spreadsheet value for every named structure -------------
CP = readtable(xlsxFile);
acrCP = strtrim(string(CP.Acronym));

% ThisStudy = vessel length density mean, Var17 = capillary length density mean
valSrc = { 'VLD', 'ThisStudy' ; ...
           'CLD', 'Var17'     };

nNamed  = numel(T.nmName);
resolved = struct();
for v = 1:size(valSrc,1)
    raw = str2double(string(CP.(valSrc{v,2})));

    own = nan(nNamed,1);
    for k = 1:nNamed
        j = find(acrCP == T.nmName(k), 1);
        if ~isempty(j), own(k) = raw(j); end
    end

    % walk up to the nearest ancestor that has a value
    val    = own;
    direct = ~isnan(own);
    for k = 1:nNamed
        if ~isnan(val(k)), continue; end
        p = T.parent(k);
        while p > 0
            if ~isnan(own(p)), val(k) = own(p); break; end
            p = T.parent(p);
        end
    end
    resolved.(valSrc{v,1}).val    = val;
    resolved.(valSrc{v,1}).direct = direct;
    fprintf('%s: %d/%d structures have their own value, %d filled from an ancestor, %d still empty\n', ...
            valSrc{v,1}, sum(direct), nNamed, sum(~direct & ~isnan(val)), sum(isnan(val)));
end

% label value -> index in the named list (0 = unnamed atlas slot)
lab2nm = zeros(max(T.labIdx)+1,1);
lab2nm(T.labIdx(T.nmRow)+1) = (1:nNamed)';

%% --------------------- 4. weighted value per mask -------------------------
VLD_est    = nan(nSel,1);  CLD_est    = nan(nSel,1);
VLD_single = nan(nSel,1);  CLD_single = nan(nSel,1);
FracDirVLD = nan(nSel,1);  FracDirCLD = nan(nSel,1);
Nvox       = zeros(nSel,1);
NLabels    = zeros(nSel,1);
Composition = strings(nSel,1);

rawVLD = str2double(string(CP.ThisStudy));
rawCLD = str2double(string(CP.Var17));

for r = 1:nSel
    f = fullfile(vessDir, MaskFile{r});
    if ~exist(f,'file'), warning('missing %s - skipped', f); continue; end
    M = load(f,'ROI');
    if ~isequal(size(M.ROI), size(Uu))
        warning('%s has size %s, expected %s - skipped', MaskFile{r}, ...
                mat2str(size(M.ROI)), mat2str(size(Uu)));
        continue;
    end

    lv = round(Uu(M.ROI == 1));
    Nvox(r) = numel(lv);
    u  = unique(lv);
    u  = u(u >= 0 & u <= max(T.labIdx));

    parts = strings(0);
    for v = {'VLD','CLD'}
        val = resolved.(v{1}).val;  dir = resolved.(v{1}).direct;
        wsum = 0; nsum = 0; ndir = 0;
        for a = u'
            k = lab2nm(a+1);
            if k == 0 || isnan(val(k)), continue; end     % unnamed slot / no value anywhere
            c    = sum(lv == a);
            wsum = wsum + c*val(k);
            nsum = nsum + c;
            if dir(k), ndir = ndir + c; end
        end
        if nsum > 0
            if strcmp(v{1},'VLD')
                VLD_est(r) = wsum/nsum;  FracDirVLD(r) = ndir/nsum;
            else
                CLD_est(r) = wsum/nsum;  FracDirCLD(r) = ndir/nsum;
            end
        end
    end

    for a = u'
        k = lab2nm(a+1);
        if k == 0, continue; end
        parts(end+1) = sprintf('%s(%d)', T.nmName(k), sum(lv == a)); %#ok<SAGROW>
    end
    NLabels(r)     = numel(parts);
    Composition(r) = strjoin(parts,' ');

    j = find(acrCP == string(Acronym{r}), 1);
    if ~isempty(j)
        VLD_single(r) = rawVLD(j);
        CLD_single(r) = rawCLD(j);
    end
end

%% ------------------------------ 5. save ----------------------------------
XSEL = table(Name, Acronym, MaskFile, Nvox, NLabels, ...
             VLD_est, VLD_single, FracDirVLD, ...
             CLD_est, CLD_single, FracDirCLD, Composition);

outBase = fullfile(vessDir,'ROIsel_xvalue');
save(outBase,'XSEL');
writetable(XSEL,[outBase '.csv']);

fprintf('\n%-8s %-8s %8s %8s %8s | %8s %8s %8s\n', ...
        'ROI','Acronym','VLDest','VLDsing','fracDir','CLDest','CLDsing','fracDir');
for r = 1:nSel
    fprintf('%-8s %-8s %8.4f %8.4f %8.2f | %8.4f %8.4f %8.2f\n', ...
        Name{r}, Acronym{r}, VLD_est(r), VLD_single(r), FracDirVLD(r), ...
        CLD_est(r), CLD_single(r), FracDirCLD(r));
end
fprintf('\nSaved %s.mat / .csv\n', outBase);

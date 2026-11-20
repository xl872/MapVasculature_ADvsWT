%% myAggregateHigherROIs.m
%
%  Fills in the higher-level ("parent") ROIs of the Allen hierarchy that have no
%  voxels carrying their own label in Allen_label_AMBMC30+orig, and therefore came
%  out as 0 in Density<tag>_ROI_WTAD.mat.
%
%  A parent is filled with the voxel-count-weighted average of every descendant
%  that does have a value:
%
%        mean_parent = sum_i ( nvox_i * mean_i ) / sum_i ( nvox_i )
%
%  Because every voxel of the atlas carries exactly one label, the descendants at
%  all depths of a sub-tree are disjoint, so they can all be summed in one go
%  without double counting (CA gets CA1, CA1slm, CA1so, CA1sp, CA1sr, CA2, ...).
%
%  Where the hierarchy comes from
%  ------------------------------
%  * 1-s2.0-S0896627321000805-mmc2.xlsx columns "Level 1" ... "Level 10" give the
%    full ancestor path of every structure the paper reports.
%  * labels.txt is a depth-first pre-order listing of the same Allen
%    ontology (verified: no illegal depth jumps, and every spreadsheet structure
%    falls inside the sub-tree block of the parent the spreadsheet declares for it).
%    That is what lets the ~235 sub-structures that are in labels.txt but absent
%    from the spreadsheet (CA1so, ENTl2a, PVHmm, MEApd-b, SUBd-sr, ...) still be
%    attached to the right parent, so none of their voxels are lost.
%
%  Where the weights come from
%  ---------------------------
%  Allen_label_AMBMC30+orig, put through exactly the same resample/crop that
%  myROIforCompare.m applied when the means were computed (imresize3 onto the
%  blocky mask grid, then U = U(:,1:end-165,:)), so the weights match the voxels
%  the means were actually taken over.
%
%  Output (written next to the input .mat files)
%  ---------------------------------------------
%        Density<tag>_ROI_WTAD_full.mat  -> variable ROI: the 6 columns of the input
%                                           plus Nvox / SubtreeVox / NvoxUsed /
%                                           NComp / IsAggregated
%        Density<tag>_ROI_WTAD_full.csv  -> same thing, for eyeballing
%
%  ../myComparePlotDensity_marked.m loads Density<tag>_ROI_WTAD_full and joins it
%  to the spreadsheet with
%        mergedTable = innerjoin(CP, ROI, 'LeftKeys','Acronym', 'RightKeys','Name');

clc; clear; close all;

%% ------------------------------- SETTINGS --------------------------------
vessDir   = '/path/to/vess';   % <-- set this: WORK_DIR with the intermediate results
xlsxFile  = fullfile(fileparts(mfilename('fullpath')),'..','data','1-s2.0-S0896627321000805-mmc2.xlsx');   % in ../data
labelFile = fullfile(vessDir,'labels.txt');
atlasBrik = fullfile(vessDir,'Allen_label_AMBMC30+orig.BRIK');
maskBrik  = fullfile(vessDir,'ambmcmask+orig.BRIK');

sr        = 15;    % same downsample factor as myROIforCompare.m
cropTail  = 165;   % same crop as myROIforCompare.m:  U = U(:,1:end-165,:)

% Tables to fill in. {input pattern , output pattern , tags}
% The hierarchy and the voxel weights are the same for all of them, so they are
% built once and reused.
JOBS = { 'Density%d_ROI_WTAD.mat'  , 'Density%d_ROI_WTAD_full'  , [1 6 11 16 21 -1 -6 -11 -16 -21] ; ...
         'Intensity%d_ROI_WTAD.mat', 'Intensity%d_ROI_WTAD_full', 1 };

% false -> keep every value that was measured directly, only fill the 0 rows.
%          Conservative; used for the tables in ../data.
% true  -> recompute EVERY parent as the weighted mean of its whole sub-tree,
%          including its own directly-labelled voxels. Worth trying: for a
%          structure like 'grey' or 'CTX' the direct value only describes the
%          leftover voxels that carry the generic label, not the whole structure.
%          Check the OwnVoxFrac printed below to see how much that matters.
RECOMPUTE_ALL = false;

% Cache for the voxel counts (delete this file to force a recount).
countFile = fullfile(vessDir,'ROI_voxelcount_AMBMC30.mat');

% Structures the spreadsheet and labels.txt spell differently. An extra row is
% appended under the spreadsheet spelling so innerjoin can find it.
% {spelling in the xlsx , spelling in labels.txt}
ALIAS = { 'VISl2/3' , 'VISl23'  ; ...
          'CUL4 5'  , 'CUL4,5'  };

% The six structures whose name gives no clue about their parent (they are not in
% the spreadsheet and are not named after their parent). Parents are the Allen
% CCF ones; each child does sit inside its parent's block in labels.txt.
% {child , parent}
PARENT_OVERRIDE = { 'isl' , 'OT'  ; ...
                    'islm', 'OT'  ; ...
                    'NC'  , 'PVZ' ; ...
                    'Mmme', 'MM'  ; ...
                    'ME'  , 'HY'  ; ...
                    'SCzo', 'SCs' };

%% -------------------- 1. voxel count per atlas label ----------------------
if exist(countFile,'file')
    fprintf('Loading cached voxel counts from %s\n', countFile);
    load(countFile,'voxCount');
else
    fprintf('Counting voxels in %s (one-off, takes a few minutes)...\n', atlasBrik);

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

    maxL     = double(max(U(:)));
    voxCount = histcounts(U, -0.5:1:(maxL+0.5))';   % voxCount(k) = #voxels labelled k-1
    clear U

    save(countFile,'voxCount','sr','cropTail','atlasBrik','maskBrik');
    fprintf('Saved voxel counts to %s\n', countFile);
end
nLabelCount = numel(voxCount);          % covers labels 0 .. nLabelCount-1

%% ------------------------- 2. read labels.txt -----------------------------
LAB     = readtable(labelFile);
labIdx  = double(LAB{:,1});
labName = strtrim(string(LAB{:,2}));
nLab    = numel(labIdx);

% "[123]" entries are unnamed atlas slots, not part of the ontology.
isNamed = cellfun(@isempty, regexp(cellstr(labName),'^\[\d+\]$','once')) & strlength(labName)>0;
nmRow   = find(isNamed);                % labels.txt rows that are real structures
nmName  = labName(nmRow);               % still in depth-first pre-order
nNamed  = numel(nmRow);
fprintf('labels.txt: %d entries, %d named structures\n', nLab, nNamed);

%% --------- 3. depth of each structure, from the spreadsheet paths ---------
CPh    = readtable(xlsxFile);           % same call myComparePlotDensity.m makes
nLevel = 10;

lvl = strings(height(CPh), nLevel);
for c = 1:nLevel
    col = CPh{:,c};
    s   = string(col);
    s(ismissing(s) | s == "NaN") = "";
    lvl(:,c) = strtrim(s);
end
acrXls = strtrim(string(CPh.Acronym));
acrXls(ismissing(acrXls)) = "";

depthMap = containers.Map('KeyType','char','ValueType','double');
for r = 1:height(CPh)
    p = lvl(r,:);
    p = p(strlength(p) > 0);
    % header / unit / footnote rows are skipped: a real row ends its path on its acronym
    if isempty(p) || acrXls(r) == "" || p(end) ~= acrXls(r), continue; end
    for k = 1:numel(p)
        key = char(p(k));
        if ~isKey(depthMap,key), depthMap(key) = k; end
    end
end
fprintf('Hierarchy: %d structures placed straight from the spreadsheet paths\n', depthMap.Count);

%% ------- 4. give every structure in labels.txt a depth in that tree -------
dep    = nan(nNamed,1);
howSet = strings(nNamed,1);
for k = 1:nNamed
    key = char(nmName(k));
    if isKey(depthMap,key)
        dep(k)    = depthMap(key);
        howSet(k) = "xlsx";
    end
end

% prefMat(k,j) = nmName(j) is a strict prefix of nmName(k).
% NB: startsWith(A,B) with two arrays tests "starts with ANY element of B", not
% element by element, so this has to be done by length group with strncmp.
lenName = strlength(nmName);
nmChar  = cellstr(nmName);
prefMat = false(nNamed,nNamed);
for k = 1:nNamed
    for L = 1:lenName(k)-1
        idx = find(lenName == L);
        if isempty(idx), continue; end
        hit = strncmp(nmChar{k}, nmChar(idx), L);
        prefMat(k, idx(hit)) = true;
    end
end

% Sub-structures the spreadsheet never lists (CA1so, PVHmm, ENTl2a, ...) are named
% after their parent, so hang them off the longest preceding name that is a prefix
% of theirs. Iterate to a fixed point (handles MEApd-a -> MEApd -> MEA).
pos     = (1:nNamed)';
changed = true;
while changed
    changed = false;
    for k = 1:nNamed
        if ~isnan(dep(k)), continue; end
        cand = find(pos < k & ~isnan(dep) & prefMat(k,:)');
        if ~isempty(cand)
            [~,b]     = max(lenName(cand));
            p         = cand(b);
            dep(k)    = dep(p) + 1;
            howSet(k) = "prefix";
            changed   = true;
        end
    end
end

% The handful whose name says nothing about their parent.
for a = 1:size(PARENT_OVERRIDE,1)
    k = find(nmName == string(PARENT_OVERRIDE{a,1}), 1);
    p = find(nmName == string(PARENT_OVERRIDE{a,2}), 1);
    if isempty(k) || isempty(p) || isnan(dep(p)), continue; end
    dep(k)    = dep(p) + 1;
    howSet(k) = "override";
end

% Safety net: anything still unplaced becomes a child of whatever structure
% immediately precedes it in the depth-first listing.
for k = 1:nNamed
    if ~isnan(dep(k)), continue; end
    p = find(pos < k & ~isnan(dep), 1, 'last');
    if isempty(p), dep(k) = 1; else, dep(k) = dep(p) + 1; end
    howSet(k) = "order";
end
fprintf('Hierarchy: %d from spreadsheet, %d by name prefix, %d by override, %d by listing order\n', ...
        sum(howSet=="xlsx"), sum(howSet=="prefix"), sum(howSet=="override"), sum(howSet=="order"));

% Sanity check: in a valid depth-first pre-order the depth never jumps by >1.
badJump = find(diff(dep) > 1);
if isempty(badJump)
    fprintf('Hierarchy check: depth-first pre-order is consistent (no illegal depth jumps)\n');
else
    warning('%d illegal depth jump(s), e.g. %s -> %s', ...
            numel(badJump), nmName(badJump(1)), nmName(badJump(1)+1));
end

%% ------------- 5. sub-tree of each structure (contiguous block) -----------
% In a depth-first pre-order the descendants of entry k are exactly the entries
% after k up to the first one whose depth is <= depth(k).
blockEnd = pos;
for k = 1:nNamed
    e = k;
    for j = k+1:nNamed
        if dep(j) <= dep(k), break; end
        e = j;
    end
    blockEnd(k) = e;
end

% End-to-end check: every structure the spreadsheet lists must land inside the
% sub-tree of the parent the spreadsheet declares for it.
nChk = 0; nBad = 0;
for r = 1:height(CPh)
    p = lvl(r,:);  p = p(strlength(p) > 0);
    if numel(p) < 2 || acrXls(r) == "" || p(end) ~= acrXls(r), continue; end
    kc = find(nmName == p(end), 1);
    kp = find(nmName == p(end-1), 1);
    if isempty(kc) || isempty(kp), continue; end     % e.g. 'root', not in labels.txt
    nChk = nChk + 1;
    if ~(kc > kp && kc <= blockEnd(kp))
        nBad = nBad + 1;
        if nBad <= 5, fprintf('  ! %s is not inside the sub-tree of %s\n', p(end), p(end-1)); end
    end
end
fprintf('Hierarchy check: %d/%d spreadsheet structures sit inside their declared parent\n', ...
        nChk-nBad, nChk);

%% ------------------------ 6. process each table --------------------------
for job = 1:size(JOBS,1)
inPattern = JOBS{job,1};  outPattern = JOBS{job,2};  tags = JOBS{job,3};
for tag = tags
    inFile = fullfile(vessDir, sprintf(inPattern, tag));
    if ~exist(inFile,'file')
        warning('missing %s - skipped', inFile);
        continue;
    end
    S    = load(inFile,'ROI');
    ROI  = S.ROI;
    nRoi = height(ROI);

    % ROI row for every labels.txt row (0 = not present)
    roiOfLab   = zeros(nLab,1);
    [tf,loc]   = ismember(labIdx, double(ROI.label));
    roiOfLab(tf) = loc(tf);
    nmRoi      = roiOfLab(nmRow);          % ROI row of each named structure

    % voxels carrying each label, indexed by ROI row
    nvox = zeros(nRoi,1);
    lb   = double(ROI.label);
    okL  = lb >= 0 & lb < nLabelCount;
    nvox(okL) = voxCount(lb(okL)+1);

    wt = double(ROI.WTmean);   ws = double(ROI.WTstd);
    ad = double(ROI.ADmean);   as = double(ROI.ADstd);
    wtNew = wt;  wsNew = ws;  adNew = ad;  asNew = as;

    subVox = zeros(nRoi,1);   % voxels in the whole sub-tree, incl. the node itself
    nUsed  = zeros(nRoi,1);   % voxels that actually went into the weighted mean
    nComp  = zeros(nRoi,1);   % how many sub-ROIs contributed
    isAgg  = false(nRoi,1);

    for k = 1:nNamed
        t = nmRoi(k);
        if t == 0, continue; end

        rr = nmRoi(k:blockEnd(k));         % the node itself + all its descendants
        rr = rr(rr > 0);
        if isempty(rr), continue; end

        n = nvox(rr);
        subVox(t) = sum(n);

        % --- WT ---
        v   = wt(rr);
        use = n > 0 & isfinite(v) & v ~= 0;
        if any(use) && (RECOMPUTE_ALL || wt(t) == 0 || ~isfinite(wt(t)))
            w        = n(use);
            wtNew(t) = sum(w .* v(use))        / sum(w);
            wsNew(t) = sum(w .* ws(rr(use)))   / sum(w);   % approximate, see note below
            nUsed(t) = sum(w);
            nComp(t) = sum(use);
            isAgg(t) = true;
        end

        % --- AD ---
        v   = ad(rr);
        use = n > 0 & isfinite(v) & v ~= 0;
        if any(use) && (RECOMPUTE_ALL || ad(t) == 0 || ~isfinite(ad(t)))
            w        = n(use);
            adNew(t) = sum(w .* v(use))        / sum(w);
            asNew(t) = sum(w .* as(rr(use)))   / sum(w);   % approximate, see note below
            if nUsed(t) == 0
                nUsed(t) = sum(w);
                nComp(t) = sum(use);
            end
            isAgg(t) = true;
        end
    end

    ROI.WTmean = wtNew;   ROI.WTstd = wsNew;
    ROI.ADmean = adNew;   ROI.ADstd = asNew;

    ownFrac = nan(nRoi,1);
    ownFrac(subVox > 0) = nvox(subVox > 0) ./ subVox(subVox > 0);

    ROI.Nvox         = nvox;      % voxels carrying this exact label
    ROI.SubtreeVox   = subVox;    % voxels in this structure incl. all descendants
    ROI.OwnVoxFrac   = ownFrac;   % Nvox / SubtreeVox
    ROI.NvoxUsed     = nUsed;     % voxels behind the aggregated value
    ROI.NComp        = nComp;     % number of sub-ROIs averaged
    ROI.IsAggregated = isAgg;

    % --- spelling aliases so the spreadsheet rows can join ---
    for a = 1:size(ALIAS,1)
        src = ALIAS{a,2};   dst = ALIAS{a,1};
        j = find(strcmp(ROI.Name, src), 1);
        if ~isempty(j) && ~any(strcmp(ROI.Name, dst))
            ROI = [ROI; ROI(j,:)];              %#ok<AGROW>
            ROI.Name{end} = dst;
        end
    end

    outBase = fullfile(vessDir, sprintf(outPattern, tag));
    save(outBase, 'ROI');
    writetable(ROI, [outBase '.csv']);

    nmRoiOk   = nmRoi(nmRoi > 0);
    nFilled   = sum(isAgg(nmRoiOk) & wt(nmRoiOk) == 0 & wtNew(nmRoiOk) ~= 0);
    nHadValue = sum(wt(nmRoiOk) ~= 0);
    nNowValue = sum(wtNew(nmRoiOk) ~= 0);
    fprintf('tag %4d : WT values %d -> %d  (+%d filled), still empty %d   -> %s.mat\n', ...
            tag, nHadValue, nNowValue, nFilled, sum(wtNew(nmRoiOk)==0), outBase);
end
end

%% ---------------------------- closing notes ------------------------------
fprintf('\nNotes\n-----\n');
fprintf(['* WTstd/ADstd on aggregated rows are the voxel-weighted average of the\n' ...
         '  component stds. That is only an approximation (the original stds are\n' ...
         '  across animals, not across voxels). The means are exact.\n']);
fprintf(['* Column OwnVoxFrac = Nvox./SubtreeVox tells you how much of a structure\n' ...
         '  its own directly-labelled voxels cover. Where that is small and the row\n' ...
         '  already had a value, the stored value describes only those leftover\n' ...
         '  voxels - set RECOMPUTE_ALL = true to replace it with the whole-structure\n' ...
         '  weighted mean instead.\n']);

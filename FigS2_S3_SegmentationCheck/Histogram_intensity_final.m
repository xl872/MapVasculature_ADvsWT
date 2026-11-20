%% =====================================================================
%  seg_agreement_metrics.m
%
%  Standalone. How much do the three segmentations agree, and how much would a
%  measured vessel density move if you had picked a different one.
%
%    Dice        2 |A n B| / (|A| + |B|)
%    Jaccard     |A n B| / |A u B|,  equal to D/(2-D), reported only by convention
%    Vd          vessel voxels divided by brain voxels
%    Vd CV       SD across the three masks divided by their mean, in percent
%
%  Everything runs on the WHOLE brain mask after the CROPEND crop, in 3D. Vd is
%  also computed per block with mysizedownd, which gives one CV per block and
%  therefore a distribution with a mean and an SD rather than a single number.
%
%  A caveat that belongs in the methods section. The three masks are
%  deterministic transforms of one acquisition, not independent repeats, so this
%  CV is sensitivity to the reconstruction choice, not measurement uncertainty.
% =====================================================================
clc; close all; clear;
% WORK_DIR = folder holding AD/, WT/ (registered volumes + segmentation .mat),
% the AMBMC/Allen BRIKs (HarvardDataverse/AMBMC) and the intermediate results.
WORK_DIR = '/path/to/vess';   % <-- set this
cd(WORK_DIR);                 % BrikLoad needs the HEAD/BRIK pair in the cwd

%% ================= parameters =================
TAG       = 16;
SUBJIDX   = 1:4;
CROPEND   = 165;          % voxels removed from the end of dim 2
ERODER    = 1;            % brain mask erosion radius in voxels, 0 disables
MASKFILE  = 'ambmcmask+orig.BRIK';

VARIANTS = struct( ...
  'key',   {'sr', ...
            'wsinc5', ...
            'nn'}, ...
  'label', {'SR segmentation', ...
            'wsinc5 segmentation', ...
            'NN segmentation'}, ...
  'glob',  {'./AD/*sr_reg1+orig_resultinrange50.mat', ...
            './AD/*regraw+orig_resultinrange50_raw.mat', ...
            './AD/*regNN+orig_resultinrange50_NN.mat'});

VOXMM     = 0.030;        % fallback voxel size in mm, used only if the header
                          % cannot be read
GRIDFILE  = '_sr_reg1+orig.HEAD';   % header consulted for the true voxel size

% Dice and Jaccard are reported for two scopes. The whole brain is the
% volumetric number. The slab is the same region the display figures show, which
% is useful for tying a number to a picture but is far smaller and therefore
% noisier, and is dominated by in plane agreement because it is only a few
% slices thick. Vd stays whole brain either way.
SLICE     = 175;          % slab centre, index along dim 3
STEP      = 1;            % slab half thickness in voxels
SLABYVIEW = [];           % [] keeps the full in plane extent. Set to 151:370 to
                          % match the display figures exactly.

VDBLOCK   = 15;           % block size in voxels for the per block Vd
VDSWEEP   = [5 10 15 20 30];   % block sizes to sweep, [] disables.
                          % CV depends strongly on block size, a smaller block
                          % holds fewer voxels and so carries more sampling
                          % noise. Publish the sweep, not one number alone.
VDMIN     = 0;            % ignore blocks whose mean Vd is at or below this

OUTDIR    = './figs';     % only the result .mat is written, no figures
if ~exist(OUTDIR,'dir'), mkdir(OUTDIR); end

nVar = numel(VARIANTS);


%% ================= brain mask =================
maskRaw = loadbrik(MASKFILE);

% Known offset between the brain mask file and the data, corrected the same way
% as in the display script so the two analyses use an identical brain mask.
m1 = zeros(size(maskRaw));
m1(1:339,:,1:243) = maskRaw(1:339,:,7:249);
maskRaw = m1 .* maskRaw;
maskRaw = maskRaw(:, 1:end-CROPEND, :);
if ERODER > 0
    maskRaw = imerode(maskRaw, strel('sphere', ERODER));
end
roi = maskRaw > 0;                       % whole brain, CROPEND already removed

roiSlab = roi;                           % the region the display figures show
roiSlab(:, :, [1:SLICE-STEP-1, SLICE+STEP+1:end]) = false;
if ~isempty(SLABYVIEW)
    keep = false(1, size(roi,2));
    keep(SLABYVIEW) = true;
    roiSlab(:, ~keep, :) = false;
end

SCOPES = {'brain', roi; 'slab', roiSlab};

fprintf('brain mask   %s\n', mat2str(size(roi)));
for sc = 1:size(SCOPES,1)
    fprintf('  %-6s %d voxels\n', SCOPES{sc,1}, nnz(SCOPES{sc,2}));
end
fprintf('\n');


%% ================= main loop =================
refList = dir(VARIANTS(1).glob);
if isempty(refList), error('no files matched %s', VARIANTS(1).glob); end

agree  = [];      % subj, scope (1 brain, 2 slab), varA, varB, Dice, Jaccard
tri    = [];      % subj, scope, three way overlap, n intersection, n union
vdGlob = [];      % subj, Vd per variant, mean, SD, CV percent
vdBlk  = [];      % subj, nBlocks, CV mean, CV SD, CV median, CV p25, CV p75
sweep  = [];      % subj, block voxels, block um, nBlocks, CV median, CV mean
cvPool = {};      % per block CV values, kept for the pooled summary

for index = SUBJIDX

    stem = refList(index).name(1:19);
    fprintf('=== %s ===\n', stem);
    root = refList(index).folder;

    % voxel size from the header when it is available
    vox = VOXMM;
    hf  = [root,'/',stem,GRIDFILE];
    if exist('BrikInfo','file') && exist(hf,'file')
        try, I = BrikInfo(hf); vox = abs(I.DELTA(1)); catch, end
    end

    % ---------- the three masks ----------
    Seg = cell(1, nVar);
    for iv = 1:nVar
        vl = dir(VARIANTS(iv).glob);
        kk = find(strncmp({vl.name}, stem, numel(stem)), 1);
        if isempty(kk)
            fprintf('  %s not found, skipped\n', VARIANTS(iv).key);
            continue;
        end
        S = load([vl(kk).folder,'/',vl(kk).name]);
        if ~isfield(S,'Vout'), error('%s has no Vout', vl(kk).name); end
        Seg{iv} = double(S.Vout{TAG}) > 0 & roi;
    end
    have = ~cellfun(@isempty, Seg);
    if nnz(have) < 2
        fprintf('  fewer than two masks available, skipping this subject\n');
        continue;
    end

    fprintf('  mask voxels ');
    for iv = 1:nVar
        if have(iv), fprintf('  %s %d', VARIANTS(iv).key, nnz(Seg{iv})); end
    end
    fprintf('\n');

    % ---------- Dice and Jaccard, both scopes side by side ----------
    fprintf('  %-22s %16s %16s\n', '', 'whole brain', 'slab');
    fprintf('  %-22s %8s %7s %8s %7s\n', 'pair', 'Dice', 'Jacc', 'Dice', 'Jacc');
    for a = 1:nVar
        for b = a+1:nVar
            if ~have(a) || ~have(b), continue; end
            v = nan(1,4);
            for sc = 1:size(SCOPES,1)
                r = SCOPES{sc,2};
                [dc, jc] = dice_jaccard(Seg{a} & r, Seg{b} & r);
                v(2*sc-1) = dc;  v(2*sc) = jc;
                agree(end+1,:) = [index sc a b dc jc];                %#ok<SAGROW>
            end
            fprintf('  %-22s %8.4f %7.4f %8.4f %7.4f\n', ...
                    sprintf('%s vs %s', VARIANTS(a).key, VARIANTS(b).key), v);
        end
    end
    % ---------- three way consensus overlap ----------
    % |A n B n C| / |A u B u C|, the generalised Jaccard. Reported alongside the
    % pairwise values because a claim about all three agreeing is not implied by
    % three pairwise numbers.
    if nnz(have) == nVar
        for sc = 1:size(SCOPES,1)
            r = SCOPES{sc,2};
            I3 = Seg{1} & Seg{2} & Seg{3} & r;
            U3 = (Seg{1} | Seg{2} | Seg{3}) & r;
            t3 = nnz(I3) / max(nnz(U3), 1);
            fprintf('  three way overlap, %-6s %8.4f   (%d of %d voxels)\n', ...
                    SCOPES{sc,1}, t3, nnz(I3), nnz(U3));
            tri(end+1,:) = [index sc t3 nnz(I3) nnz(U3)];              %#ok<SAGROW>
        end
    end

    fprintf('  slab mask voxels');
    for iv = 1:nVar
        if have(iv), fprintf('  %s %d', VARIANTS(iv).key, nnz(Seg{iv} & roiSlab)); end
    end
    fprintf('\n');

    % ---------- global Vd ----------
    vdG = nan(1, nVar);
    for iv = 1:nVar
        if have(iv), vdG(iv) = nnz(Seg{iv}) / nnz(roi); end
    end
    mG = mean(vdG,'omitnan');
    sG = std( vdG, 0, 'omitnan');            % sample SD, n-1
    fprintf('  global Vd  ');
    for iv = 1:nVar
        if have(iv), fprintf('  %s %.5f', VARIANTS(iv).key, vdG(iv)); end
    end
    fprintf('\n             mean %.5f   SD %.5f   CV %.2f%%   volume %.3f mm3\n', ...
            mG, sG, 100*sG/mG, mG*nnz(roi)*vox^3);
    vdGlob(end+1,:) = [index, vdG, mG, sG, 100*sG/mG];                %#ok<SAGROW>

    % ---------- per block Vd, and the CV distribution ----------
    [cv, nb] = block_cv(Seg, have, roi, VDBLOCK, VDMIN);
    st = cv_stats(cv);
    fprintf('  per block, %.0f um cube, n = %d\n', VDBLOCK*vox*1000, nb);
    fprintf('    CV   mean %.2f +/- %.2f %%   median %.2f   IQR %.2f to %.2f\n', ...
            st.mean, st.sd, st.med, st.p25, st.p75);
    vdBlk(end+1,:) = [index, nb, st.mean, st.sd, st.med, st.p25, st.p75]; %#ok<SAGROW>
    cvPool{end+1}  = cv(:);                                            %#ok<SAGROW>

    % The mean of a CV distribution is pulled by its right tail, so the median
    % is the more stable summary. Both are printed and both are stored.

    % ---------- block size sweep ----------
    if ~isempty(VDSWEEP)
        fprintf('  block size sweep\n');
        for bs = VDSWEEP
            [c2, n2] = block_cv(Seg, have, roi, bs, VDMIN);
            s2 = cv_stats(c2);
            fprintf('    %4.0f um   n %7d   CV median %6.2f%%   mean %6.2f%%\n', ...
                    bs*vox*1000, n2, s2.med, s2.mean);
            sweep(end+1,:) = [index, bs, bs*vox*1000, n2, s2.med, s2.mean]; %#ok<SAGROW>
        end
    end
    fprintf('\n');
end


%% ================= pooled =================
disp('================ pooled ================');

if ~isempty(agree)
    disp(' ');
    disp('Dice and Jaccard, mean +/- SD across subjects');
    for sc = 1:size(SCOPES,1)
        fprintf('  scope %s\n', SCOPES{sc,1});
        for a = 1:nVar
            for b = a+1:nVar
                k = agree(:,2) == sc & agree(:,3) == a & agree(:,4) == b;
                if ~any(k), continue; end
                fprintf('    %-20s Dice %.4f +/- %.4f   Jaccard %.4f +/- %.4f   (n = %d)\n', ...
                        sprintf('%s vs %s', VARIANTS(a).key, VARIANTS(b).key), ...
                        mean(agree(k,5)), std(agree(k,5)), ...
                        mean(agree(k,6)), std(agree(k,6)), nnz(k));
            end
        end
    end
end

if ~isempty(tri)
    fprintf('\n  three way overlap, mean +/- SD across subjects\n');
    for sc = 1:size(SCOPES,1)
        k = tri(:,2) == sc;
        if ~any(k), continue; end
        fprintf('    %-6s %.4f +/- %.4f   (n = %d)\n', ...
                SCOPES{sc,1}, mean(tri(k,3)), std(tri(k,3)), nnz(k));
    end
end

if ~isempty(vdGlob)
    fprintf('\n  global Vd CV across subjects   %.2f +/- %.2f %%   (n = %d)\n', ...
            mean(vdGlob(:,end)), std(vdGlob(:,end)), size(vdGlob,1));
end

if ~isempty(vdBlk)
    fprintf('  per block CV, subject medians  %.2f +/- %.2f %%   (n = %d)\n', ...
            mean(vdBlk(:,5)), std(vdBlk(:,5)), size(vdBlk,1));
    cvAll = cell2mat(cvPool(:));
    st = cv_stats(cvAll);
    fprintf('  per block CV, all blocks       mean %.2f +/- %.2f %%   median %.2f\n', ...
            st.mean, st.sd, st.med);
end

disp(' ');
disp('agree    subj  scope  varA  varB    Dice  Jaccard');
disp('         scope 1 brain, 2 slab.  variants 1 sr, 2 wsinc5, 3 nn');
disp(round(agree, 4));
disp('vdGlob   subj  Vd per variant, mean, SD, CV percent');
disp(round(vdGlob, 5));
disp('vdBlk    subj  nBlocks  CVmean  CVsd  CVmed  CVp25  CVp75');
disp(round(vdBlk, 3));

save(fullfile(OUTDIR,'seg_agreement_metrics.mat'), ...
     'agree','tri','vdGlob','vdBlk','sweep','cvPool','VARIANTS','TAG', ...
     'VDBLOCK','VDMIN','CROPEND','ERODER','SLICE','STEP','SLABYVIEW');


%% =====================================================================
%  local functions
% =====================================================================
function [V, I] = loadbrik(fn)
% Tolerates both BrikLoad output conventions, [err,V,Info] and [V,Info].
[a, b, c] = BrikLoad(fn);
if isstruct(b), V = a; I = b; else, V = b; I = c; end
V = double(V);
end


function [d, j] = dice_jaccard(A, B)
inter = nnz(A & B);
d = 2*inter / max(nnz(A) + nnz(B), 1);
j = inter    / max(nnz(A | B), 1);
end


function [cv, nb] = block_cv(Seg, have, roi, blk, vdmin)
% Vessel density per block for each mask, then the coefficient of variation
% across masks, block by block. All masks share the same brain denominator, so
% the blocks are comparable. Blocks that mysizedownd rejects for holding too
% little brain, and blocks with no vessel in any mask, return NaN rather than
% zero so they cannot drag the summary down.
maps = {};
for iv = 1:numel(Seg)
    if ~have(iv), continue; end
    maps{end+1} = mysizedownd(double(Seg{iv}), blk, double(roi));   %#ok<AGROW>
end
valid = mysizedownd(double(roi), blk, double(roi)) > 0;

V4 = cat(4, maps{:});
mu = mean(V4, 4);
sd = std( V4, 0, 4);            % sample SD, n-1
cv = 100 * sd ./ mu;
cv(~valid | mu <= vdmin) = NaN;
nb = nnz(~isnan(cv));
end


function s = cv_stats(cv)
v = cv(~isnan(cv));
if isempty(v)
    s = struct('mean',NaN,'sd',NaN,'med',NaN,'p25',NaN,'p75',NaN);
    return;
end
s = struct('mean', mean(v), 'sd', std(v), 'med', median(v), ...
           'p25', prctile(v,25), 'p75', prctile(v,75));
end


function output = mysizedownd(input, sr, mask0)
% Block mean of input over mask0, in sr by sr by sr blocks. Same as
% _common/mysizedownd.m, included so this script runs standalone.
input1 = zeros(ceil(size(input,1)/sr)*sr, ceil(size(input,2)/sr)*sr, ceil(size(input,3)/sr)*sr);
mask   = input1;
input1(1:size(input,1), 1:size(input,2), 1:size(input,3)) = input;
mask(  1:size(input,1), 1:size(input,2), 1:size(input,3)) = mask0;

Vdv = squeeze(sum(reshape(input1, sr, size(input1,1)/sr, size(input1,2), size(input1,3)), 1, 'omitnan'));
Vdv = squeeze(sum(reshape(Vdv, size(input1,1)/sr, sr, size(input1,2)/sr, size(input1,3)), 2, 'omitnan'));
Vdv = squeeze(sum(reshape(Vdv, size(input1,1)/sr, size(input1,2)/sr, sr, size(input1,3)/sr), 3, 'omitnan'));
mask = squeeze(sum(reshape(mask, sr, size(input1,1)/sr, size(input1,2), size(input1,3)), 1, 'omitnan'));
mask = squeeze(sum(reshape(mask, size(input1,1)/sr, sr, size(input1,2)/sr, size(input1,3)), 2, 'omitnan'));
mask = squeeze(sum(reshape(mask, size(input1,1)/sr, size(input1,2)/sr, sr, size(input1,3)/sr), 3, 'omitnan'));
output = Vdv ./ mask;
output(isnan(output)) = 0;
output(mask < sr^3*0.1) = 0;
end
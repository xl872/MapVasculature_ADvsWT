%% =================================================================
%  vessel_seg_check.m
%
%  Question this script answers
%    Vessels are DARK in the raw data. Does the segmentation mask spill onto
%    the BRIGHT rim around the vessels instead of sitting on the vessel lumen?
%    The brightest QBRIGHT of the mask voxels are flagged and drawn in a second
%    colour so the spill can be seen directly.
%
%  The 75 um data is *_75um_regonly+orig, on its own native 75 um grid.
%  Everything upstream of that file is preprocessing and is not touched here.
%  The 30 um view of it, used for the histograms and the line profiles, is
%  produced inside this script by a nearest neighbour lookup through the
%  headers, so the picture and the numbers can never come from different files
%  or different resampling settings. No separate regNN volume is read.
%
%  Figures 2, 3 and 3a draw TWO image layers in one axes, the 75 um underlay on
%  its own grid and the 30 um mask on its own grid, each placed by physical
%  position. Neither layer is resampled for display, so the underlay blocks are
%  exactly 75 um wide and the mask voxels are exactly where they belong. The
%  axes stay in 30 um index units so the detail box, the line coordinates and
%  the scale bar all keep their existing meaning.
%
%  The line profiles are not interpolated either. Every plotted value is a
%  value that was actually sampled. Profiles are aligned on the vessel centre
%  and the samples are binned onto each modality's own native step, 75 um for
%  the acquired data and 30 um for the SR.
%
%  Three segmentations are compared against the SAME underlay and the SAME z
%  score reference, so the only thing that differs between them is the mask
%    sr      segmented on the super resolution volume
%    wsinc5  segmented on the wsinc5 interpolated volume
%    NN      segmented on the nearest neighbour interpolated volume
%  The manual line profiles run once, on the reference variant, because they
%  compare the two IMAGES rather than the masks.
% =================================================================
clc; close all; clear;
% WORK_DIR = folder holding AD/, WT/ (registered volumes + segmentation .mat),
% the AMBMC/Allen BRIKs (HarvardDataverse/AMBMC) and the intermediate results.
WORK_DIR = '/path/to/vess';   % <-- set this
cd(WORK_DIR);                 % BrikLoad needs the HEAD/BRIK pair in the cwd
%% ================= parameters =================
QBRIGHT     = 0.05;         % flag this fraction of mask voxels, brightest first.
                            % The threshold comes from the mask's own
                            % distribution, so nothing outside the mask has to
                            % be classified as tissue for this to work.
TAG         = 16;
SUBJIDX     = 1:1;
FIGSUBJ     = 1;            % draw figures only for these; set to SUBJIDX for all
% The three segmentations. The first entry is the reference. It fixes the detail
% box for every variant so the three zoomed figures show the same region, and it
% is the only one that runs the manual line profiles.
VARIANTS = struct( ...
  'key',   {'sr', ...
            'wsinc5', ...
            'nn'}, ...
  'label', {'SR segmentation', ...
            'wsinc5 segmentation', ...
            'NN segmentation'}, ...
  'glob',  {'./AD/*sr_reg1+orig_resultinrange50.mat', ...
            './AD/*regraw+orig_resultinrange50_raw.mat', ...
            './AD/*regNN+orig_resultinrange50_NN.mat'}, ...
  'lines', {true, false, false});
SLICE       = 175;          % display slab centre, 30 um grid index
STEP        = 1;            % slab half thickness in 30 um voxels
CROPEND     = 165;          % voxels removed from the end of dim 2
YVIEW       = 151:370;      % columns drawn in the figures
ERODER      = 1;            % brain mask erosion radius in 30 um voxels, 0 disables
ZOOMSIZE    = [90 90];      % detail box size in voxels [rows cols]
ZOOMBOX     = [196 245 116 160];%[106 250 121 215]; %          % [c0 c1 r0 r1], empty means find it automatically
DISPLO      = 1;            % lower percentile for display contrast
DISPHI      = 99.5;         % upper percentile for display contrast
ROTCW       = true;         % rotate every image figure 90 degrees clockwise.
                            % Applied to the display arrays and to the 75 um
                            % coordinate vectors together, so both layers turn
                            % as one and the axes stay in 30 um index units.
                            % ZOOMBOX and LINES are therefore in the ROTATED
                            % frame; a set saved with ROTCW on is not valid with
                            % it off.
% Colours. Identity follows the entity, so the vessel mask is blue and the
% outlier group is red in every figure they appear in. The variant comparison
% needs no extra hues because all three series ARE masks, so they stay blue and
% the panel title carries the identity instead.
COLREST     = [0 0 1];               % 'b'  the vessel mask
COLFLAG     = [1 0 0];               % 'r'  the outliers, brightest QBRIGHT
COLRAW      = [0.051 0.580 0.533];   % #0D9488 teal,   raw 75 um profile
COLSR       = [0.482 0.173 0.749];   % #7B2CBF purple, SR 30 um profile
% All four were checked together on every pair. Worst case is purple against
% blue at dE 9.3 under protanopia, above the dE 8 floor, and blue against red is
% dE 35.7, so the two colours that matter most are far apart. The profile series
% also carry secondary encoding, square error bar markers on a coarse grid for
% the acquired data against round ones on a fine grid for the SR, so identity
% does not rest on hue alone.
ALPHAREST   = 0.50;
ALPHAFLAG   = 0.50;
OVERLAYMODE = 'fill';       % 'fill' or 'outline'; outline keeps the raw signal
                            % under the mask fully visible, which is often the
                            % better way to judge whether the mask spills
OUTLINEW    = 2;            % outline thickness in voxels
FONTPT        = 20;         % one font size for every piece of text in every figure
SCALEBAR_MAIN = 2;          % mm, figures 1, 2 and 5
SCALEBAR_ZOOM = 0.5;        % mm, figures 3 and 3a
SCALEBARPOS   = 'northwest';   % 'northwest' or 'southwest'
% Histogram. Bin edges are explicit so the same edges can be reused across
% subjects, and the counts are normalised by hand rather than through the
% Normalization option, which leaves no ambiguity about the denominator.
HISTNORM    = 'percentage'; % 'percentage' | 'probability' | 'pdf' | 'count'
HISTNBIN    = 100;
HISTRANGE   = [0.2 99.8];   % percentiles of the data used for the bin span
HISTYLOG    = false;        % log y makes the flagged tail readable
% z score reference. The slab brain mask supplies mean and std, computed once
% per subject so all three variants share it.
ZROBUST     = false;
% Manual line profiles, drawn on figure 3 of the reference variant. At that
% magnification a segment can be placed across a single vessel with voxel level
% precision. Along the covered voxels the 75 um data and the SR are compared.
% The voxels are taken exactly, with no interpolation, so the 75 um trace comes
% out at the acquisition step and that step IS the resolution limit.
DOLINES     = true;
NLINES      = 10;           % a SEM from three lines is close to meaningless
LINEMODE    = 'load';       % 'load' reuse the saved set, and draw only if there
                            %        is none or it was saved in a different
                            %        display frame
                            % 'draw' always draw a new set and overwrite the
                            %        saved one
LINEFILE    = '';           % where the set lives, empty means
                            % OUTDIR/lines_<subject>.mat. Endpoints are stored in
                            % PHYSICAL mm, not in pixels, so the same set stays
                            % valid after rotating the view, moving the detail
                            % box or changing the display window.
LINEFIG     = 'auto';       % when to write the annotated detail figure.
                            % 'auto' only after a fresh drawing session,
                            % 'always', or 'never'
LINES = [];         % explicit [x1 y1 x2 y2] per row in CURRENT display
                            % pixels. Overrides both modes and is never saved
                            % over. Frame dependent, unlike the saved file.
% Profile grids. Nothing is interpolated. Each covered voxel is assigned to the
% nearest bin centre on its own modality's grid and the samples in a bin are
% averaged, so every plotted value is a value that was actually sampled. The
% grids match the native steps, 75 um for the acquired data and 30 um for the
% SR, which is why the two traces carry different point densities.
PROFGRID75  = -0.3:0.075:0.3;
PROFGRID30  = -0.3:0.030:0.3;
PROFXLIM    = [];           % x limits of figure 6, always symmetric about the
                            % alignment origin so the two sides of the vessel can
                            % be compared by eye. Empty picks the widest
                            % symmetric window both grids actually cover. Set it
                            % explicitly, e.g. [-0.3 0.3], to keep the axis
                            % identical across subjects and runs.
ALIGNMODE   = 'mask';       % where zero sits on that axis.
                            % 'mask'   centroid of the mask voxels the line
                            %          crosses. Neutral, favours neither image.
                            % 'rawmin' or 'srmin' align on one trace's minimum,
                            %          which makes that trace look sharper and
                            %          deeper than it is. Diagnostic only.
ERRSTYLE    = 'both';       % 'bar' | 'band' | 'both'
MASKVOTE    = 0.5;          % a grid point counts as inside the mask when at
                            % least this fraction of lines has mask there. This
                            % sets the width of the blue band, so it is a
                            % choice, not a measurement. State it in the caption.
DO75        = true;         % run the native 75 um de-duplicated statistics
SLAB75MM    = [];           % slab thickness in mm for the 75 um stats,
                            % empty means match the display slab exactly
OCCMIN      = 0;            % keep 75 um voxels with occupancy strictly above this
OUTDIR      = './figs';
if ~exist(OUTDIR,'dir'), mkdir(OUTDIR); end
nVar = numel(VARIANTS);
%% ================= brain mask =================
maskRaw = loadbrik('ambmcmask+orig.BRIK');
% Correct the known offset between the brain mask file and the data. After this
% the mask header does not describe where the mask sits, so all physical
% coordinate work below uses the headers of the DATA, never the mask.
m1 = zeros(size(maskRaw));
m1(1:339,:,1:243) = maskRaw(1:339,:,7:249);
maskRaw = m1 .* maskRaw;
maskRaw = maskRaw(:, 1:end-CROPEND, :);
if ERODER > 0
    maskRaw = imerode(maskRaw, strel('sphere', ERODER));
end
maskFull = maskRaw > 0;                                    % whole brain, 3D
maskSlab = maskRaw;                                        % display slab only
maskSlab(:, [1:YVIEW(1)-1, YVIEW(end)+1:end], :) = 0;
maskSlab(:, :, [1:SLICE-STEP-1, SLICE+STEP+1:end])   = 0;
maskSlab = maskSlab > 0;
fprintf('brain mask voxels   whole brain %d   display slab %d\n\n', ...
        nnz(maskFull), nnz(maskSlab));
%% ================= main loop =================
refList = dir(VARIANTS(1).glob);
if isempty(refList), error('no files matched %s', VARIANTS(1).glob); end
Ins_all_c  = repmat({{}}, 1, nVar);   % z scored 75 um signal under each mask
Occ_flag_c = repmat({{}}, 1, nVar);
Occ_rest_c = repmat({{}}, 1, nVar);
summary30  = [];     % subj, variant, thrZ, nSeg, nFlag, median, mean, sd
summary75  = [];     % subj, variant, nInd, thr75, occFlag, occRest, p
lineData   = struct('p',{},'d',{},'raw',{},'sr',{},'inMask',{},'inFlag',{},'rho',{});
for index = SUBJIDX
    stem = refList(index).name(1:19);
    fprintf('================ %s ================\n', stem);
    root = refList(index).folder;
    % ---------- the 75 um data, on its own native grid ----------
    f75 = [root,'/',stem,'_75um_regonly+orig.BRIK'];
    if ~exist(f75,'file'), error('%s not found', f75); end
    [r75, I75] = loadbrik(f75);
    sz75 = size(r75);
    % ---------- the SR volume, 30 um ----------
    [srv, Isr] = loadbrik([root,'/',stem,'_sr_reg1+orig.BRIK']);
    srv  = srv(:, 1:end-CROPEND, :);   % crop is from the END of dim 2, so
    sz30 = size(srv);                  % ORIGIN and DELTA stay valid
    if any(Isr.ORIENT_SPECIFIC ~= I75.ORIENT_SPECIFIC)
        warning('orientation codes differ, array axes may not correspond');
    end
    % ---------- the 30 um view of the 75 um data ----------
    % Every 30 um voxel takes the exact value of the 75 um sample that contains
    % it. Derived here rather than read from a separate file, so the underlay
    % and the histogram can never disagree.
    raw = nn_from(r75, I75, Isr, sz30);
    [mu,  sd ] = ref_stats(raw(maskSlab), ZROBUST);
    [muS, sdS] = ref_stats(srv(maskSlab), ZROBUST);
    raw_z = (raw - mu ) ./ sd;
    sr_z  = (srv - muS) ./ sdS;
    fprintf('z score reference   75um mu %.5g sd %.5g   sr mu %.5g sd %.5g   n %d\n', ...
            mu, sd, muS, sdS, nnz(maskSlab));
    % ---------- prep for the de-duplicated 75 um statistics ----------
    have75 = false;
    if DO75
        vfrac  = prod(abs(Isr.DELTA)) / prod(abs(I75.DELTA));
        occBrn = min(reshape(accum_into(maskFull, Isr, I75, sz75), sz75)*vfrac, 1);
        kslab = (SLICE-STEP):(SLICE+STEP);
        zc    = Isr.ORIGIN(3) + (kslab-1)*Isr.DELTA(3);
        if isempty(SLAB75MM)
            zlo = min(zc) - abs(Isr.DELTA(3))/2;
            zhi = max(zc) + abs(Isr.DELTA(3))/2;
        else
            zlo = mean(zc) - SLAB75MM/2;
            zhi = mean(zc) + SLAB75MM/2;
        end
        zax = I75.ORIGIN(3) + (0:sz75(3)-1)*I75.DELTA(3);
        inZ = false(sz75);  inZ(:,:, zax >= zlo & zax <= zhi) = true;
        yc  = Isr.ORIGIN(2) + ([YVIEW(1) YVIEW(end)]-1)*Isr.DELTA(2);
        yax = I75.ORIGIN(2) + (0:sz75(2)-1)*I75.DELTA(2);
        inY = false(sz75);  inY(:, yax >= min(yc) & yax <= max(yc), :) = true;
        inBr = occBrn >= 0.5 & inZ & inY;
        z75v = (r75 - mean(r75(inBr))) ./ std(r75(inBr));
        have75 = true;
    end
    % ---------- display geometry ----------
    drawFigs = ismember(index, FIGSUBJ);
    if drawFigs
        % One rotation applied everywhere. rt turns a display array, and the
        % 75 um coordinate vectors are turned with it, so the two layers stay
        % registered to each other.
        if ROTCW, rt = @(A) rot90(A, -1); else, rt = @(A) A; end
        nRow0 = sz30(1);            % before the rotation
        nCol0 = numel(YVIEW);
        if ROTCW, nRow = nCol0; nCol = nRow0; else, nRow = nRow0; nCol = nCol0; end
        mmv = abs(Isr.DELTA(1));    % in plane isotropic, so rotating keeps it
        Asr = rt(squeeze(srv(:, YVIEW, SLICE)));
        loS = prctile(Asr(:), DISPLO);
        hiS = prctile(Asr(:), DISPHI);
        % The 75 um grid expressed in 30 um sub image index units. Both volumes
        % already live in the same base space, so this is a per axis scale and
        % shift with no rotation of its own, which is what image() can place.
        y75a = (I75.ORIGIN(1) + (0:sz75(1)-1)*I75.DELTA(1) - Isr.ORIGIN(1))/Isr.DELTA(1) + 1;
        x75a = (I75.ORIGIN(2) + (0:sz75(2)-1)*I75.DELTA(2) - Isr.ORIGIN(2))/Isr.DELTA(2) + 1 ...
               - (YVIEW(1) - 1);
        % The 75 um slice at the physical position of the display slice
        z30  = Isr.ORIGIN(3) + (SLICE-1)*Isr.DELTA(3);
        k75  = min(max(round((z30 - I75.ORIGIN(3))/I75.DELTA(3)) + 1, 1), sz75(3));
        A75a = squeeze(r75(:,:,k75));
        % Contrast from the part of that slice actually inside the window,
        % measured before the rotation so the window test is straightforward
        wr   = y75a >= 1 & y75a <= nRow0;
        wc   = x75a >= 1 & x75a <= nCol0;
        sub  = A75a(wr, wc);
        lo75 = prctile(sub(:), DISPLO);
        hi75 = prctile(sub(:), DISPHI);
        % Turn the underlay and its coordinates together. Clockwise means the
        % old top edge becomes the right edge, so new row = old column and
        % new column = nRow0 - old row + 1.
        A75 = rt(A75a);
        if ROTCW
            x75 = [nRow0 - y75a(end) + 1, nRow0 - y75a(1) + 1];
            y75 = [x75a(1), x75a(end)];
        else
            x75 = x75a;
            y75 = y75a;
        end
        fprintf('display   75 um slice %d of %d,  %d x %d blocks in the window,  rot90cw %d\n', ...
                k75, sz75(3), nnz(wr), nnz(wc), ROTCW);
        fw = 620;                                  % keep the window near the
        fh = round(fw * nRow / nCol) + 80;         % aspect of the rotated image
        % ---- figure 1, underlay only, no mask so it is drawn once ----
        f1 = figure('Name','raw only','Position',[60 520 fw fh],'Color','w');
        draw_grids(A75, x75, y75, lo75, hi75, [], [], nRow, nCol);
        add_scalebar(gca, mmv, SCALEBAR_MAIN, [nRow nCol], [0 0], FONTPT, SCALEBARPOS);
        title(sprintf('%s  75 um (regonly), slice %d', stem, SLICE), 'Interpreter','none');
        fontsize(f1, FONTPT, 'points');
        exportgraphics(f1, fullfile(OUTDIR, sprintf('%s_fig1_raw.png', stem)), ...
                       'Resolution', 300);
    end
    histAx = gobjects(1, nVar);
    edgesShared = [];    % fixed by the reference variant, reused by the others
    %% ============ inner loop over the three segmentations ============
    for iv = 1:nVar
        V  = VARIANTS(iv);
        vl = dir(V.glob);
        kk = find(strncmp({vl.name}, stem, numel(stem)), 1);
        if isempty(kk)
            fprintf('--- %s   no file matching %s, skipped\n', V.key, V.glob);
            continue;
        end
        S = load([vl(kk).folder,'/',vl(kk).name]);
        if ~isfield(S,'Vout')
            error('%s contains no variable named Vout', vl(kk).name);
        end
        Vout = S.Vout;
        fprintf('--- %s  (%s)\n', V.key, vl(kk).name);
        % ---------- segmentation ----------
        Seg3D = double(Vout{TAG}) > 0 & maskFull;   % full 3D, needed for occupancy
        Seg   = double(Vout{TAG}) > 0 & maskSlab;   % slab only, everything else
        if ~any(Seg(:))
            fprintf('    empty mask in the display slab, skipped\n');
            continue;
        end
        % ---------- brightness threshold ----------
        % From this mask's own distribution, so no voxel outside the mask has to
        % be labelled as tissue. thrZ is comparable across variants and subjects
        % because raw_z uses the same reference throughout.
        thrZ  = prctile(raw_z(Seg), 100*(1-QBRIGHT));
        flagB = Seg & raw_z >= thrZ;
        restB = Seg & ~flagB;
        fprintf('    bright threshold  z = %+.4f\n', thrZ);
        fprintf('    mask voxels in slab %d   flagged %d (%.1f%%)\n', ...
                nnz(Seg), nnz(flagB), 100*nnz(flagB)/max(nnz(Seg),1));
        fprintf('    median z   flagged %+.3f   rest %+.3f\n', ...
                median(raw_z(flagB)), median(raw_z(restB)));
        Ins_raw = raw_z(Seg);
        Ins_all_c{iv}{end+1} = Ins_raw(:);                                %#ok<SAGROW>
        summary30(end+1,:) = [index, iv, thrZ, nnz(Seg), nnz(flagB), ...
                              median(Ins_raw), mean(Ins_raw), std(Ins_raw)]; %#ok<SAGROW>
        % ---------- de-duplicated statistics on the native 75 um grid ----------
        % Occupancy uses the FULL 3D mask. Restricting to the slab first would
        % clip the vessel in z and bias the occupancy downwards.
        if have75
            occSeg = min(reshape(accum_into(Seg3D, Isr, I75, sz75), sz75)*vfrac, 1);
            selM   = inBr & occSeg > OCCMIN;
            nInd   = nnz(selM);
            if nInd < 50
                fprintf(['    [75 um] only %d independent voxels, the slab is ', ...
                         'thinner than one 75 um voxel. Set SLAB75MM to widen it.\n'], nInd);
            end
            if nInd > 0
                zm75  = z75v(selM);
                occm  = occSeg(selM);
                thr75 = prctile(zm75, 100*(1-QBRIGHT));
                isF   = zm75 >= thr75;
                occF  = occm(isF);
                occR  = occm(~isF);
                Occ_flag_c{iv}{end+1} = occF(:);                          %#ok<SAGROW>
                Occ_rest_c{iv}{end+1} = occR(:);                          %#ok<SAGROW>
                pv = NaN;
                if exist('ranksum','file') && numel(occF) > 3 && numel(occR) > 3
                    pv = ranksum(occF, occR);
                end
                fprintf('    [75 um] 30 um mask voxels %d -> independent 75 um %d (%.1f x)\n', ...
                        nnz(Seg3D), nInd, nnz(Seg3D)/max(nInd,1));
                fprintf('    [75 um] median occupancy   flagged %.3f   rest %.3f   p %.3g\n', ...
                        median(occF), median(occR), pv);
                summary75(end+1,:) = [index, iv, nInd, thr75, ...
                                      median(occF), median(occR), pv];    %#ok<SAGROW>
            end
        end
        if ~drawFigs, continue; end
        % ---------- mask layer for THIS variant, on the 30 um grid ----------
        Aflg = rt(squeeze(flagB(:, YVIEW, SLICE)));
        Arst = rt(squeeze(restB(:, YVIEW, SLICE)));
        Lrst = ov_layer(Arst, OVERLAYMODE, OUTLINEW);
        Lflg = ov_layer(Aflg, OVERLAYMODE, OUTLINEW);
        [rgbM, aM] = mask_layer(Lrst, Lflg, COLREST, COLFLAG, ALPHAREST, ALPHAFLAG);
        % The detail box is fixed by the reference variant and reused, otherwise
        % the three zoomed figures would show three different regions.
        if iv == 1
            if isempty(ZOOMBOX)
                zbox = auto_zoom_box(Aflg, ZOOMSIZE);
            else
                zbox = ZOOMBOX;
            end
            fprintf('    detail box [c0 c1 r0 r1] = %s\n', mat2str(zbox));
            fprintf('    to fix it by hand, zoom figure 3 then type round([xlim ylim])\n');
            % ---- figure 3a, plain underlay detail, same box, drawn once ----
            f3a = figure('Name','zoom plain','Position',[60 60 620 620],'Color','w');
            draw_grids(A75, x75, y75, lo75, hi75, [], [], nRow, nCol);
            xlim([zbox(1) zbox(2)]); ylim([zbox(3) zbox(4)]);
            add_scalebar(gca, mmv, SCALEBAR_ZOOM, [zbox(4)-zbox(3), zbox(2)-zbox(1)], ...
                         [zbox(1) zbox(3)], FONTPT, SCALEBARPOS);
            title(sprintf('%s  detail, 75 um only', stem), 'Interpreter','none');
            fontsize(f3a, FONTPT, 'points');
            exportgraphics(f3a, fullfile(OUTDIR, sprintf('%s_fig3a_zoom_plain.png', stem)), ...
                           'Resolution', 300);
        end
        pre = sprintf('%s_%s', stem, V.key);     % every file carries the variant
        %% ---- figure 2, overlay ----
        f2 = figure('Name',['overlay ' V.key],'Position',[700 480 fw fh+70],'Color','w');
        draw_grids(A75, x75, y75, lo75, hi75, rgbM, aM, nRow, nCol);
        %%
        rectangle('Position',[zbox(1) zbox(3) zbox(2)-zbox(1) zbox(4)-zbox(3)], ...
                  'EdgeColor',[1 1 0],'LineWidth',1.4);
        add_scalebar(gca, mmv, SCALEBAR_MAIN, [nRow nCol], [0 0], FONTPT, SCALEBARPOS);
        legend_outside(gca, {COLREST, COLFLAG}, ...
            {'vessel mask', sprintf('brightest %.0f%%  (z >= %+.2f)', 100*QBRIGHT, thrZ)});
        title(sprintf('%s   %s', stem, V.label), 'Interpreter','none');
        fontsize(f2, FONTPT, 'points');
        exportgraphics(f2, fullfile(OUTDIR, [pre '_fig2_overlay.png']), 'Resolution', 300);
        % ---- figure 3, detail with the overlay ----
        f3 = figure('Name',['zoom ' V.key],'Position',[700 60 620 620],'Color','w');
        draw_grids(A75, x75, y75, lo75, hi75, rgbM, aM, nRow, nCol);
        xlim([zbox(1) zbox(2)]); ylim([zbox(3) zbox(4)]);
        add_scalebar(gca, mmv, SCALEBAR_ZOOM, [zbox(4)-zbox(3), zbox(2)-zbox(1)], ...
                     [zbox(1) zbox(3)], FONTPT, SCALEBARPOS);
        title(sprintf('%s   %s, detail', stem, V.label), 'Interpreter','none');
        fontsize(f3, FONTPT, 'points');
        exportgraphics(f3, fullfile(OUTDIR, [pre '_fig3_zoom.png']), 'Resolution', 300);
        % ---- figure 4, distribution ----
        % The reference variant fixes the bin edges and every later panel reuses
        % them, so the three figures can be laid side by side and compared.
        f4 = figure('Name',['hist ' V.key],'Position',[1240 60 680 440],'Color','w');
        [edgesShared, fracOut] = plot_dist(Ins_raw, thrZ, QBRIGHT, HISTNORM, HISTNBIN, ...
                  HISTRANGE, HISTYLOG, COLREST, COLFLAG, FONTPT, edgesShared);
        title(sprintf('%s   %s', stem, V.label), 'Interpreter','none');
        exportgraphics(f4, fullfile(OUTDIR, [pre '_fig4_hist.png']), 'Resolution', 300);
        fprintf('    histogram   %d bins, z %+.2f to %+.2f, %.2f%% of voxels outside\n', ...
                HISTNBIN, edgesShared(1), edgesShared(end), 100*fracOut);
        fprintf('    distinct values %d  (each 75 um sample repeats up to %.0f times)\n', ...
                numel(unique(Ins_raw)), (abs(I75.DELTA(1))/abs(Isr.DELTA(1)))^3);
        % ---- figure 5, SR underlay. Both layers are 30 um here, so one
        %      composite image is exact and no two layer trick is needed.
        rgbSR = gray_rgb(Asr, loS, hiS);
        rgbSR = blend(rgbSR, Lrst, COLREST, ALPHAREST);
        rgbSR = blend(rgbSR, Lflg, COLFLAG, ALPHAFLAG);
        f5 = figure('Name',['sr underlay ' V.key],'Position',[1340 520 fw fh],'Color','w');
        imshow(rgbSR, 'InitialMagnification','fit'); hold on;
        add_scalebar(gca, mmv, SCALEBAR_MAIN, [nRow nCol], [0 0], FONTPT, SCALEBARPOS);
        title(sprintf('%s   %s on the SR underlay', stem, V.label), 'Interpreter','none');
        fontsize(f5, FONTPT, 'points');
        exportgraphics(f5, fullfile(OUTDIR, [pre '_fig5_sr.png']), 'Resolution', 300);
        %% ---- manual line profiles, reference variant only ----
        if DOLINES && V.lines
            AraZ = rt(squeeze(raw_z(:, YVIEW, SLICE)));
            AsrZ = rt(squeeze(sr_z( :, YVIEW, SLICE)));
            Amsk = rt(squeeze(Seg(  :, YVIEW, SLICE)));
            % Endpoints live in physical mm, so a saved set survives any change
            % to how the slice happens to be displayed. They are converted into
            % whatever the current frame is only at the moment they are used.
            lineFile = LINEFILE;
            if isempty(lineFile)
                lineFile = fullfile(OUTDIR, sprintf('lines_%s.mat', stem));
            end
            P = [];
            if ~isempty(LINES)
                P = LINES;
                fprintf('    lines taken from the LINES parameter\n');
            elseif strcmpi(LINEMODE, 'load')
                L = load_lines(lineFile);
                if isempty(L)
                    fprintf('    no saved set at %s, drawing a new one\n', lineFile);
                else
                    P = mm2disp(L.mm, Isr, nRow0, YVIEW(1), ROTCW);
                    fprintf('    %d lines loaded from %s\n', size(P,1), lineFile);
                    if abs(L.z - z30) > abs(Isr.DELTA(3))/2
                        fprintf(['    note, picked on the slice at z = %.3f mm, ', ...
                                 'now showing %.3f mm\n'], L.z, z30);
                    end
                    oob = any(P(:,[1 3]) < 1 | P(:,[1 3]) > nCol, 2) | ...
                          any(P(:,[2 4]) < 1 | P(:,[2 4]) > nRow, 2);
                    if any(oob)
                        fprintf('    %d fall outside the current window and are dropped\n', ...
                                nnz(oob));
                        P = P(~oob, :);
                    end
                end
            end
            drewNew = isempty(P);
            P = get_lines(f3, NLINES, P);
            fprintf('    LINES = %s\n', mat2str(round(P, 1)));
            if strcmpi(LINEFIG,'always') || (strcmpi(LINEFIG,'auto') && drewNew)
                fontsize(f3, FONTPT, 'points');   % the line labels came after
                exportgraphics(f3, fullfile(OUTDIR, [pre '_fig3_zoom_lines.png']), ...
                               'Resolution', 300);
            end
            % Every line contributes its covered voxels as they are. Zero is the
            % centroid of the mask voxels the line crosses, not the minimum of
            % either trace, because aligning on one image's minimum would make
            % that image look artificially sharp. Samples are then binned onto
            % each modality's own native step, with no interpolation, so every
            % value on the plot is a value that was actually sampled.
            g75 = PROFGRID75;
            g30 = PROFGRID30;
            GR = nan(size(P,1), numel(g75));   % acquired 75 um, coarse grid
            GS = nan(size(P,1), numel(g30));   % SR, fine grid
            GM = nan(size(P,1), numel(g30));   % mask, fine grid like the SR
            perLine = nan(size(P,1), 4);   % rho, min raw, min SR, offset in mm
            for iL = 1:size(P,1)
                [cc, rr, dmm] = line_voxels(P(iL,1:2), P(iL,3:4), mmv, [nRow nCol]);
                lin = sub2ind(size(AraZ), rr, cc);
                pR  = AraZ(lin);   pS = AsrZ(lin);
                inM = Amsk(lin) > 0;
                inF = Aflg(lin) > 0;
                [mR, iR] = min(pR);   [mS, iS] = min(pS);
                perLine(iL,:) = [corr(pR(:), pS(:)), mR, mS, dmm(iR)-dmm(iS)];
                switch lower(ALIGNMODE)
                    case 'mask',   d0 = mean(dmm(inM));
                    case 'rawmin', d0 = dmm(iR);
                    case 'srmin',  d0 = dmm(iS);
                end
                if isnan(d0), d0 = dmm(iR); end
                dc = dmm - d0;
                % Binning, not interpolation. Averaging within a line first keeps
                % the line as the unit of replication, so the SEM below is still
                % across lines and is not inflated by voxels that share a value.
                GR(iL,:) = bin_profile(dc, pR,          g75);
                GS(iL,:) = bin_profile(dc, pS,          g30);
                GM(iL,:) = bin_profile(dc, double(inM), g30);
                lineData(index, iL).p      = P(iL,:);   %#ok<SAGROW>
                lineData(index, iL).d      = dmm;
                lineData(index, iL).raw    = pR;
                lineData(index, iL).sr     = pS;
                lineData(index, iL).inMask = inM;
                lineData(index, iL).inFlag = inF;
                lineData(index, iL).rho    = perLine(iL,1);
            end
            % Geometry is the authority and is what took human effort. The
            % profiles are stored alongside as a record of this run, but they
            % are recomputed from the geometry every time, so a change to the
            % volumes or the z scoring can never leave a stale trace behind.
            if isempty(LINES)
                L = struct('mm',       disp2mm(P, Isr, nRow0, YVIEW(1), ROTCW), ...
                           'z',        z30, ...
                           'subject',  stem, ...
                           'variant',  V.key, ...
                           'ORIGIN',   Isr.ORIGIN, ...
                           'DELTA',    Isr.DELTA, ...
                           'profiles', lineData(index, 1:size(P,1)), ...
                           'note',     'mm endpoints are authoritative, profiles are a record');
                save(lineFile, 'L');
                fprintf('    %d lines saved to %s  (mm coordinates)\n', size(P,1), lineFile);
            end
            % ---- figure 6, averaged comparison, mean with SEM ----
            % Each modality is averaged on its own grid, so one point on the
            % teal trace is one acquired 75 um sample and one point on the
            % purple trace is one SR voxel.
            nR  = sum(~isnan(GR), 1);
            okR = nR >= max(2, ceil(0.5*size(P,1)));
            gkR = g75(okR);
            mRk = mean(GR,1,'omitnan');                       mRk = mRk(okR);
            eRk = std(GR,0,1,'omitnan')./sqrt(max(nR,1));     eRk = eRk(okR);
            % The SR trace and the mask band share the fine grid, so they share
            % one validity mask and cannot drift apart.
            nS  = sum(~isnan(GS), 1);
            okS = nS >= max(2, ceil(0.5*size(P,1)));
            gkS = g30(okS);
            mSk = mean(GS,1,'omitnan');                       mSk = mSk(okS);
            eSk = std(GS,0,1,'omitnan')./sqrt(max(nS,1));     eSk = eSk(okS);
            mMk = mean(GM,1,'omitnan');                       mMk = mMk(okS);
            if isempty(gkR) || isempty(gkS)
                warning('no grid point is reached by enough lines, figure 6 skipped');
            else
            wR = fwhm_dip(gkR, mRk);
            wS = fwhm_dip(gkS, mSk);
            % argmin of the coarse trace can only land on an acquired sample, so
            % the centroid below half depth is the comparable position estimate.
            cR = dip_centroid(gkR, mRk);
            cS = dip_centroid(gkS, mSk);
            dPos  = perLine(:,4);                     % per line argmin offset
            sePos = std(dPos,'omitnan')/sqrt(sum(~isnan(dPos)));
            f6 = figure('Name','average profile','Position',[1240 60 720 460],'Color','w');
            ax = gca; hold(ax,'on');
            yl = [min([mRk-eRk, mSk-eSk])-0.25, max([mRk+eRk, mSk+eSk])+0.55];
            hM = shade_runs(ax, gkS, mMk >= MASKVOTE, COLREST, 0.12, yl);
            if any(strcmpi(ERRSTYLE, {'band','both'}))
                sem_band(ax, gkR, mRk, eRk, COLRAW);
                sem_band(ax, gkS, mSk, eSk, COLSR);
            end
            hR = plot(ax, gkR, mRk, '-', 'Color',COLRAW, 'LineWidth',2.2);
            hS = plot(ax, gkS, mSk, '-', 'Color',COLSR,  'LineWidth',2.2);
            if any(strcmpi(ERRSTYLE, {'bar','both'}))
                % Both grids are short enough that every point carries its bar.
                errorbar(ax, gkR, mRk, eRk, 's', 'Color',COLRAW, ...
                         'LineStyle','none','LineWidth',1.4,'MarkerSize',6, ...
                         'MarkerFaceColor',COLRAW,'CapSize',4,'HandleVisibility','off');
                errorbar(ax, gkS, mSk, eSk, 'o', 'Color',COLSR, ...
                         'LineStyle','none','LineWidth',1.4,'MarkerSize',5, ...
                         'MarkerFaceColor',COLSR,'CapSize',4,'HandleVisibility','off');
            end
            % Symmetric about the alignment origin, and restricted to the window
            % both grids actually cover.
            if isempty(PROFXLIM)
                xm  = min([abs(gkR(1)) abs(gkR(end)) abs(gkS(1)) abs(gkS(end))]);
                xl2 = [-xm xm];
            else
                xl2 = PROFXLIM;
            end
            ylim(ax, yl); xlim(ax, xl2);
            xlabel(ax, 'Distance from the mask center (mm)');
            ylabel(ax, 'Zscored signal(mean \pm SEM)');
            % Explicit handles, so a mask band that breaks into several runs
            % cannot shift the labels onto the wrong series.
            legend(ax, [hM hR hS], {'Vessel mask','75 um','SR 30 um'}, ...
                   'Location','southwest', 'Box','off');
            title(ax, sprintf('FWHM  75um %.0f um   SR %.0f um      shift %+.0f um', ...
                              1000*wR, 1000*wS, 1000*(cR - cS)));
            fontsize(f6, FONTPT, 'points');
            exportgraphics(f6, fullfile(OUTDIR, [pre '_fig6_avgprofile.png']), ...
                           'Resolution', 300);
            dd = perLine(:,2) - perLine(:,3);
            fprintf('    averaged over %d lines\n', size(P,1));
            fprintf('    min depth   75um %+.3f   SR %+.3f   paired diff %+.3f +/- %.3f\n', ...
                    mean(perLine(:,2)), mean(perLine(:,3)), mean(dd), std(dd)/sqrt(numel(dd)));
            fprintf('    FWHM of the mean dip   75um %.0f um   SR %.0f um\n', 1000*wR, 1000*wS);
            fprintf('    dip centroid           75um %+.0f um   SR %+.0f um   shift %+.0f um\n', ...
                    1000*cR, 1000*cS, 1000*(cR - cS));
            fprintf('    per line argmin offset %+.0f +/- %.0f um   (n = %d)\n', ...
                    1000*mean(dPos,'omitnan'), 1000*sePos, sum(~isnan(dPos)));
            fprintf(['    NN quantisation alone is +/- %.0f um, so an offset smaller ', ...
                     'than that is not evidence of anything\n'], ...
                    1000*abs(I75.DELTA(1))/2);
            fprintf('    per line r   %s\n', mat2str(round(perLine(:,1)',3)));
            fprintf('    bins reached by enough lines   75um %d of %d   SR %d of %d\n', ...
                    numel(gkR), numel(g75), numel(gkS), numel(g30));
            if size(P,1) < 6
                fprintf('    SEM from %d lines is very noisy, raise NLINES to about 10\n', ...
                        size(P,1));
            end
            end
        end
    end   % variants
    %% ---- figure 7, the three histograms as small multiples ----
    % Same quantity on the same x axis, one panel per segmentation. No new hues
    % are needed because all three series are masks, so they stay blue and the
    % panel title carries the identity.
    if drawFigs
        f7 = figure('Name','variant comparison','Position',[100 40 720 900],'Color','w');
        t7 = tiledlayout(f7, nVar, 1, 'TileSpacing','compact', 'Padding','compact');
        for iv = 1:nVar
            if isempty(Ins_all_c{iv}), continue; end
            v = Ins_all_c{iv}{end};
            histAx(iv) = nexttile(t7);
            thrV = prctile(v, 100*(1-QBRIGHT));
            plot_dist(v, thrV, QBRIGHT, HISTNORM, HISTNBIN, HISTRANGE, HISTYLOG, ...
                      COLREST, COLFLAG, FONTPT, edgesShared);
            title(histAx(iv), sprintf('%s   n = %d,  flagged z >= %+.2f', ...
                  VARIANTS(iv).label, numel(v), thrV), 'Interpreter','none');
            if iv < nVar, xlabel(histAx(iv), ''); end
        end
        good = isgraphics(histAx);
        if nnz(good) > 1, linkaxes(histAx(good), 'x'); end
        title(t7, sprintf('%s   75 um signal under each mask, shared underlay', stem), ...
              'Interpreter','none');
        fontsize(f7, FONTPT, 'points');
        exportgraphics(f7, fullfile(OUTDIR, sprintf('%s_fig7_variants.png', stem)), ...
                       'Resolution', 300);
    end
end
%% ================= pooled results =================
% One set of edges for every pooled panel, taken from the reference variant.
edgesPool = [];
for iv = 1:nVar
    if isempty(Ins_all_c{iv}), continue; end
    va   = cell2mat(Ins_all_c{iv}(:));
    thrA = prctile(va, 100*(1-QBRIGHT));
    fg = figure('Name',['group hist ' VARIANTS(iv).key], ...
                'Position',[100 100 680 440],'Color','w');
    edgesPool = plot_dist(va, thrA, QBRIGHT, HISTNORM, HISTNBIN, HISTRANGE, HISTYLOG, ...
              COLREST, COLFLAG, FONTPT, edgesPool);
    title(sprintf('all subjects, %s, %d mask voxels', VARIANTS(iv).label, numel(va)), ...
          'Interpreter','none');
    exportgraphics(fg, fullfile(OUTDIR, ...
        sprintf('group_hist_%s.png', VARIANTS(iv).key)), 'Resolution', 300);
end
disp(' ');
disp('30 um grid   variant 1 = sr, 2 = wsinc5, 3 = nn');
disp('  subj   var     thrZ      nSeg    nFlag    median     mean       sd');
disp(round(summary30, 4));
if ~isempty(summary75)
    disp(' ');
    disp('native 75 um grid, independent voxels');
    disp('  subj   var    nInd     thr75    occFlag   occRest      p');
    disp(round(summary75, 4));
end
disp(' ');
disp('variant comparison, pooled over subjects');
fprintf('  %-22s %10s %10s %10s %10s\n', 'segmentation','nVoxels','median z','thr z','frac > 0');
for iv = 1:nVar
    if isempty(Ins_all_c{iv}), continue; end
    va = cell2mat(Ins_all_c{iv}(:));
    fprintf('  %-22s %10d %10.3f %10.3f %10.3f\n', VARIANTS(iv).label, numel(va), ...
            median(va), prctile(va, 100*(1-QBRIGHT)), mean(va > 0));
end
fprintf(['  A mask whose voxels sit on darker signal has a lower median z and a\n', ...
         '  lower fraction above zero. thr z says how far the brightest %.0f%% reaches.\n'], ...
        100*QBRIGHT);
save(fullfile(OUTDIR,'seg_check_result.mat'), ...
     'summary30','summary75','Ins_all_c','Occ_flag_c','Occ_rest_c','lineData', ...
     'VARIANTS','QBRIGHT','SLICE','STEP','TAG','ZROBUST','MASKVOTE', ...
     'PROFGRID75','PROFGRID30');
%% =================================================================
%  local functions
% =================================================================
function [V, I] = loadbrik(fn)
% Tolerates both BrikLoad output conventions, [err,V,Info] and [V,Info].
[a, b, c] = BrikLoad(fn);
if isstruct(b), V = a; I = b; else, V = b; I = c; end
V = double(V);
end
function V30 = nn_from(V75, I75, I30, sz30)
% Nearest neighbour view of a coarse volume on a fine grid, through the headers.
% Each output voxel receives the exact value of the coarse sample that contains
% it, so no value that was not acquired can appear. Indexing a 3D array with
% three index vectors does the whole resample in one operation.
sz75 = size(V75);
j1 = min(max(axmap(1:sz30(1), 1, I30, I75), 1), sz75(1));
j2 = min(max(axmap(1:sz30(2), 2, I30, I75), 1), sz75(2));
j3 = min(max(axmap(1:sz30(3), 3, I30, I75), 1), sz75(3));
V30 = V75(j1, j2, j3);
end
function j = axmap(i, k, Isrc, Idst)
% Index along axis k of one grid mapped to the nearest index of another grid.
j = round((Isrc.ORIGIN(k) + (i-1)*Isrc.DELTA(k) - Idst.ORIGIN(k)) ./ Idst.DELTA(k)) + 1;
end
function v = bin_profile(dc, vals, centres)
% Average of the samples that land in each bin. No interpolation, so every
% returned number is either a value that was sampled or the mean of a few that
% were. Bins that no sample reaches stay NaN rather than being filled in.
%
% Each sample is assigned to the nearest bin centre, which is why the bin width
% has to match the modality's own step. On the 75 um grid the samples inside one
% bin nearly all carry the same acquired value, so the average is that value.
v = nan(1, numel(centres));
if numel(centres) < 2 || isempty(dc), return; end
w  = centres(2) - centres(1);
k  = round((dc(:) - centres(1)) / w) + 1;
ok = k >= 1 & k <= numel(centres);
if ~any(ok), return; end
s = accumarray(k(ok), double(vals(ok)), [numel(centres) 1], @mean, NaN);
v = s(:).';
end
function [m, s] = ref_stats(v, robust)
% Centre and spread of the z score reference.
if robust
    m = median(v);
    s = 1.4826 * median(abs(v - m));
else
    m = mean(v);
    s = std(v);
end
end
function rgb = gray_rgb(A, lo, hi)
% Grayscale RGB. Percentile limits, never max, so one outlier voxel cannot
% crush the whole image to black.
g = (double(A) - lo) ./ max(hi - lo, eps);
rgb = repmat(min(max(g,0),1), 1, 1, 3);
end
function [rgbM, aM] = mask_layer(Lrst, Lflg, cRest, cFlag, aRest, aFlag)
% The mask as its own truecolour image plus a per pixel alpha, so it can share
% an axes with an underlay on a different grid without either being resampled.
sz   = size(Lrst);
rgbM = zeros([sz 3]);
aM   = zeros(sz);
for c = 1:3
    ch = zeros(sz);
    ch(Lrst) = cRest(c);
    ch(Lflg) = cFlag(c);
    rgbM(:,:,c) = ch;
end
aM(Lrst) = aRest;
aM(Lflg) = aFlag;
end
function ax = draw_grids(A75, x75, y75, lo75, hi75, rgbM, aM, nRow, nCol)
% Two image layers in one axes. The underlay sits on its native 75 um grid and
% the mask on the 30 um grid, each placed by physical position, so neither is
% resampled for display. The axes stay in 30 um index units, which is why the
% detail box, the line coordinates and the scale bar all still mean what they
% meant before. image() reads only the first and last of the coordinate vectors
% and spaces the rest linearly, which is exact for a regular grid.
ax = gca;
image(ax, [x75(1) x75(end)], [y75(1) y75(end)], gray_rgb(A75, lo75, hi75));
hold(ax, 'on');
if ~isempty(rgbM)
    hm = image(ax, [1 nCol], [1 nRow], rgbM);
    set(hm, 'AlphaData', aM);
end
set(ax, 'YDir', 'reverse');
daspect(ax, [1 1 1]);
xlim(ax, [1 nCol]); ylim(ax, [1 nRow]);
set(ax, 'XTick', [], 'YTick', [], 'Box', 'off');
end
function rgb = blend(rgb, m, col, a)
% Alpha blend one binary layer into an RGB image. Only used where both layers
% are already on the same grid.
w = a * double(m > 0);
for c = 1:3
    rgb(:,:,c) = (1-w).*rgb(:,:,c) + w.*col(c);
end
rgb = min(max(rgb,0),1);
end
function bx = auto_zoom_box(F, sz)
% Centre the detail box on the densest cluster of flagged voxels.
S = imgaussfilt(double(F), 6);
if ~any(S(:)), S = ones(size(F)); end
[~, k]  = max(S(:));
[r, c]  = ind2sub(size(F), k);
c0 = max(1, c-round(sz(2)/2)); c1 = min(size(F,2), c0+sz(2)); c0 = max(1, c1-sz(2));
r0 = max(1, r-round(sz(1)/2)); r1 = min(size(F,1), r0+sz(1)); r0 = max(1, r1-sz(1));
bx = [c0 c1 r0 r1];
end
function [edges, fracOut] = plot_dist(vSeg, thrZ, q, nrm, nbin, qrng, ylog, cSeg, cFlag, fpt, edgesIn)
% Distribution of the 75 um signal under the mask. One series, so no legend box,
% the title names it. The flagged tail is shaded.
%
% edgesIn lets several panels share one set of bin edges. That is not cosmetic.
% With 'percentage' the bar height is counts divided by N, so a panel whose bins
% happen to be wider gets taller bars for the same underlying shape. Panels meant
% to be compared must share edges.
%
% The denominator is always numel(vSeg), including values outside the drawn
% range, so the bars sum to slightly less than 100. fracOut reports how much.
if nargin < 11 || isempty(edgesIn)
    lo = prctile(vSeg, qrng(1));
    hi = prctile(vSeg, qrng(2));
    edges = linspace(lo, hi, nbin+1);
else
    edges = edgesIn;
end
lo = edges(1);  hi = edges(end);
fracOut = mean(vSeg < lo | vSeg > hi);
ax = gca; hold(ax,'on');
stairs(ax, edges, histprofile(vSeg, edges, nrm), '-', 'Color',cSeg, 'LineWidth',2.2);
if ylog, set(ax,'YScale','log'); end
yl = ylim(ax);
patch(ax, [thrZ hi hi thrZ], [yl(1) yl(1) yl(2) yl(2)], cFlag, ...
      'FaceAlpha',0.15, 'EdgeColor','none');
plot(ax, [thrZ thrZ], yl, '--', 'Color',cFlag, 'LineWidth',1.8);
xlim(ax, [lo hi]); ylim(ax, yl);
xlabel(ax, 'Zscored 75 um signal intensity');
ylabel(ax, norm_label(nrm));
box(ax,'off');
fontsize(gcf, fpt, 'points');
end
function y = histprofile(v, edges, nrm)
% Counts on the given edges, normalised, returned ready for stairs().
n = histcounts(v, edges);
switch lower(nrm)
    case 'percentage',  y = 100 * n / numel(v);
    case 'probability', y = n / numel(v);
    case 'pdf',         y = n / numel(v) ./ diff(edges);
    otherwise,          y = n;
end
y = [y y(end)];
end
function s = norm_label(nrm)
switch lower(nrm)
    case 'percentage',  s = 'Percentage of voxels (%)';
    case 'probability', s = 'Fraction of voxels';
    case 'pdf',         s = 'Probability density';
    otherwise,          s = 'Number of voxels';
end
end
function P = get_lines(fh, n, preset)
% Collect n line segments interactively from an existing figure, or replay a
% preset so a figure can be regenerated exactly. The clean copy of the
% figure has already been written to disk by the caller before this runs, so
% whatever gets drawn here only affects the annotated copy.
if ~isempty(preset)
    P = preset;
    figure(fh); ax = gca; hold(ax,'on');
    for k = 1:size(P,1)
        plot(ax, P(k,[1 3]), P(k,[2 4]), '-', 'Color',[1 1 1], 'LineWidth',1.5);
        text(ax, P(k,1), P(k,2), sprintf(' %d', k), ...
             'Color',[1 1 1], 'FontWeight','bold', 'FontSize',14);
    end
    return;
end
figure(fh);
ax  = gca;
old = get(get(ax,'Title'), 'String');
P   = zeros(n, 4);
hasROI = ~isempty(which('drawline'));
for k = 1:n
    title(ax, sprintf('line %d of %d   drag the ends, then double click to accept', k, n));
    if hasROI
        h = drawline(ax, 'Color',[1 1 1], 'LineWidth',1.5, 'Label',sprintf('%d',k));
        wait(h);
        pos = h.Position;
    else
        pos = ginput(2);
        plot(ax, pos(:,1), pos(:,2), '-w', 'LineWidth',1.5);
        text(ax, pos(1,1), pos(1,2), sprintf(' %d', k), ...
             'Color',[1 1 1], 'FontWeight','bold', 'FontSize',14);
    end
    P(k,:) = [pos(1,1) pos(1,2) pos(2,1) pos(2,2)];
end
title(ax, old);
end
function L = load_lines(fn)
% The saved set, or empty if there is none. No frame check is needed because the
% endpoints are stored in physical mm rather than in display pixels.
L = [];
if ~exist(fn, 'file'), return; end
S = load(fn);
if ~isfield(S, 'L') || ~isfield(S.L, 'mm') || isempty(S.L.mm), return; end
L = S.L;
end
function mm = disp2mm(P, I30, nRow0, ycol0, rot)
% Display pixel endpoints to physical mm along array dims 1 and 2. x is the
% display column and y is the display row. Clockwise rotation maps the display
% column onto dim 1 reversed and the display row onto dim 2.
mm = zeros(size(P,1), 4);
for k = 1:2
    x = P(:, 2*k-1);  y = P(:, 2*k);
    if rot
        i1 = nRow0 - x + 1;   i2 = ycol0 + y - 1;
    else
        i1 = y;               i2 = ycol0 + x - 1;
    end
    mm(:, 2*k-1) = I30.ORIGIN(1) + (i1-1)*I30.DELTA(1);
    mm(:, 2*k)   = I30.ORIGIN(2) + (i2-1)*I30.DELTA(2);
end
end
function P = mm2disp(mm, I30, nRow0, ycol0, rot)
% Inverse of disp2mm, into whatever the current display frame happens to be.
P = zeros(size(mm,1), 4);
for k = 1:2
    i1 = (mm(:, 2*k-1) - I30.ORIGIN(1)) ./ I30.DELTA(1) + 1;
    i2 = (mm(:, 2*k)   - I30.ORIGIN(2)) ./ I30.DELTA(2) + 1;
    if rot
        P(:, 2*k-1) = nRow0 - i1 + 1;   P(:, 2*k) = i2 - ycol0 + 1;
    else
        P(:, 2*k-1) = i2 - ycol0 + 1;   P(:, 2*k) = i1;
    end
end
end
function [c, r, d] = line_voxels(p1, p2, mmPerVox, imgSize)
% Every 30 um voxel the segment actually passes through, in order, with no
% interpolation. Oversample the segment, round to voxel indices, then drop
% repeats so each covered voxel appears exactly once.
len = hypot(p2(1)-p1(1), p2(2)-p1(2));
n   = max(2, ceil(4*len));
cc  = round(linspace(p1(1), p2(1), n));
rr  = round(linspace(p1(2), p2(2), n));
cc  = min(max(cc, 1), imgSize(2));
rr  = min(max(rr, 1), imgSize(1));
k   = [true, (diff(cc) ~= 0) | (diff(rr) ~= 0)];
c   = cc(k);
r   = rr(k);
d   = hypot(c - c(1), r - r(1)) * mmPerVox;
end
function h = shade_runs(ax, x, flag, col, a, yl)
% Shade the x ranges where flag is true, one patch per contiguous run. Only the
% first patch stays handle visible, so a band that breaks into several runs
% still contributes exactly one legend entry.
h = gobjects(1);
flag = logical(flag(:)');
d = diff([false flag false]);
s = find(d == 1);
e = find(d == -1) - 1;
for k = 1:numel(s)
    x0 = x(max(s(k)-1, 1));
    x1 = x(min(e(k)+1, numel(x)));
    p = patch(ax, [x0 x1 x1 x0], [yl(1) yl(1) yl(2) yl(2)], col, ...
              'FaceAlpha', a, 'EdgeColor','none');
    if k == 1, h = p; else, p.HandleVisibility = 'off'; end
end
if isempty(s)
    h = patch(ax, [NaN NaN NaN NaN], [NaN NaN NaN NaN], col, ...
              'FaceAlpha', a, 'EdgeColor','none');
end
end
function sem_band(ax, x, m, e, col)
% Faint mean plus or minus SEM band behind the mean line.
x = x(:)'; m = m(:)'; e = e(:)';
k = ~isnan(m) & ~isnan(e);
if nnz(k) < 2, return; end
x = x(k); m = m(k); e = e(k);
patch(ax, [x fliplr(x)], [m+e fliplr(m-e)], col, ...
      'FaceAlpha',0.16, 'EdgeColor','none', 'HandleVisibility','off');
end
function w = fwhm_dip(x, y)
% Full width at half depth of the central dip, measured against a baseline taken
% from the two ends of the profile. Returns NaN if the dip does not come back up
% to half depth inside the window, which is itself informative.
w = NaN;
n = numel(y);
if n < 5, return; end
e = max(1, round(0.12*n));
base = mean([y(1:e), y(end-e+1:end)], 'omitnan');
[ymin, k] = min(y);
if ~(ymin < base), return; end
half = base + 0.5*(ymin - base);
iL = find(y(1:k)   >= half, 1, 'last');
iR = k - 1 + find(y(k:end) >= half, 1, 'first');
if isempty(iL) || isempty(iR) || iL >= n || iR <= 1, return; end
w = lin_cross(x(iR-1), y(iR-1), x(iR), y(iR), half) ...
  - lin_cross(x(iL),   y(iL),   x(iL+1), y(iL+1), half);
end
function c = dip_centroid(x, y)
% Centre of mass of the part of the profile below half depth, weighted by how
% far below it goes. Unlike argmin this is not locked to the sampling grid, so
% it can be compared between a coarse trace and a fine one without the coarse
% grid contributing a spurious offset.
c = NaN;
n = numel(y);
if n < 5, return; end
e = max(1, round(0.12*n));
base = mean([y(1:e), y(end-e+1:end)], 'omitnan');
ymin = min(y);
if ~(ymin < base), return; end
half = base + 0.5*(ymin - base);
w = max(0, half - y);
w(isnan(w)) = 0;
if sum(w) == 0, return; end
c = sum(w(:).' .* x(:).') / sum(w);
end
function xc = lin_cross(x1, y1, x2, y2, yt)
if y2 == y1, xc = x1; else, xc = x1 + (yt - y1) * (x2 - x1) / (y2 - y1); end
end
function m = ov_layer(m, omode, w)
% Turn a binary layer into what actually gets drawn on the image.
if strcmpi(omode, 'outline')
    m = imdilate(bwperim(m > 0), strel('disk', max(w-1, 0)));
end
end
function add_scalebar(ax, mmPerVox, barMM, imgSize, offset, fpt, pos)
% White bar plus label. imgSize is [rows cols] of the visible region in 30 um
% voxels and offset shifts it when the axes are zoomed to a sub region. The axes
% run with YDir reverse, so a smaller y is higher on screen.
if nargin < 5 || isempty(offset), offset = [0 0]; end
if nargin < 6 || isempty(fpt),    fpt = 16;       end
if nargin < 7 || isempty(pos),    pos = 'northwest'; end
L  = barMM / mmPerVox;
x0 = offset(1) + imgSize(2)*0.06;
switch lower(pos)
    case 'northwest'
        y0 = offset(2) + imgSize(1)*0.09;
        ty = y0 + imgSize(1)*0.030;      % label sits below the bar
        va = 'top';
    otherwise                            % southwest
        y0 = offset(2) + imgSize(1)*0.93;
        ty = y0 - imgSize(1)*0.030;      % label sits above the bar
        va = 'bottom';
end
plot(ax, [x0 x0+L], [y0 y0], '-','Color',[1 1 1],'LineWidth',10);
text(ax, x0, ty, sprintf('%g mm', barMM), ...
     'Color',[1 1 1],'FontSize',fpt,'FontWeight','bold','VerticalAlignment',va);
end
function legend_outside(ax, cols, labs)
% Legend below the image rather than on top of it, so no data is hidden. The
% entries are proxy handles with no data, which is the only way to attach a real
% legend object to an image axes.
h = gobjects(1, numel(cols));
for k = 1:numel(cols)
    h(k) = plot(ax, NaN, NaN, 's', 'MarkerSize', 12, ...
                'MarkerFaceColor', cols{k}, 'MarkerEdgeColor','none');
end
lg = legend(ax, h, labs, 'Location','southoutside', 'Orientation','horizontal');
lg.Box = 'off';
lg.TextColor = [0.15 0.15 0.15];
end
function cnt = accum_into(maskSrc, Isrc, Idst, dstSize)
% Drop every nonzero source voxel into the destination voxel that contains it,
% using physical coordinates. Matrix sizes, FOV and origin may all differ. The
% count per destination voxel times the voxel volume ratio is the occupancy
% fraction, and the result is already de-duplicated.
idx = find(maskSrc);
if isempty(idx), cnt = zeros(prod(dstSize),1); return; end
[a1,a2,a3] = ind2sub(size(maskSrc), idx);
j1 = axmap(a1, 1, Isrc, Idst);
j2 = axmap(a2, 2, Isrc, Idst);
j3 = axmap(a3, 3, Isrc, Idst);
ok = j1>=1 & j1<=dstSize(1) & j2>=1 & j2<=dstSize(2) & j3>=1 & j3<=dstSize(3);
cnt = accumarray(sub2ind(dstSize, j1(ok), j2(ok), j3(ok)), 1, [prod(dstSize) 1]);
end
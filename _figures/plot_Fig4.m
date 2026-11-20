function plot_Fig4(outdir)
%PLOT_FIG4  Fig 4A: VD differential map at the top-50 % threshold on the five
%coronal slices of the 12 ROIs (drawing as plot_FD_slices.m: template + jet
%overlay, alpha 0.6 where significant: slices.Ah, BH-FDR q < 0.01 after
%smoothing, see README), drawn as AD - WT (red = higher VD in AD); vertical
%colour bar.
%Fig 4B: one box plot per ROI with the drawing code of ROIbar.m (boxplot,
%colours, linewidth 1, ttest2 stars, fontsize 10) plus the ROI name as title;
%only the figure size is set (see figpos). Values from Fig4B_ROI_values.mat.
%Writes Fig4A_ROI_targets.json (arrow targets for the ROI labels, inches from
%the top-left of panel A, with the number of significant voxels of each ROI
%mask in the map, n_sig = n_ADgtWT + n_ADltWT) and Fig4B_stats.csv.
%The target is the centroid of the ROI mask on its slice (left hemisphere for
%bilateral ROIs). For the ROIs in SIG_TARGET (DG, PAG: no significant voxel at
%the mask centroid) it is the centroid of the left-hemisphere AD > WT
%significant voxels of the mask, or the significant voxel nearest to it when
%that centroid is not on one.
S = figstyle();
here = fileparts(mfilename('fullpath'));
dd = fullfile(here,'..','Fig4_Fig5_FigS4_VesselDensityMaps','data');
D = load(fullfile(dd,'FD_VD_tag1.mat'),'slices','info');
sl = D.slices; crop = mapcrop(sl);
ys = [210 260 290 320 370];

%% A: map ---------------------------------------------------------------
gap = 0.02; tileW = (6.5-4*gap)/5;
[f, tiles] = mappanel(sl, D.info.climB, -1, D.info.climF, ys, 5, tileW, gap, S, crop);
exportpanel(f,outdir,'Fig4A_VD_map_AD-WT_50pct','raster'); close(f);

f = cbpanel(6.5/7,1.45,D.info.climF,[-0.1 0.1], ...
    {'Vessel Density (AD-WT),','p<0.01'},'vertical','normal');
exportpanel(f,outdir,'Fig4A_colorbar'); close(f);

%% ROI arrow targets (centroid of each ROI on its slice; left hemisphere for
%% bilateral ROIs) --------------------------------------------------------
rois   = {'Cg','Septal','PIR','SS','CM','DG','Amygdalar','LEn','MGN','MEn','IC','PAG5'};
labels = {'ACC','MS','Pir','S2','CM','DG','Amy-HP','lEC','MGN','mEC','IC','PAG'};
midline = {'Cg','Septal','CM','PAG5'};
cache = fullfile(dd,'Fig4A_ROI_centroids.mat');
if exist(cache,'file')
    C = load(cache); C = C.C;
else
    C = struct('roi',{},'y',{},'z',{},'x',{});
    for r = 1:numel(rois)
        M = load(fullfile(dd,'roi_masks',['ROI' rois{r} '.mat']));
        R = permute(M.ROI,[3 2 1]) > 0;            % [z y x] like FD
        y = find(squeeze(any(any(R,1),3)));
        [z,x] = find(squeeze(R(:,y(1),:)));
        if ~ismember(rois{r},midline), keep = x < 173; z = z(keep); x = x(keep); end
        C(end+1) = struct('roi',rois{r},'y',y(1),'z',mean(z),'x',mean(x)); %#ok<AGROW>
    end
    save(cache,'C');
end
nz = size(sl.U,1); u = S.mapUpsample;
SIG_TARGET = {'DG','PAG5'};
T = struct('label',{},'roi',{},'tile',{},'x_in',{},'y_in',{},'n_sig',{},'n_ADgtWT',{},'n_ADltWT',{},'target',{});
for r = 1:numel(rois)
    c = C(strcmp({C.roi},rois{r}));
    t = find(ys == c.y);
    % significant voxels of the ROI mask on its slice (AD - WT = -dF)
    M = load(fullfile(dd,'roi_masks',['ROI' rois{r} '.mat']));
    R = squeeze(permute(M.ROI,[3 2 1]) > 0); R = squeeze(R(:,c.y,:));
    j = find(sl.yIdx == c.y);
    sg = squeeze(sl.Ah(:,j,:)) & R; d = -squeeze(sl.dF(:,j,:));
    np = nnz(sg & d > 0); nn = nnz(sg & d < 0);
    tz = c.z; tx = c.x; how = 'mask centroid';
    if ismember(rois{r}, SIG_TARGET)
        Q = sg & d > 0; Q(:,173:end) = false;    % left hemisphere, AD > WT, significant
        [qz,qx] = find(Q);
        tz = mean(qz); tx = mean(qx); how = 'centroid of left AD>WT significant voxels';
        if ~Q(round(tz),round(tx))
            [~,k] = min((qz-tz).^2 + (qx-tx).^2);
            fprintf('Fig 4A %-7s centroid of the %d left AD>WT voxels (z %.1f, x %.1f) is not on one: nearest is (z %d, x %d)\n', ...
                labels{r}, numel(qz), tz, tx, qz(k), qx(k));
            tz = qz(k); tx = qx(k); how = 'left AD>WT significant voxel nearest to their centroid';
        end
        [qz0,qx0] = find(sg);
        fprintf('Fig 4A %-7s target (z %.1f, x %.1f) -> (z %.1f, x %.1f); nearest significant voxel to the old target: %.1f voxels\n', ...
            labels{r}, c.z, c.x, tz, tx, sqrt(min((qz0-c.z).^2 + (qx0-c.x).^2)));
    end
    row = nz - tz + 1 - crop(1) + 1;             % display row in the cropped tile
    col = tx - crop(3) + 1;
    nr = (crop(2)-crop(1)+1); nc = (crop(4)-crop(3)+1);
    T(end+1) = struct('label',labels{r},'roi',rois{r},'tile',t, ...
        'x_in',tiles(t).left + (col-0.5)/nc*tiles(t).width, ...
        'y_in',tiles(t).top  + (row-0.5)/nr*tiles(t).height, ...
        'n_sig',np+nn,'n_ADgtWT',np,'n_ADltWT',nn,'target',how); %#ok<AGROW>
    fprintf('Fig 4A %-7s y=%d: %d significant voxels (AD>WT %d, AD<WT %d) of %d\n', labels{r}, c.y, np+nn, np, nn, nnz(R));
end
fid = fopen(fullfile(outdir,'Fig4A_ROI_targets.json'),'w');
fprintf(fid,'%s',jsonencode(struct('tiles',tiles,'targets',T,'crop',crop,'upsample',u),'PrettyPrint',true)); fclose(fid);

%% B: ROI box plots (ROIbar.m), in the order of the figure ----------------
V = load(fullfile(dd,'Fig4B_ROI_values.mat'));
order = {'Cg','Septal','PIR','SS','DG','CM','LEn','MGN','MEn','Amygdalar','IC','PAG5'};
labs  = {'ACC','MS','Pir','S2','DG','CM','lEC','MGN','mEC','Amy-HP','IC','PAG'};
Tb = table();
for r = 1:numel(order)
    i = find(string(V.names) == order{r});
    MdenWT = V.MdenWT; MdenAD = V.MdenAD;
    f = figure(Position=figpos(6.5/7,1.45,10));
    h = boxplot([MdenWT(:,i),MdenAD(:,i)],'Labels',{'WT','AD'}, ...
        'colors', [0, 0.4470, 0.7410;0.8500, 0.3250, 0.0980]);
    set(h,'linewidth',1)
    xt = get(gca, 'XTick');
    yl=ylim;
    yl=yl(1)+0.95*(yl(2)-yl(1));
    [~,p]=ttest2(MdenWT(:,i),MdenAD(:,i));
    hold on
    st = '';
    if p<0.001
        st = '***';
    elseif p<0.01
        st = '**';
    elseif p<0.05
        st = '*';
    end
    if ~isempty(st), text(mean(xt([1 2])), yl, st,HorizontalAlignment='center'), end
    ylabel('Vessel density')
    title(labs{r},'FontWeight','normal')
    fontsize(10,"points")
    exportpanel(f,outdir,sprintf('Fig4B_%02d_%s',r,labs{r})); close(f);
    Tb = [Tb; table(string(order{r}),string(labs{r}),sum(~isnan(MdenWT(:,i))),sum(~isnan(MdenAD(:,i))), ...
        mean(MdenWT(:,i),'omitnan'),mean(MdenAD(:,i),'omitnan'),p,string(st), ...
        'VariableNames',{'roi','label','nWT','nAD','meanWT','meanAD','p','stars'})]; %#ok<AGROW>
end
disp(Tb);
writetable(Tb, fullfile(outdir,'Fig4B_stats.csv'));
end

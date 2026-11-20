function plot_FD_slices(fdFile, outFile)
%PLOT_FD_SLICES  Redraw a final differential (FD) map figure (Fig 4A, 5A-E, 6A)
%from the small data/FD_*.mat.
%
%   plot_FD_slices('data/FD_VD_tag1.mat')                 % Fig 4A
%   plot_FD_slices('data/FD_intensity_tag1.mat','x.jpg')  % Fig 6A
%
% The drawing code is the plot loop of loadandplot.m /
% Random_group_result_double_*.m (18 coronal slices, every 10th voxel,
% template underlay, jet overlay with alpha 0.6 where significant & in brain).
% FD.slices holds exactly the arrays that loop reads; FD.FD is the full 3D
% map (dF .* Ah .* mask) for anyone who wants other slices.
S = load(fdFile,'slices','info');
sl = S.slices; climB = S.info.climB; climF = S.info.climF;
figure('Position',[100,100,1800,700]);
tiledlayout(3,6,'TileSpacing','none','Padding','tight');
for j = 1:size(sl.U,2)
    nexttile;
    B = squeeze(sl.U(:,j,:));
    B = repmat(mat2gray(double(B),double(climB)),[1,1,3]);
    imagesc(B); axis xy; axis image off; hold on;
    F = squeeze(sl.dF(:,j,:));
    hF = imagesc(F,climF); axis xy;
    colormap('jet');
    set(hF,'AlphaData',0.6*double(squeeze(sl.Ah(:,j,:))).*double(squeeze(sl.mask(:,j,:))));
end
if nargin > 1 && ~isempty(outFile)
    print(gcf,outFile,'-djpeg','-r300');
end
end

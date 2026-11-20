function plot_Fig5(outdir)
%PLOT_FIG5  Fig 5A-E: VD differential maps for the threshold ranges 2-10,
%10-20, 20-30, 30-40 and 40-50 % (tag -21, -16, -11, -6, -1), AD - WT, on the
%six coronal slices y = 210 ... 360 (drawing as plot_FD_slices.m); horizontal
%jet colour bar. Tile placement is written to Fig5_tiles.json (for the region
%arrows).
S = figstyle();
here = fileparts(mfilename('fullpath'));
dd = fullfile(here,'..','Fig4_Fig5_FigS4_VesselDensityMaps','data');
tags  = [-21 -16 -11 -6 -1];
names = {'2-10','10-20','20-30','30-40','40-50'};
ys = [210 230 270 300 330 360];
gap = 0.02; panelW = (6.5-2*0.16)/3; tileW = (panelW-gap)/2;
J = struct();
for k = 1:5
    D = load(fullfile(dd,sprintf('FD_VD_tag%d.mat',tags(k))),'slices','info');
    if k == 1, crop = mapcrop(D.slices); end
    [f, tiles] = mappanel(D.slices, D.info.climB, -1, D.info.climF, ys, 2, tileW, gap, S, crop);
    exportpanel(f,outdir,sprintf('Fig5%s_VD_map_AD-WT_%spct',char('A'+k-1),names{k}),'raster'); close(f);
    J.(char('A'+k-1)) = tiles;
    climF = D.info.climF;
end
J.crop = crop; J.upsample = S.mapUpsample; J.frame = [size(D.slices.U,1) size(D.slices.U,3)];
fid = fopen(fullfile(outdir,'Fig5_tiles.json'),'w'); fprintf(fid,'%s',jsonencode(J,'PrettyPrint',true)); fclose(fid);
f = cbpanel(1.9,0.45,climF,[-0.1 0.1],'Vessel Density (AD-WT), p<0.01','horizontal','normal');
exportpanel(f,outdir,'Fig5_colorbar'); close(f);
end

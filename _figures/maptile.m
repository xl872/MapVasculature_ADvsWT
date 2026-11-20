function [RGB, info] = maptile(sl, j, climB, signF, climF, S, crop)
%MAPTILE  One coronal slice of a differential map as an RGB image: template
%underlay (grey) + jet overlay with alpha 0.6 where significant and in brain,
%exactly as plot_FD_slices.m / loadandplot.m, then cropped to the brain and
%upsampled S.mapUpsample x for display.
%   sl     : slices struct of data/FD_*.mat (U, dF, Ah, mask : [z 18 x])
%   j      : slice index (1..18)
%   signF  : +1 keeps dF (WT - AD), -1 shows AD - WT
%   crop   : [r0 r1 c0 c1] in display rows/cols (rows top->bottom)
B = mat2gray(double(squeeze(sl.U(:,j,:))), double(climB));
F = signF*double(squeeze(sl.dF(:,j,:)));
A = S.mapAlpha*double(squeeze(sl.Ah(:,j,:))).*double(squeeze(sl.mask(:,j,:)));
% axis xy in plot_FD_slices.m -> flip rows for image order
B = flipud(B); F = flipud(F); A = flipud(A);
B = B(crop(1):crop(2),crop(3):crop(4));
F = F(crop(1):crop(2),crop(3):crop(4));
A = A(crop(1):crop(2),crop(3):crop(4));
u = S.mapUpsample;
if u > 1
    B = min(max(imresize(B,u,'bicubic'),0),1);
    F = imresize(F,u,'bilinear');
    A = min(max(imresize(A,u,'bilinear'),0),1);
end
cm = jet(256);
k = round((F-climF(1))/diff(climF)*255)+1; k = min(max(k,1),256);
C = reshape(cm(k,:),[size(F) 3]);
RGB = (1-A).*repmat(B,[1 1 3]) + A.*C;
info.size = size(B);
end

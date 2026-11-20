function crop = mapcrop(sl, margin)
%MAPCROP  Common crop box [r0 r1 c0 c1] (display rows top->bottom, cols) that
%contains the brain mask of all 18 slices, plus a margin in voxels. Using one
%box for every map keeps all tiles (Fig 4A, 5, 6A) at the same scale.
if nargin < 2, margin = 6; end
M = squeeze(any(sl.mask > 0.05, 2));   % [z x]
M = flipud(M);
rows = find(any(M,2)); cols = find(any(M,1));
crop = [max(rows(1)-margin,1), min(rows(end)+margin,size(M,1)), ...
        max(cols(1)-margin,1), min(cols(end)+margin,size(M,2))];
end

function [f, tiles] = mappanel(sl, climB, signF, climF, yList, nCols, tileW, gap, S, crop)
%MAPPANEL  Grid of map tiles (row-major over yList) in a new figure whose size
%is set by the tile width; tiles are separated by thin white gaps.
%   tiles : struct array with the placement of every tile in inches from the
%           top-left of the panel (for adding arrows/labels afterwards)
nT = numel(yList); nRows = ceil(nT/nCols);
asp = (crop(2)-crop(1)+1)/(crop(4)-crop(3)+1);
tileH = tileW*asp;
W = nCols*tileW + (nCols-1)*gap; H = nRows*tileH + (nRows-1)*gap;
f = newpanel(W,H);
tiles = struct('y',{},'col',{},'row',{},'left',{},'top',{},'width',{},'height',{});
for t = 1:nT
    j = find(sl.yIdx == yList(t), 1);
    if isempty(j), error('slice y=%d not in the saved slices', yList(t)); end
    RGB = maptile(sl, j, climB, signF, climF, S, crop);
    r = ceil(t/nCols); c = t - (r-1)*nCols;
    left = (c-1)*(tileW+gap); top = (r-1)*(tileH+gap);
    ax = axes(f,'Units','inches','Position',[left H-top-tileH tileW tileH]);
    image(ax,RGB);
    set(ax,'XLim',[0.5 size(RGB,2)+0.5],'YLim',[0.5 size(RGB,1)+0.5],'YDir','reverse','Visible','off');
    tiles(end+1) = struct('y',yList(t),'col',c,'row',r,'left',left,'top',top, ...
        'width',tileW,'height',tileH); %#ok<AGROW>
end
end

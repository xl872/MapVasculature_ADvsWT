function layoutScatterKey(ax, lgFit, axKey, lgKey, varargin)
%LAYOUTSCATTERKEY  Lay out a scatter axes with a fit legend underneath and a
%   marker key to the right, using the legends' and the axes' own measured
%   sizes instead of hard-coded positions. Survives changes of figure size and
%   font size.
%
%   layoutScatterKey(ax, lgFit, axKey, lgKey)
%   layoutScatterKey(..., 'Name', value)
%
%   ax     main axes
%   lgFit  legend with the fit equations (create it with Location 'southoutside')
%   axKey  invisible axes holding the dummy markers of the key, [] if there is none
%   lgKey  legend of that key, [] if there is none
%
%   Pass axKey = [] and lgKey = [] when the figure has no marker key: the axes
%   then simply take the full width.
%
%   Call it last, after every label, title and font size is set - it measures
%   them. Room for the x/y labels and the tick labels comes from ax.TightInset,
%   so the margins below are only the outer padding of the figure.
%
%   Options (normalized figure units):
%     PadL 0.015  PadR 0.015  PadT 0.025  PadB 0.015   outer padding
%     GapKey 0.02   gap between the axes and the key legend
%     GapFit 0.02   gap between the x label and the fit legend
%     MinKeyFontSize 6
%
%   The key is made to fit vertically by first trying one column, then two,
%   then shrinking its font down to MinKeyFontSize.

p = inputParser;
p.addParameter('PadL', 0.015);
p.addParameter('PadR', 0.015);
p.addParameter('PadT', 0.025);
p.addParameter('PadB', 0.015);
p.addParameter('GapKey', 0.02);
p.addParameter('GapFit', 0.02);
p.addParameter('MinKeyFontSize', 6);
p.parse(varargin{:});
o = p.Results;

hasKey = ~isempty(lgKey) && isgraphics(lgKey);

ax.Units = 'normalized';
lgFit.Units = 'normalized';
if hasKey
    axKey.Units = 'normalized';
    lgKey.Units = 'normalized';
end

% --- make the key fit the available height -------------------------------
if hasKey
    availH = 1 - o.PadT - o.PadB;
    lgKey.NumColumns = 1;
    drawnow;
    if lgKey.Position(4) > availH
        lgKey.NumColumns = 2;
        drawnow;
    end
    while lgKey.Position(4) > availH && lgKey.FontSize > o.MinKeyFontSize
        lgKey.FontSize = max(o.MinKeyFontSize, lgKey.FontSize - 0.5);
        drawnow;
    end
    kw = lgKey.Position(3);
    kh = lgKey.Position(4);
else
    kw = -o.GapKey;   % cancels the gap, so the axes take the full width
    kh = 0;
end

% --- size of the fit legend, measured while it is still auto-placed ------
drawnow;
fw = lgFit.Position(3);
fh = lgFit.Position(4);

% --- the axes box --------------------------------------------------------
% TightInset is [left bottom right top]: what the labels need outside the box.
% It shifts as the box moves, so settle it with a couple of passes.
axPos = ax.Position;
for pass = 1:3
    ti  = ax.TightInset;
    axL = o.PadL + ti(1);
    axB = o.PadB + fh + o.GapFit + ti(2);
    axW = 1 - axL - kw - o.GapKey - o.PadR - ti(3);
    axH = 1 - axB - o.PadT - ti(4);
    if axW <= 0.05 || axH <= 0.05
        warning('layoutScatterKey:tooSmall', ...
                ['the labels and legends leave almost no room for the axes ' ...
                 '(w %.2f, h %.2f) - use a bigger figure or a smaller font'], axW, axH);
        axW = max(axW, 0.05);
        axH = max(axH, 0.05);
    end
    axPos = [axL axB axW axH];
    ax.Position = axPos;
    drawnow;
end

% Setting Position detaches a legend from its Location, so place the fit legend
% explicitly - otherwise it would keep re-shrinking the axes.
lgFit.Position = [axPos(1) + (axPos(3) - fw)/2, o.PadB, fw, fh];
if hasKey
    axKey.Position = axPos;
    lgKey.Position = [axPos(1) + axPos(3) + o.GapKey, ...
                      axPos(2) + (axPos(4) - kh)/2, kw, kh];
end
end

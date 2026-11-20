function f = cbpanel(w_in, h_in, cl, ticks, label, orient, direction)
%CBPANEL  MATLAB jet colorbar as its own panel (placed w_in x h_in inches,
%10-pt text printed at 7 pt, see figpos).
%   orient    'vertical'   tick labels and (rotated) label on the right
%             'horizontal' bold label above the bar and the end values on
%                          either side of it (Fig 5/6)
%   direction 'normal' or 'reverse'
f = figure(Position=figpos(w_in,h_in,10));
ax = axes(f,'Visible','off');
colormap(ax,jet); clim(ax,cl);
cb = colorbar(ax);
cb.Ticks = ticks;
cb.Direction = direction;
fontsize(f,10,'points')
if strcmp(orient,'vertical')
    cb.Position = [0.14 0.10 0.13 0.80];
    cb.AxisLocation = 'in';                  % tick labels and label on the right
    cb.Label.String = label;
    cb.Label.FontWeight = 'bold';
else
    cb.Location = 'southoutside';
    cb.Position = [0.17 0.12 0.66 0.30];
    cb.TickLabels = repmat({''},1,numel(ticks));
    ends = {sprintf('%g',cl(1)), sprintf('%g',cl(2))};
    if strcmp(direction,'reverse'), ends = ends([2 1]); end
    tb = {'EdgeColor','none','VerticalAlignment','middle','FontSize',10,'Margin',0,'Interpreter','none'};
    annotation(f,'textbox',[0 0.12 0.155 0.30],'String',ends{1},'HorizontalAlignment','right',tb{:});
    annotation(f,'textbox',[0.845 0.12 0.155 0.30],'String',ends{2},'HorizontalAlignment','left',tb{:});
    annotation(f,'textbox',[0 0.50 1 0.48],'String',label,'EdgeColor','none', ...
        'HorizontalAlignment','center','VerticalAlignment','bottom', ...
        'FontSize',10,'FontWeight','bold','Margin',0,'Interpreter','none');
end
end

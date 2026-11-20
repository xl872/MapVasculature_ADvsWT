function exportpanel(f, outdir, name, varargin)
%EXPORTPANEL  Save one panel as vector (PDF, SVG, EMF) and 600-dpi PNG at its
%exact figure size (no tight cropping, so panels line up when assembled).
%   exportpanel(f, outdir, name)              % vector + png
%   exportpanel(f, outdir, name, 'raster')   % image-only panel: png + pdf/svg
%                                             % (pdf/svg embed the image)
%Uses exportgraphics(...,'Padding','figure') where available (R2025a+), else
%print with the paper set to the figure size (R2023b: the panels are made
%with R2023b, whose figures boxplotGroup supports).
if ~exist(outdir,'dir'), mkdir(outdir); end
base = fullfile(outdir, name);
drawnow;
isRaster = any(strcmpi(varargin,'raster'));
if ~isMATLABReleaseOlderThan('R2025a')
    if isRaster, ct = 'image'; else, ct = 'vector'; end
    eg = {'Padding','figure','BackgroundColor','white'};
    exportgraphics(f,[base '.pdf'],'ContentType',ct,'Resolution',600,eg{:});
    exportgraphics(f,[base '.svg'],'ContentType',ct,'Resolution',600,eg{:});
    if ~isRaster, exportgraphics(f,[base '.emf'],'ContentType','vector',eg{:}); end
    exportgraphics(f,[base '.png'],'Resolution',600,eg{:});
else
    u = f.Units; f.Units = 'inches'; p = f.Position; f.Units = u;
    set(f,'Color','w','InvertHardcopy','off','PaperUnits','inches', ...
        'PaperSize',p(3:4),'PaperPositionMode','manual','PaperPosition',[0 0 p(3:4)]);
    if isRaster, ct = '-image'; else, ct = '-vector'; end
    print(f,[base '.pdf'],'-dpdf',ct,'-r600');
    print(f,[base '.svg'],'-dsvg',ct,'-r600');
    if ~isRaster, print(f,[base '.emf'],'-dmeta','-vector'); end
    print(f,[base '.png'],'-dpng','-r600');
end
fprintf('  saved %s.{pdf,svg,%spng}\n', base, repmat('emf,',1,~isRaster));
end

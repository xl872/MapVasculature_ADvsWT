function p = figpos(w_in, h_in, fsCode)
%FIGPOS  figure Position (pixels) for a panel that is placed w_in x h_in
%inches wide/high in the assembled figure. The drawing code of the panel sets
%its own fontsize(fsCode); the figure size is chosen so that every panel
%prints with 7-pt text:
%   figure size = placed size * fsCode/7.
%   figure(Position=figpos(2.25,1.75,11))
k = 7/fsCode;
ppi = get(groot,'ScreenPixelsPerInch');
p = [100 100 round(w_in/k*ppi) round(h_in/k*ppi)];
end

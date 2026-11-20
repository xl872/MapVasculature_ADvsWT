function make_all_panels(outRoot, which)
%MAKE_ALL_PANELS  Regenerate the data panels of the main figures and of
%Fig S4/S5.
%   make_all_panels('<out>')                    % all
%   make_all_panels('<out>', {'Fig3','FigS4'})  % a subset
% Fig 1 has no data panels here; of Fig 2 only panel B is drawn. Each panel is
% written as .pdf/.svg/.emf (vector) and .png (600 dpi) under <out>/FigN/, at
% the size it is placed in the figure (text prints at 7 pt).
if nargin < 2, which = {'Fig2','Fig3','Fig4','Fig5','Fig6','FigS4','FigS5'}; end
addpath(fileparts(mfilename('fullpath')));
if any(strcmp(which,'Fig2')),  plot_Fig2B(fullfile(outRoot,'Fig2'));   end
if any(strcmp(which,'Fig3')),  plot_Fig3(fullfile(outRoot,'Fig3'));    end
if any(strcmp(which,'Fig4')),  plot_Fig4(fullfile(outRoot,'Fig4'));    end
if any(strcmp(which,'Fig5')),  plot_Fig5(fullfile(outRoot,'Fig5'));    end
if any(strcmp(which,'Fig6')),  plot_Fig6(fullfile(outRoot,'Fig6'));    end
if any(strcmp(which,'FigS4')), plot_FigS4(fullfile(outRoot,'FigS4'));  end
if any(strcmp(which,'FigS5')), plot_FigS5(fullfile(outRoot,'FigS5'));  end
end

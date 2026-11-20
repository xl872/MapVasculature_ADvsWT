function S = figstyle()
%FIGSTYLE  Shared settings of the main-figure panels. The plots use the
%colours, plot types and font sizes of the analysis code; panel sizes are set
%with figpos.m.

% voxel size of the grid on which vessel profiles / FWHM were measured:
% the registered volumes (*_regraw+orig, *_sr_reg1+orig) are on the 30 um
% AMBMC template grid (AFNI header DELTA = 0.03 mm).
S.voxel_mm = 0.030;

% statistical maps (Fig 4A, 5, 6A), as plot_FD_slices.m
S.mapAlpha = 0.6;
S.mapUpsample = 3;    % display interpolation of each 30 um map slice
end

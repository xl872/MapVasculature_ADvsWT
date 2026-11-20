# _figures: data panels of the main figures (Fig 2B, 3–6) and Fig S4, S5

Draws the data panels with the drawing code of the analysis scripts (same
colours, plot types, line widths and font sizes). Only the figure size is set,
`figure(Position=figpos(w,h,fs))`, so that every panel prints with 7-pt text
when placed at w × h inches in the 6.5-in-wide figure.

Run with **MATLAB R2023b**. In R2025b every figure counts as a UIFigure, so
`boxplotGroup` falls back to vertical group labels, and `exportgraphics`
behaves differently:

```matlab
make_all_panels('<out>')            % Fig2 ... Fig6, FigS4, FigS5 -> <out>/FigN
make_all_panels('<out>', {'Fig3'})  % only some figures
```

Each panel is written as `.pdf`, `.svg`, `.emf` (vector) and `.png` (600 dpi),
at the full figure size. Map panels (Fig 4A, 5A–E, 6A) are images, so they
have `.png` plus `.pdf`/`.svg` that contain the image.

Fig 1 is not drawn here (Fig 1C: `../Fig1C_Fig3_VesselProfiles/`).

| Script | Panels | Drawing code from | Data |
|---|---|---|---|
| `plot_Fig2B.m` | 2B histogram + thresholds | `newfigures.m` (Fig 2B section) | `Fig2_VesselSegmentation/data/Fig2B_histograms_all_animals.mat` (WT06) |
| `plot_Fig3.m` | 3A/B profiles | `vesselsize.m` | `Fig1C_Fig3_VesselProfiles/data/vsallWT.mat`, `vsallAD.mat` |
| | 3C/D FWHM distributions | `newfigures.m` (`distributionPlot`) | same |
| | 3E WT vs AD | `newfigures.m` (`boxplotGroup`) | same |
| `plot_Fig4.m` | 4A map, colour bar, ROI targets | `plot_FD_slices.m` | `Fig4_Fig5_FigS4_VesselDensityMaps/data/FD_VD_tag1.mat`, `roi_masks/` |
| | 4B box plots (12 ROIs) | `ROIbar.m` | `Fig4B_ROI_values.mat` |
| `plot_Fig5.m` | 5A–E maps + colour bar | `plot_FD_slices.m` | `.../data/FD_VD_tag{-21,-16,-11,-6,-1}.mat` |
| `plot_Fig6.m` | 6A map + colour bar | `plot_FD_slices.m` | `Fig6_FigS5_CapillaryIntensityMaps/data/FD_intensity_tag1.mat` |
| | 6B/C box plots | `ROIbarIntensity.m` | `Fig6BC_ROI_values.mat` |
| `plot_FigS4.m` | S4: control-permutation null distributions of the 12 Fig 4B ROIs | `findROI.m`, via `../_common/plot_control_null.m` | `Fig4_Fig5_FigS4_VesselDensityMaps/data/S4_null_<ROI>.mat` |
| `plot_FigS5.m` | S5: HP, VISp null distributions | `findROI_intensity.m`, via `../_common/plot_control_null.m` | `Fig6_FigS5_CapillaryIntensityMaps/data/S5_null_<ROI>.mat` |

Helpers:

| File | Role |
|---|---|
| `figpos.m` | figure size from the placed size and the code's font size |
| `cbpanel.m` | jet colour bar as its own panel |
| `exportpanel.m` | export at the full figure size |
| `mappanel.m`, `maptile.m`, `mapcrop.m`, `newpanel.m` | map tiles |
| `figstyle.m` | voxel size, map alpha |

`distributionPlot` (BSD, J. Dorn) is in `../_third_party/distributionPlot`.

`../Fig4_Fig5_FigS4_VesselDensityMaps/data/Fig4A_ROI_centroids.mat` is a cache that `plot_Fig4.m` writes: the ROI centroids on their slices, from `roi_masks/`. Delete it to recompute.

Fig S4 panels are 1.625 in wide (4 columns in 6.5 in), in the order ACC, MS, Pir, S2 / CM, DG, Amy-HP, MGN / lEC, mEC, IC, PAG. Fig S5 panels (HP, VISp) are 2.42 in wide, the width at which Fig S5 sits in the manuscript (4.85 in). `plot_FigS4.m` and `plot_FigS5.m` also write `FigS4_values_AD-WT.csv` / `FigS5_values_AD-WT.csv` with the plotted observed value and 95 % line of each ROI. `plot_Fig4.m` writes `Fig4B_stats.csv` (n, means, p and stars of each Fig 4B ROI).

Besides the panels, `plot_Fig4.m` writes `Fig4A_ROI_targets.json` and `plot_Fig5.m` writes `Fig5_tiles.json`. They give the tile positions, and for Fig 4A the arrow target of each ROI label, in inches from the top-left of the panel. Use them to place the ROI labels and region arrows when laying out the figure. The Fig 4A target is the centroid of the ROI mask on its slice (left hemisphere for bilateral ROIs). For DG and PAG (`SIG_TARGET`), whose mask centroids are not on a significant voxel, it is the centroid of the left-hemisphere AD > WT significant voxels of the mask, or the nearest such voxel if the centroid is not on one. The json also gives the number of significant voxels of each mask (`n_sig`, split into AD > WT and AD < WT) and which target rule was used (`target`).

## Conventions

- **Differential maps are AD - WT.** `data/FD_VD_*.mat` and `data/FD_intensity_tag1.mat` store WT - AD, so the scripts flip the sign.
  - VD maps (Fig 4A, 5): red/yellow = higher VD in AD.
  - Residual-signal map (Fig 6A): red/yellow = higher residual signal in AD, that is, lower capillary density in AD. Blue = higher capillary density in AD.
- **Significance mask (Fig 4A, 5, 6A).** `slices.Ah` in `FD_VD_*.mat` / `FD_intensity_tag1.mat` is the Benjamini-Hochberg FDR mask: BH q of the p map inside the brain mask (mask > 0.5), Gaussian-smoothed (sigma = 1 voxel), q < 0.01, then |difference| >= 0.005 (VD) or 0.01 (signal intensity). The panels show the difference where `slices.Ah` is set; voxels outside the brain mask are never marked.
- **FWHM** is converted from voxels of the registered 30 µm AMBMC grid to mm with `S.voxel_mm = 0.030`. Profile x-axes use the same 30 µm spacing.
- **Fig 2B** is the histogram of the image that the segmentation thresholds (`myseg_tissue_top2to50.m`): `exp(-(V-B)/1000+1)` inside the mask, inverted so that vessels are bright, then NLM-filtered. Its counts are those of `histcounts(raw(mask1==1))`.
- **Fig S4, S5.** `S4_null_*.mat` and `S5_null_*.mat` store WT - AD (histogram, Gaussian fit, 95 % line, observed value). `plot_control_null.m` draws them as AD - WT by mirroring the x-axis at plot time (edges and counts reversed and negated; observed value, 95 % line and Gaussian centre negated; x-limits mirrored). Nothing is refitted. The Fig S5 x-axis is the z-scored residual signal intensity difference (AD - WT).
- **Fig 4B stars** come from `ttest2` over the rows of `Fig4B_ROI_values.mat`, as in `ROIbar.m`.

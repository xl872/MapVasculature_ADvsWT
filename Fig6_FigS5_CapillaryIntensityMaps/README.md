# Fig 6 and Fig S5: capillary (residual-signal) maps

Quick plots from `data/` (`addpath ../_common`):
```matlab
plot_FD_slices('data/FD_intensity_tag1.mat')                     % Fig 6A
plot_ROI_boxplots('data/Fig6BC_ROI_values.mat', {'VISp','HP'})   % Fig 6B, 6C
plot_control_null('data/S5_null_HP.mat','HP')                    % Fig S5 (also VISp), AD - WT
```

Full pipeline:
1. `Random_group_result_double_intensity_z.m`: z-scored image, voxels of the top-50 % vessel mask excluded, 450 µm block means, permutation maps.
   Each epoch draws 10 of the 11 animals per group at random; each animal's map and its left-right mirror fill consecutive slots, and the slots are cleared before each group. The intensity is read from each animal's `*_75um.nii_regraw` volume. The seed is fixed (`SEED`) and saved with the draws.
   Significance mask of the maps: Benjamini-Hochberg FDR of the voxel-wise p map (one-sample t-test over the 100 permutation maps), computed over the voxels inside the brain mask (mask > 0.5); the q map is smoothed (`imgaussfilt3`, sigma = 1 voxel) and thresholded at q < 0.01; voxels with |difference| < 0.01 (signal intensity) are then excluded. Voxels outside the brain mask are never marked.
2. `ROIbarIntensity.m`: Fig 6B/C box plots (`data/forROIbarplotINT.mat`, `data/roi_masks/`).
3. `Random_group_result_double_intensity_z_CON.m` + `findROI_intensity.m`: control permutations (100×) and the Fig S5 histograms (`data/S5_null_<ROI>.mat`).
   These values are WT - AD. `../_common/plot_control_null.m` draws them as AD - WT, like Fig 6: it mirrors the x-axis at plot time and leaves the `.mat` files unchanged. The x-axis is labelled "Residual signal intensity difference (AD - WT)": the quantity is the z-scored residual signal intensity.

`data/`:
- `FD_intensity_tag1.mat`: the map slices of Fig 6A: template, difference (WT - AD), significance mask `Ah`, brain mask.
- `Fig6BC_ROI_values.mat`, `.csv`: the Fig 6B/C values of VISp and HP (22 rows per group = 11 animals × 2 hemisphere copies).
- `forROIbarplotINT.mat`: the per-animal block maps that `ROIbarIntensity.m` reads.
- `roi_masks/INTROI<name>.mat`: the HP and VISp masks.
- `S5_null_<ROI>.mat`: the Fig S5 histogram, Gaussian fit, 95 % line and observed value of each ROI.

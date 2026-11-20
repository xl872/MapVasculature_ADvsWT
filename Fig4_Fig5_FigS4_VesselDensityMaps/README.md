# Fig 4, Fig 5, Fig S4: vessel-density (VD) differential maps

Quick plots from `data/` (`addpath ../_common`):
```matlab
plot_FD_slices('data/FD_VD_tag1.mat')                       % Fig 4A (top 50 %)
plot_ROI_boxplots('data/Fig4B_ROI_values.mat', {'Cg','Septal','PIR','SS','CM','DG','Amygdalar','LEn','MGN','MEn','IC','PAG5'})  % Fig 4B
plot_FD_slices('data/FD_VD_tag-21.mat')                     % Fig 5A (2-10 %); -16, -11, -6, -1 = 5B-E
plot_control_null('data/S4_null_Cg.mat','ACC')              % one Fig S4 panel (AD - WT)
```

Full pipeline (needs the per-animal segmentations from `../Fig2_VesselSegmentation`):
1. `Random_group_result_double_density.m`: block VD (450 µm), 100 permutation epochs over random subsets of animals (hemispheres folded), voxel-wise t-test maps. `tag` selects the threshold (1 = top 50 %; −21/−16/−11/−6/−1 = the Fig 5 ranges; 6/11/16/21 = Fig S6).
   Each epoch draws 10 of the 11 animals per group at random; each animal's map and its left-right mirror fill consecutive slots, and the slots are cleared before each group. The seed is fixed (`SEED`) and saved with the draws.
2. `loadandplot.m`: final map figure.
   Significance mask of the maps: Benjamini-Hochberg FDR of the voxel-wise p map (one-sample t-test over the 100 permutation maps), computed over the voxels inside the brain mask (mask > 0.5); the q map is smoothed (`imgaussfilt3`, sigma = 1 voxel) and thresholded at q < 0.01; voxels with |difference| < 0.005 (vessel density) are then excluded. Voxels outside the brain mask are never marked. The same lines are in `Random_group_result_double_density.m` (plot part).
3. `ROIbar.m`: Fig 4B box plots (`data/forROIbarplot.mat`, `data/roi_masks/`).
4. `Random_group_result_double_density_Contral.m` + `findROI.m`: control permutations (5 WT + 5 AD vs 5 WT + 5 AD, 1000×, same sampling) and the Fig S4 histograms (50 bins).
   These values are WT - AD, as `findROI.m` computes them. `../_common/plot_control_null.m` draws them as AD - WT, like Fig 4 and 5: it mirrors the x-axis at plot time (`'Direction','WT-AD'` draws the stored values).

`data/`:
- `FD_VD_tag<tag>.mat`: the map slices of Fig 4A (tag 1) and Fig 5 (tags −21 … −1): template, difference (WT - AD), significance mask `Ah`, brain mask.
- `Fig4B_ROI_values.mat`, `.csv`: the Fig 4B values of the 12 ROIs (`MdenWT`, `MdenAD`: 22 rows per group = 11 animals × 2 hemisphere copies).
- `forROIbarplot.mat`: the per-animal block VD maps that `ROIbar.m` reads.
- `roi_masks/ROI<name>.mat`: the 12 ROI masks.
- `S4_null_<ROI>.mat`: the Fig S4 histogram, Gaussian fit, 95 % line and observed value of each ROI. `S4_observed_VD_tag1.mat`, `.csv`: observed difference (WT - AD) and voxel count of each ROI.
- `Fig4A_ROI_centroids.mat`: ROI centroids, a cache written by `../_figures/plot_Fig4.m`.

ROI names: Cg = ACC, Septal = MS, PIR = Pir, SS = S2, Amygdalar = Amy-HP, LEn = lEC, MEn = mEC, PAG5 = PAG.

# Code and figure data: brain-wide CBV-MRI vascular mapping in APP/PS1 vs WT mice

Each folder holds the scripts for one or more figures and, in `data/`, the processed data that the plotting step reads. Figures that share scripts share a folder. The raw MRI data are on Harvard Dataverse https://doi.org/10.7910/DVN/5YEOZH.

| Folder | Figures |
|---|---|
| `SuperResolution/` | Super-resolution network: inference (used for every SR volume) and training |
| `_preprocessing/` | SR, registration and segmentation of a single animal |
| `Fig1C_Fig3_VesselProfiles/` | Fig 1C (vessel profiles, raw vs SR) and Fig 3 (vessel profiles and size by threshold, WT vs AD) |
| `Fig2_VesselSegmentation/` | Fig 2: AFNI preprocessing, top 2–50 % segmentation, histogram |
| `Fig4_Fig5_FigS4_VesselDensityMaps/` | Fig 4 (VD map at 50 % + ROI box plots), Fig 5 (VD maps by threshold range), Fig S4 (control permutations) |
| `Fig6_FigS5_CapillaryIntensityMaps/` | Fig 6 (residual-signal map + VISp/HP box plots), Fig S5 (control permutations) |
| `FigS1_RawVsSR/` | Fig S1: signal before vs after SR |
| `FigS2_S3_SegmentationCheck/` | Fig S2 and S3: SR segmentation checked against raw data; agreement between SR, raw and interpolated segmentations |
| `FigS6_S7_S8_AtlasComparison/` | Fig S6–S8: regional VD and residual signal vs the Ji et al. 2021 vascular atlas |
| `_figures/` | Data panels of the main figures (Fig 2B, 3–6) and of Fig S4, S5, drawn with the analysis plotting code at print size |
| `_common/` | Shared helpers (block averaging, plotting from `data/`) |
| `_third_party/` | FastNonLocalMeans3D, boxplotGroup, pval_adjust, distributionPlot (with their licenses) |

## Animals and inputs
All 22 animals (11 WT, 11 AD) are analysed from their in vivo 75 µm scans after super-resolution to 37.5 µm (`SuperResolution/`). Each SR volume is registered to the 30 µm AMBMC template and segmented (`_preprocessing/`, `Fig2_VesselSegmentation/`). Harvard Dataverse also holds in vivo 37.5 µm scans of five animals as additional data; the analysis does not use them.

## Group analysis
- **Block maps.** For each threshold, each animal's vessel mask is turned into a map of vessel density in 450 µm blocks (15³ voxels of the 30 µm grid; Fig 4, 5, S6). For Fig 6, the image is z-scored inside the brain, the voxels of the top-50 % vessel mask are excluded and the remaining signal is averaged in the same blocks.
- **Sampling.** In each of 100 permutation epochs, 10 of the 11 animals of each group are drawn at random. Each drawn animal's map and its left-right mirror go into consecutive slots (20 maps per group), and the slot array is cleared before each group. The seed is fixed (`SEED = 20261118`, `rng(SEED,'twister')`), and the draws are saved with the maps.
- **Statistics.** Per epoch, a voxel-wise two-sample t-test (AD vs WT) gives a p map, and the WT - AD difference map is kept where p < 0.05 (for vessel density then smoothed, sigma = 1 voxel). Over the epochs, a one-sample t-test of these maps gives the p map of the group difference.
- **Control permutations (Fig S4, S5).** 5 WT + 5 AD are compared with another 5 WT + 5 AD, with the same block maps and slots: 1000 epochs for vessel density, 100 for the residual signal.
- Scripts: `Random_group_result_double_density.m`, `Random_group_result_double_density_Contral.m` (`Fig4_Fig5_FigS4_VesselDensityMaps/`), `Random_group_result_double_intensity_z.m`, `Random_group_result_double_intensity_z_CON.m` (`Fig6_FigS5_CapillaryIntensityMaps/`). The Fig 3 data are made by `Fig1C_Fig3_VesselProfiles/fig3_vsall.m`.

## Significance of the differential maps
In Fig 4A, 5 and 6A a voxel is coloured where the group difference is significant after Benjamini-Hochberg FDR correction: BH q of the voxel-wise p map within the brain mask (mask > 0.5), Gaussian smoothing of q (sigma = 1 voxel), q < 0.01, and an absolute difference of at least 0.005 (vessel density) or 0.01 (residual signal intensity). The step is in `loadandplot.m`, `Random_group_result_double_density.m` and `Random_group_result_double_intensity_z.m`; the map data in the `data/` folders are exported with it.

## Requirements
- MATLAB (R2023b or later) with the Image Processing, Statistics and Curve Fitting toolboxes. The `_figures` panels are made with R2023b; see `_figures/README.md`.
- AFNI, and its MATLAB reader `BrikLoad` (https://github.com/afni/afni_matlab)
- Python 3 with TensorFlow (SR inference) or PyTorch (SR training)

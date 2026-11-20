# Fig 2: preprocessing and vessel segmentation

- `afni_preprocessing.sh`: `3dUnifize` and two-pass `3dAllineate` (lpa) to the 30 µm AMBMC template, with the per-subject settings stored in the AFNI header history of each registered dataset.
- `myseg_tissue_top2to50.m`: Gaussian background removal, sigmoid normalization, fast non-local means, then masks of the top 50 % … 2 % of voxels in 2 % steps. Writes one `*_resultinrange50.mat` per animal, used by all later analyses.
- `plot_Fig2B_histogram.m`: makes panel B.

`data/Fig2B_histograms_all_animals.mat`: intensity histograms and threshold values for all 22 animals; panel B shows WT06. `data/AD02_threshold_level_map.mat`: threshold-level map of AD02, the candidate subject for the panel C/D overlays.

Needs `../_third_party/FastNonLocalMeans3D` and `BrikLoad`.

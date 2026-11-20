# _common

- `mysizedown.m`, `mysizedownd.m`, `mysizedown_SI.m` — 15³-voxel (450 µm) block averaging used by the group scripts.
- `labels.txt` — Allen label table for `Allen_label_AMBMC30+orig` (depth-first pre-order of the Allen ontology).
- `plot_FD_slices.m`, `plot_ROI_boxplots.m` — draw Fig 4A/5/6A and Fig 4B/6B-C from `data/`.
- `plot_control_null.m` — draws one Fig S4/S5 panel from `data/S4_null_<ROI>.mat` / `S5_null_<ROI>.mat`. The files store WT - AD; the panel is drawn as AD - WT by default (x-axis mirrored at plot time).

External dependencies (not bundled): AFNI MATLAB reader `BrikLoad` / `BrikInfo` (https://github.com/afni/afni_matlab, GPL); MATLAB Image Processing, Statistics and Curve Fitting toolboxes. Bundled third-party code with licenses: `../_third_party/` (FastNonLocalMeans3D — Tristán-Vega et al. 2012, BSD; pval_adjust — Nuno Fachada, MIT; boxplotGroup — Adam Danz, BSD; distributionPlot — J. Dorn, BSD). The Benjamini-Hochberg step of `loadandplot.m`, `Random_group_result_double_density.m` and `Random_group_result_double_intensity_z.m` is written out in the scripts and computes the same q values as `pval_adjust(p,'fdr')`.

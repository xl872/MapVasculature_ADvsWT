# Fig S6, S7, S8: comparison with the Ji et al. 2021 vascular atlas

Put `1-s2.0-S0896627321000805-mmc2.xlsx` (Ji et al. 2021, *Neuron*, supplementary table) in `data/` first.

| Script | Figure | Main data |
|---|---|---|
| `myComparePlotDensity_marked.m` | S6 A–E (tags 1/6/11/16/21 = top 50/40/30/20/10 %) | `Density<tag>_ROI_WTAD_full`, `Density<tag>_ROIsel_WTAD`, `ROIsel_xvalue` |
| `myPlotFigureS6SlopesSTE.m` | S6 F | `Density<tag>WTAD_mice`, `labels.txt`, plus the AMBMC atlas labels and mask (`ATLAS_DIR`: `Allen_label_AMBMC30+orig`, `ambmcmask+orig`) |
| `myComparePlotIntensity_marked.m` | S7 | `IntensityLeft50_SIvalid5_ROI_WTAD_full`, `Intensity1_ROIsel_WTAD`, `ROIsel_xvalue` |
| `myComparePlotIntensity_less50.m` | S8 A–F | `IntensityLeft{10..50}_SIvalid5_ROI_WTAD_full` |

For S6 A, set `SHOW_MARKED = true` in `myComparePlotDensity_marked.m`. For B–E, set it to `false`.

`data/`:
- `Density<tag>_ROI_WTAD`: mean and std over the 100 permutation group maps of every Allen atlas label (`labels.txt`), WT and AD; `Density<tag>_ROI_WTAD_full` adds the parent structures of the Allen hierarchy (voxel-weighted means of their sub-structures).
- `Density<tag>_ROIsel_WTAD`: the same for the 12 hand-drawn ROI masks of Fig 4B (the markers of S6 A). `ROIsel_xvalue`: the atlas-table x value of each marked ROI, averaged over the atlas labels the mask covers.
- `Density<tag>WTAD_mice`: the group maps of each mouse (S6 F).
- `Intensity1_ROI_WTAD_full`, `Intensity1_ROIsel_WTAD`, `IntensityLeft<percent>_SIvalid5_ROI_WTAD_full`: the residual-signal tables of S7 and S8.
- `ROI_voxelcount_AMBMC30.mat`: voxel count of each atlas label (weights of the parent structures).

`upstream/` shows how these tables are made from the per-animal maps; those scripts need the full-size intermediate files:
- `myROIforCompare.m`: `Density<tag>_ROI_WTAD`, `Intensity1_ROI_WTAD`.
- `myAggregateHigherROIs.m`: the `*_full` tables.
- `myROIforCompare_selected.m`, `myROIforCompare_selected_intensity.m`, `myROIxLengthDensity_selected.m`: `*_ROIsel_WTAD`, `ROIsel_xvalue`.
- `myseg_tissue_less50.m`, `Random_group_result_double_intensity_z_less50.m`, `findROI_intensity_less50.m`: the residual 10–50 % masks, maps and `IntensityLeft*` tables.

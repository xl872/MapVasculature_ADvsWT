# Fig 1C and Fig 3: vessel signal profiles

| Script | Makes |
|---|---|
| `fig3_vsall.m` | `data/vsall<WT|AD>.mat` from the segmentations of all 22 animals (`../Fig2_VesselSegmentation`): slice `sn = 413`, all 25 thresholds, all 11 animals of each group |
| `vesselsize.m`, section `vessExp` | Fig 1C profile: mean normalized signal across vessels, 75 µm raw vs 37.5 µm SR (`data/before.mat`, `data/after.mat`) |
| `beforeafter_vesselsize.m` | Fig 1C images: raw and SR slices with the zoom box; extracts the raw/SR profiles of the vessels inside it |
| `vesselsize.m`, sections `Vesselpercent...` | Fig 3A/B: mean profile per threshold (50/30/20/10/2 %), WT and AD (`data/vsallWT.mat`, `data/vsallAD.mat`) |
| `newfigures.m`, `vesselsize.m` (last section) | Fig 3C/D: violin plots of vessel size per threshold (`distributionPlot`) |
| `newfigures.m`, section `Percentbox` | Fig 3E: WT vs AD vessel size per threshold |

`vsall<WT|AD>.mat`: `VC` (threshold × animal × 15-point profile, in voxels of the registered 30 µm grid), `VW` (threshold × animal × vessel: width at half prominence, in voxels), `numPixels` (threshold × animal: pixel count of each vessel cross-section). `data/before.mat`, `data/after.mat`: the raw and SR vessel profiles of Fig 1C. `data/AD_ROI.mat`, `data/WT_ROI.mat`: mean signal of the registered SR (`af_ROI`) and raw (`bf_ROI`) volumes in four brain divisions, read by the first section of `newfigures.m`.

`fig3_vsall.m` makes `data/vsall*.mat`; the other scripts draw from `data/`.

Needs `../_third_party/distributionPlot` and `../_third_party/boxplotGroup`.

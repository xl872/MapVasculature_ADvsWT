# Fig S2 and S3: segmentation checks (AD02)

- `Histogram_intensity_new.m`: Fig S2 A–D, and the Fig S3 panels for the three masks (SR, raw, non-learned interpolation).
- `Histogram_intensity_final.m`: the Fig S3 table (Dice, Jaccard, three-way consensus, CV of VD; AD02–AD05).
- `myseg_tissue_raw.m`, `myreg_regonly_regNN.sh`: segmentation and registration of the raw and interpolated inputs.

`data/`:
- `seg_check_result.mat`: masked z-values, line profiles and summaries.
- `seg_agreement_metrics.mat`: the S3 table.
- `lines_AZ_Vasculature_AD02.mat`: the 10 hand-drawn lines.
- `AD02_S2S3_arrays.mat`: the image and mask slices of the panels.

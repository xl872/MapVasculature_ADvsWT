# Fig S1: signal intensity before vs after super-resolution

- `FigS1_plot_from_data.m`: draws all panels from `data/`.
- `FigS1_raw_vs_SR.m`: computation. It needs the raw 75 µm `*.nii` and their `*.sr.nii.gz` in `DATA_DIR/AD` and `DATA_DIR/WT`. Each volume is masked at intensity > 1000 and averaged into 2×2×2 blocks (`imresize3`, scale 0.01). It fits SR = a·raw + b per animal, then compares a and b with a two-sample t-test.

`data/FigS1_block_means.csv`: 8 block means per animal before and after SR. `data/FigS1_linear_fit.csv`: a and b per animal. The t-test on these gives p = 0.125 (a) and p = 0.107 (b). Animals: AD01–AD11 incl. AD07_1/AD07_2 (12), and WT01–09, WT11, WT12 (11).

# _preprocessing: SR, registration and segmentation of one animal

These scripts turn one animal's in vivo 75 µm scan into the inputs of the group scripts: the SR volume, the registered volumes and the vessel masks.

Order, per animal:

1. **Reorient (WT10 only):** `python reorient_WT10.py <WT10_75um.nii> <WT10_75um_reoriented.nii>`.
   The SR input needs a diagonal affine. The header of WT10's 75 µm scan is oblique (16.8°) with axis codes L,I,P, while the other scans are R,I,A. The script turns the volume to R,I,A: array axes 0 and 2 are flipped (no permutation, no resampling, values unchanged), and the affine is set to orthogonal 0.075 mm, keeping the world position of the array centre.
2. **SR:** `python run_sr_one.py <in_75um.nii> <out.sr.nii.gz> <out.linear.nii.gz>`.
   The computation of `../SuperResolution/inference_XinSR/run_sr_inference.py` with the weights `XinSR.h5`; the input and output paths come from the command line. Needs TensorFlow (tested with Python 3.9, TensorFlow 2.9.1, nibabel 5.0.0).
3. **Registration (AFNI):** `bash reg_one.sh <ID> <sr.nii.gz> <raw_75um.nii>` (set `VESS`).
   The SR volume is registered to the WT06 SR base (`3dUnifize -clfrac 0.5`, `3dAllineate -lpa -twopass -cmass`, 30 µm grid). The raw 75 µm volume is then registered to the animal's SR volume (`*_75um.nii_regraw`, used for the residual signal intensity).
4. **Segmentation:** `seg_one('<...>_reg1+orig.HEAD', outDir, 'top2to50', 1000)` for the vessel masks (`*_resultinrange50.mat`), and `seg_one(..., 'less50', 500)` (AD) or `1000` (WT) for the residual 10–50 % masks of Fig S7/S8. The computation is that of `../Fig2_VesselSegmentation/myseg_tissue_top2to50.m` and `../FigS6_S7_S8_AtlasComparison/upstream/myseg_tissue_less50.m`.

The group scripts then read `AD/` and `WT/` in the work folder (see the figure folders).
